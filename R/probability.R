# ---------------------------------------------------------------------------
# Posterior probability of selection: how often each level of a term ranks in
# the top (or bottom) fraction across the MCMC draws.
# ---------------------------------------------------------------------------

#' Posterior probability of ranking in the top fraction
#'
#' For a model term, compute the posterior probability that each level ranks
#' among the best `threshold` fraction of levels. Working from the pooled
#' posterior draws (see [solution()]), each MCMC draw is ranked independently and
#' the top \eqn{k = round(threshold \times n)} levels are flagged; the reported
#' probability is the fraction of draws in which a level is flagged. This is the
#' Bayesian "probability of being in the top X\%" used to rank selection
#' candidates by their whole posterior, not just their point estimate.
#'
#' Because every draw is ranked on its own, adding the intercept `mu` would not
#' change the ranking, so `pr()` is unaffected by centering.
#'
#' @param fit A `breedRB_fit` (single-trait).
#' @param term Term identifier (label as written in the formula, e.g. `"gen"`, or
#'   the internal key). If `NULL`, the first random term is used.
#' @param type Optional; one of `"random"` or `"fixed"`, validated against the
#'   term's actual role. For a fixed term fitted with treatment contrasts the
#'   reference level is not among the coefficients and is therefore excluded.
#' @param threshold Selected fraction in `(0, 1)`; `0.20` = top 20%. Ignored when
#'   `pair = TRUE`.
#' @param higher Logical (default `TRUE`). If `TRUE` the "top" is the highest
#'   values (e.g. yield); set `FALSE` when smaller values are better. Ignored when
#'   `pair = TRUE`.
#' @param pair Logical (default `FALSE`). If `TRUE`, return the pairwise table of
#'   posterior probabilities \eqn{P(A > B)} for every pair of levels instead of
#'   the top-fraction summary.
#' @return If `pair = FALSE`, a data frame with `effect` (level name), `solution`
#'   (posterior mean), and `prob` (posterior probability of ranking in the
#'   selected fraction), ordered by decreasing `prob`, with `k` and `threshold`
#'   attached as attributes. If `pair = TRUE`, a data frame with `A`, `B`, and
#'   `prob` \eqn{= P(A > B)} for every unordered pair, where `A` is the member
#'   with the larger posterior mean (so `prob >= 0.5` for well-separated pairs),
#'   ordered by decreasing `prob`.
#' @examples
#' \donttest{
#' # probability each genotype is in the top 20%
#' pr(fit, term = "gen", type = "random", threshold = 0.20)
#' # pairwise probabilities P(A > B)
#' pr(fit, term = "gen", type = "random", pair = TRUE)
#' }
#' @seealso [solution()] for the underlying posterior draws.
#' @export
pr <- function(fit, term = NULL, type = NULL, threshold = 0.20, higher = TRUE,
               pair = FALSE) {
  stopifnot(inherits(fit, "breedRB_fit"))
  if (!isTRUE(pair) && (!is.numeric(threshold) || length(threshold) != 1L ||
                        threshold <= 0 || threshold >= 1)) {
    stop("`threshold` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }

  sol   <- solution(fit, term = term, type = type)
  draws <- attr(sol, "draws")                       # [nDraws x nEffect]
  n     <- ncol(draws)

  if (isTRUE(pair)) {
    means <- colMeans(draws)
    ord   <- order(-means)                          # A = higher posterior mean
    D     <- draws[, ord, drop = FALSE]
    nm    <- colnames(D)
    # P(A > B): compare each column against all others in one sweep per column.
    P <- matrix(NA_real_, n, n, dimnames = list(nm, nm))
    for (i in seq_len(n)) P[i, ] <- colMeans(D[, i] > D)
    ij  <- utils::combn(n, 2L)                       # unordered pairs, A-index < B-index
    out <- data.frame(A    = nm[ij[1L, ]],
                      B    = nm[ij[2L, ]],
                      prob = P[cbind(ij[1L, ], ij[2L, ])],
                      row.names = NULL)
    out <- out[order(-out$prob), ]
    rownames(out) <- NULL
    attr(out, "comparison") <- "P(A > B)"
    return(out)
  }

  k     <- max(1L, round(threshold * n))            # levels flagged per draw

  # Rank each draw independently; flag the top-k (or bottom-k) levels.
  sgn  <- if (isTRUE(higher)) -1 else 1
  flag <- apply(draws, 1L, function(r) rank(sgn * r, ties.method = "min") <= k)
  probs <- rowMeans(flag)                            # names carried from column names

  out <- data.frame(
    effect   = names(probs),
    solution = colMeans(draws)[names(probs)],
    prob     = as.numeric(probs),
    row.names = NULL
  )
  out <- out[order(-out$prob, -out$solution), ]
  attr(out, "k")         <- k
  attr(out, "threshold") <- threshold
  out
}

#' Posterior probability that new genotypes exceed a training-population threshold
#'
#' Scores genotypes that were **not** in the training data (from their markers,
#' via [predict.breedRB_fit()]) and reports, for each, the posterior probability
#' that its genomic value clears a bar defined by the **training population** —
#' e.g. "probability this new line beats the top 10% of the training set".
#'
#' The calculation is fully Bayesian and done per posterior draw: in each MCMC
#' draw the bar is the corresponding quantile of the *training* genotypes'
#' breeding values in that draw (so the threshold itself carries uncertainty),
#' and the new genotype's predicted value in that same draw is compared against
#' it. The reported probability is the fraction of draws the genotype clears the
#' bar. Because the whole posterior is used, a new genotype with **more marker
#' information** (a tighter predictive posterior) that sits above the bar earns a
#' higher probability than an equally-ranked but more uncertain one — two lines
#' with the same point prediction can have very different probabilities.
#'
#' Requires a marker-based `mrk()` term (a `vm()` fit cannot score new
#' genotypes). The probability is invariant to `add_mu`: the intercept shifts the
#' new value and the training bar equally, so it cancels in the comparison.
#'
#' @param fit A `breedRB_fit` (single-trait) with a `mrk()` genomic term.
#' @param M_new Marker matrix for the genotypes to score: genotype IDs in the row
#'   names, markers in columns (same requirements as [predict.breedRB_fit()];
#'   missing calls are imputed with the training mean).
#' @param threshold Training-population fraction in `(0, 1)` defining the bar;
#'   `0.10` = the top 10% of the training population (its 90th percentile when
#'   `higher = TRUE`).
#' @param term Genomic term identifier (label or key). Defaults to the first
#'   genomic term.
#' @param higher Logical (default `TRUE`). If `TRUE` the bar is the upper
#'   `threshold` tail (larger is better, e.g. yield); set `FALSE` when smaller
#'   values are better (then it is the lower tail).
#' @param add_mu Logical (default `FALSE`). Whether the reported `prediction`
#'   column is on the overall-mean scale (`+ mu`) or the deviation scale. Does not
#'   affect `prob` (the intercept cancels in the comparison).
#' @param prob Central credible-interval mass for the underlying prediction
#'   (default 0.95); affects only the internal prediction, not `prob`.
#' @return A data frame with `ID`, `prediction` (posterior mean predicted value),
#'   `sd` (posterior SD of the prediction — the genotype's predictive
#'   uncertainty), and `prob` (posterior probability of clearing the
#'   training-population bar), ordered by decreasing `prob`. The pooled
#'   new-genotype draws are attached as attribute `"draws"`, and the
#'   posterior-mean bar as attribute `"cutoff"`.
#' @examples
#' \donttest{
#' fit <- bbglr(yield ~ 1, random = ~ mrk(gen, M), data = dat, relmat = list(M = M))
#' # P(each new line beats the top 10% of the training population)
#' predict_pr(fit, M_new, threshold = 0.10)
#' }
#' @seealso [predict.breedRB_fit()] for the point predictions, [pr()] for the
#'   in-sample analogue.
#' @export
predict_pr <- function(fit, M_new, threshold = 0.10, term = NULL,
                       higher = TRUE, add_mu = FALSE, prob = 0.95) {
  stopifnot(inherits(fit, "breedRB_fit"))
  if (isTRUE(fit$response$multitrait)) {
    stop("predict_pr() currently supports single-trait fits.", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      threshold <= 0 || threshold >= 1) {
    stop("`threshold` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  if (is.null(term)) {
    vk <- .vm_keys(fit)
    if (!length(vk)) {
      stop("predict_pr() needs a genomic mrk() term to score new genotypes.",
           call. = FALSE)
    }
    term <- vk[1]
  }
  key <- .resolve_term(fit, term)

  # new-genotype prediction draws and the training-population draws, from the
  # same term and same intercept handling -> aligned draw-by-draw.
  pred   <- predict(fit, M_new, term = key, add_mu = add_mu, prob = prob)
  newd   <- attr(pred, "draws")                         # [nDraws x nNew]
  strain <- solution(fit, term = key, type = "random", add_mu = add_mu)
  traind <- attr(strain, "draws")                       # [nDraws x nTrain]
  if (nrow(newd) != nrow(traind)) {
    stop("Internal error: prediction and training draws are not aligned.",
         call. = FALSE)
  }

  # Per-draw bar = quantile of the TRAINING population within that draw.
  q      <- if (isTRUE(higher)) 1 - threshold else threshold
  cutoff <- apply(traind, 1L, stats::quantile, probs = q)   # length nDraws
  flags  <- if (isTRUE(higher)) newd > cutoff else newd < cutoff  # cutoff recycles per row
  probs  <- colMeans(flags)

  out <- data.frame(
    ID         = colnames(newd),
    prediction = colMeans(newd),
    sd         = apply(newd, 2L, stats::sd),
    prob       = as.numeric(probs),
    row.names  = NULL
  )
  out <- out[order(-out$prob, -out$prediction), ]
  attr(out, "threshold") <- threshold
  attr(out, "higher")    <- higher
  attr(out, "cutoff")    <- mean(cutoff)                # posterior-mean bar, for reference
  attr(out, "draws")     <- newd
  out
}
