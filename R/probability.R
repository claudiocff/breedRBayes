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
