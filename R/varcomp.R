# ---------------------------------------------------------------------------
# Variance components: posterior summaries and full posterior draws.
# ---------------------------------------------------------------------------

#' Posterior summary of a numeric vector of draws
#' @keywords internal
.post_summary <- function(z, prob = 0.95) {
  a <- (1 - prob) / 2
  data.frame(
    mean   = mean(z),
    median = median(z),
    sd     = sd(z),
    lower  = unname(quantile(z, a)),
    upper  = unname(quantile(z, 1 - a))
  )
}

#' Posterior variance components of a fitted model
#'
#' Returns the posterior distribution and summaries of every random-term variance
#' and the residual variance, pooling all chains. For a **random-regression**
#' model the across-genotype intercept--slope (co)variance is appended as extra
#' `cov(...)` rows: these are the realised coefficient covariances estimated from
#' the posterior draws (see [rr_gradient()]), reported even though the intercept
#' and slope terms are fitted with independent priors.
#'
#' @param fit A `breedRB_fit`.
#' @param prob Central credible-interval mass (default 0.95).
#' @param draws Logical; if `TRUE` also return the pooled posterior draws matrix.
#' @return A data frame of per-component summaries (`term`, `mean`, `median`,
#'   `sd`, `lower`, `upper`) — the random-term variances, `varE`, and, for a
#'   random regression, the reaction-norm `cov(intercept, slope)` rows — with the
#'   pooled variance draws attached as attribute `"draws"` when `draws = TRUE`
#'   (the draws matrix covers the variance components only).
#' @examples
#' \donttest{ varcomp(fit) }
#' @export
varcomp <- function(fit, prob = 0.95, draws = FALSE) {
  stopifnot(inherits(fit, "breedRB_fit"))
  if (isTRUE(fit$response$multitrait)) return(.varcomp_mt(fit, prob))
  vc <- .read_varchains(fit)
  pooled <- do.call(rbind, vc)
  out <- do.call(rbind, lapply(colnames(pooled), function(nm) {
    cbind(data.frame(term = nm), .post_summary(pooled[, nm], prob))
  }))
  rownames(out) <- NULL
  # Reaction-norm intercept-slope (co)variance(s), if any (appended after varE).
  cov_rows <- .rr_cov_summaries(fit, prob)
  if (!is.null(cov_rows)) {
    out <- rbind(out, cov_rows)
    rownames(out) <- NULL
  }
  if (draws) attr(out, "draws") <- pooled
  out
}
