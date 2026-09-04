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
#' model the per-coefficient variances already appear as first-class variance
#' components — with the per-degree split each reaction-norm coefficient (the
#' intercept and every Legendre degree) is fitted as its own \pkg{BGLR} kernel
#' with its own variance component (e.g. `Entry`, `Entry_leg_gradient__3__deg1`,
#' ...). Only the **off-diagonal** entries of the across-genotype coefficient
#' (co)variance matrix \eqn{K} are appended after `varE`: a `cov(...)` row for
#' each coefficient pair (intercept--slope and slope--slope). These are the
#' realised coefficient covariances estimated from the posterior draws (see
#' [rr_gradient()]), reported even though the coefficient terms are fitted with
#' independent priors. The variance diagonal is not duplicated as `var(...)`
#' rows; the full \eqn{K} for the across-gradient genetic-correlation surface is
#' assembled by [rr_gradient()].
#'
#' @param fit A `breedRB_fit`.
#' @param prob Central credible-interval mass (default 0.95).
#' @param draws Logical; if `TRUE` also return the pooled posterior draws matrix.
#' @return A data frame of per-component summaries (`term`, `mean`, `median`,
#'   `sd`, `lower`, `upper`) — the random-term variance components, `varE`, and,
#'   for a random regression, the reaction-norm coefficient `cov(...)` rows (the
#'   off-diagonal entries of \eqn{K}) — with the pooled variance draws
#'   attached as attribute `"draws"` when `draws = TRUE` (the draws matrix covers
#'   the BGLR variance components only).
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
  # Reaction-norm coefficient covariances (off-diagonal of K), if any (after varE).
  cov_rows <- .rr_cov_summaries(fit, prob)
  if (!is.null(cov_rows)) {
    out <- rbind(out, cov_rows)
    rownames(out) <- NULL
  }
  if (draws) attr(out, "draws") <- pooled
  out
}
