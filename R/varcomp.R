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
#' and the residual variance, pooling all chains.
#'
#' @param fit A `breedRB_fit`.
#' @param prob Central credible-interval mass (default 0.95).
#' @param draws Logical; if `TRUE` also return the pooled posterior draws matrix.
#' @return A data frame of per-component summaries (`term`, `mean`, `median`,
#'   `sd`, `lower`, `upper`), with the pooled draws attached as attribute
#'   `"draws"` when `draws = TRUE`.
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
  if (draws) attr(out, "draws") <- pooled
  out
}
