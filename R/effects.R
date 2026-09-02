# ---------------------------------------------------------------------------
# Genetic values (GEBV): back-map BGLR's ridge coefficients (in PC space) to the
# original genotype scale, with posterior summaries.
# ---------------------------------------------------------------------------

#' Posterior genomic estimated breeding values (GEBV)
#'
#' Reconstructs per-level genetic values for a `vm()` term from the saved BGLR
#' effect samples, mapping the ridge coefficients back through the PC rotation
#' (`genetic value = PC %*% b`), pooled across chains.
#'
#' @param fit A `breedRB_fit` (single-trait).
#' @param term Genetic term identifier (label or key). Defaults to the first
#'   `vm()` term.
#' @param prob Central credible-interval mass (default 0.95).
#' @return A data frame with `ID`, posterior `gebv` (mean), `sd`, `lower`, `upper`,
#'   ordered by decreasing `gebv`, with the pooled draws matrix
#'   (`nDraws x nLevel`) attached as attribute `"draws"`.
#' @examples
#' \donttest{ head(gebv(fit)) }
#' @export
gebv <- function(fit, term = NULL, prob = 0.95) {
  stopifnot(inherits(fit, "breedRB_fit"))
  if (isTRUE(fit$response$multitrait)) {
    stop("gebv() currently supports single-trait fits; use the per-trait chains directly.",
         call. = FALSE)
  }
  if (is.null(term)) {
    vk <- .vm_keys(fit)
    if (!length(vk)) stop("No vm() term found; pass `term=` explicitly.", call. = FALSE)
    term <- vk[1]
  }
  key  <- .resolve_term(fit, term)
  meta <- fit$meta[[key]]
  vmc  <- Filter(function(c) identical(c$kind, "vm"), meta$components)
  if (!length(vmc)) stop("Term '", term, "' has no vm() genomic component.", call. = FALSE)
  PC   <- vmc[[1]]$pc
  levs <- rownames(PC)

  draws <- do.call(rbind, lapply(fit$paths, function(prefix) {
    b <- BGLR::readBinMat(paste0(prefix, "ETA_", key, "_b.bin"))
    b %*% t(PC)                                   # [nSamples x nLevel]
  }))
  colnames(draws) <- levs

  a <- (1 - prob) / 2
  out <- data.frame(
    ID     = levs,
    gebv   = colMeans(draws),
    sd     = apply(draws, 2, sd),
    lower  = apply(draws, 2, quantile, probs = a),
    upper  = apply(draws, 2, quantile, probs = 1 - a),
    row.names = NULL
  )
  out <- out[order(-out$gebv), ]
  attr(out, "draws") <- draws
  out
}
