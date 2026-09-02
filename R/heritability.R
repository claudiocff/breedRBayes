# ---------------------------------------------------------------------------
# Heritability as a full posterior distribution (computed per MCMC draw).
# ---------------------------------------------------------------------------

#' Map a user-supplied term identifier (label or key) to its on-disk key
#' @keywords internal
.resolve_term <- function(fit, id) {
  keys   <- names(fit$meta)
  labels <- vapply(fit$meta, `[[`, character(1), "label")
  if (id %in% keys)   return(id)
  hit <- keys[labels == id]
  if (length(hit)) return(hit[1])
  stop("Term '", id, "' not found. Available: ",
       paste(labels, collapse = ", "), call. = FALSE)
}

#' Keys of the random terms that carry a genomic (`vm`) component
#' @keywords internal
.vm_keys <- function(fit) {
  names(fit$meta)[vapply(fit$meta, function(m)
    any(vapply(m$components, function(c) identical(c$kind, "vm"), logical(1))),
    logical(1))]
}

#' Posterior distribution of heritability
#'
#' Computes narrow-sense (genomic) heritability for every MCMC draw, returning the
#' full posterior distribution rather than a single point estimate:
#' \deqn{h^2 = \sigma^2_{genetic} / (\sigma^2_{genetic} + \sum \sigma^2_{other} + \sigma^2_e).}
#'
#' @param fit A `breedRB_fit`.
#' @param genetic Character vector of genetic term identifiers (label or key)
#'   forming the numerator. Defaults to the `vm()` term(s).
#' @param denominator Character vector of term identifiers summed in the
#'   denominator alongside `varE`. Defaults to all random terms.
#' @param prob Central credible-interval mass (default 0.95).
#' @return A list of class `breedRB_h2`: `summary` (mean/median/sd/CI) and
#'   `draws` (posterior vector of \eqn{h^2}).
#' @examples
#' \donttest{ h2 <- heritability(fit); h2$summary; hist(h2$draws) }
#' @export
heritability <- function(fit, genetic = NULL, denominator = NULL, prob = 0.95) {
  stopifnot(inherits(fit, "breedRB_fit"))
  if (isTRUE(fit$response$multitrait)) return(.heritability_mt(fit, prob))
  if (is.null(genetic)) {
    genetic <- .vm_keys(fit)
    if (!length(genetic)) stop("No vm() genetic term found; pass `genetic=` explicitly.",
                               call. = FALSE)
  }
  gkeys <- vapply(genetic, .resolve_term, character(1), fit = fit)
  vc <- .read_varchains(fit)
  pooled <- do.call(rbind, vc)

  rand_cols <- setdiff(colnames(pooled), "varE")
  dkeys <- if (is.null(denominator)) rand_cols else
    vapply(denominator, .resolve_term, character(1), fit = fit)

  num   <- rowSums(pooled[, gkeys, drop = FALSE])
  denom <- rowSums(pooled[, union(dkeys, gkeys), drop = FALSE]) + pooled[, "varE"]
  h2    <- num / denom

  structure(list(summary = cbind(data.frame(quantity = "h2"), .post_summary(h2, prob)),
                 draws = h2),
            class = "breedRB_h2")
}

#' @export
print.breedRB_h2 <- function(x, ...) {
  cat("Posterior heritability\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}
