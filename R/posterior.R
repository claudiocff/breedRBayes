# ---------------------------------------------------------------------------
# Posterior extraction: read BGLR's on-disk MCMC output into coda objects.
# ---------------------------------------------------------------------------

#' Names of the random (BRR) terms in a fit, matching their on-disk file keys
#' @keywords internal
.random_keys <- function(fit) {
  names(fit$meta)[vapply(fit$meta, function(m) identical(m$model, "BRR"), logical(1))]
}

#' Read the per-chain variance-component traces from BGLR output
#'
#' @return A list (one element per chain) of matrices `[nSamples x nVC]` whose
#'   columns are the random-term variances plus `varE`.
#' @keywords internal
.read_varchains <- function(fit) {
  keys <- .random_keys(fit)
  lapply(fit$paths, function(prefix) {
    cols <- lapply(keys, function(k) scan(paste0(prefix, "ETA_", k, "_varB.dat"),
                                          quiet = TRUE))
    varE <- scan(paste0(prefix, "varE.dat"), quiet = TRUE)
    m <- do.call(cbind, c(cols, list(varE)))
    colnames(m) <- c(keys, "varE")
    m
  })
}

#' Convert a fitted model's posterior draws to a `coda` object
#'
#' Reads the variance-component chains (and optionally the intercept `mu`) written
#' by BGLR and returns them as a [coda::mcmc.list] (one component per chain), ready
#' for [mcmc_diag()] or `coda`/`bayesplot` functions.
#'
#' @param x A `breedRB_fit`.
#' @param what Which quantities to return: any of `"varcomp"` and `"mu"`.
#' @param ... Unused.
#' @return A [coda::mcmc.list].
#' @export
as_mcmc <- function(x, ...) UseMethod("as_mcmc")

#' @rdname as_mcmc
#' @export
as_mcmc.breedRB_fit <- function(x, what = "varcomp", ...) {
  vc <- .read_varchains(x)
  chains <- lapply(seq_along(vc), function(i) {
    m <- vc[[i]]
    if ("mu" %in% what) {
      mu <- scan(paste0(x$paths[i], "mu.dat"), quiet = TRUE)
      m <- cbind(m, mu = mu[seq_len(nrow(m))])
    }
    if (!"varcomp" %in% what) m <- m[, setdiff(colnames(m), c(colnames(vc[[i]]))), drop = FALSE]
    coda::mcmc(m, thin = x$control$thin)
  })
  coda::mcmc.list(chains)
}
