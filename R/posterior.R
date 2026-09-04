# ---------------------------------------------------------------------------
# Posterior extraction: read BGLR's on-disk MCMC output into coda objects.
# ---------------------------------------------------------------------------

#' Names of the random (BRR) *logical terms* in a fit (meta keys)
#'
#' A split random-regression term counts once here (its logical key); use
#' [.varB_keys()] for the per-block on-disk variance-file keys.
#' @keywords internal
.random_keys <- function(fit) {
  names(fit$meta)[vapply(fit$meta, function(m) identical(m$model, "BRR"), logical(1))]
}

#' On-disk variance-file keys of the random terms (one per BGLR ETA block)
#'
#' Expands each logical random term to its `eta_keys`, so a per-degree-split
#' random regression contributes one key per Legendre degree (each with its own
#' `varB` trace). Non-split terms map to themselves.
#' @keywords internal
.varB_keys <- function(fit) {
  unlist(lapply(.random_keys(fit), function(k) fit$meta[[k]]$eta_keys %||% k),
         use.names = FALSE)
}

#' Read the per-chain variance-component traces from BGLR output
#'
#' @return A list (one element per chain) of matrices `[nSamples x nVC]` whose
#'   columns are the random-term variances (one per BGLR ETA block) plus `varE`.
#' @keywords internal
.read_varchains <- function(fit) {
  keys <- .varB_keys(fit)
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
