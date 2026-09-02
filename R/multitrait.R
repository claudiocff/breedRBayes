# ---------------------------------------------------------------------------
# Multi-trait post-processing: reconstruct per-draw genetic and residual
# covariance matrices from BGLR::Multitrait output, and derive per-trait
# heritabilities and genetic correlations.
# ---------------------------------------------------------------------------

#' Column (i, j) index pairs of BGLR's covariance trace layout
#'
#' BGLR writes the upper triangle of the (symmetric) covariance matrix in
#' row-major order: (1,1), (1,2), ..., (1,t), (2,2), ..., (t,t).
#' @keywords internal
.vech_index <- function(t) {
  idx <- list()
  for (i in seq_len(t)) for (j in i:t) idx[[length(idx) + 1]] <- c(i, j)
  idx
}

#' Reconstruct per-draw covariance matrices from a `vech`-encoded trace file
#' @return An array `[nDraws x t x t]`.
#' @keywords internal
.read_cov_file <- function(file, t) {
  V <- as.matrix(read.table(file))
  idx <- .vech_index(t)
  n <- nrow(V)
  out <- array(0, dim = c(n, t, t))
  for (c in seq_along(idx)) {
    i <- idx[[c]][1]; j <- idx[[c]][2]
    out[, i, j] <- V[, c]; out[, j, i] <- V[, c]
  }
  out
}

#' ETA index (position in the BGLR ETA list) of a term key
#' @keywords internal
.eta_index <- function(fit, key) match(key, names(fit$meta))

#' Is the term identified by `key` a factor-analytic term?
#' @keywords internal
.is_fa_term <- function(fit, key) {
  kinds <- vapply(fit$meta[[key]]$components, `[[`, character(1), "kind")
  any(kinds == "fa")
}

#' Reconstruct per-draw FA genetic covariance from loadings and specific variances
#'
#' BGLR writes the loadings `W` (`[nDraws x (t*nF)]`) and specific variances
#' `PSI` (`[nDraws x t]`) rather than a `vech`-encoded `Omega`. The genetic
#' covariance for each draw is `Omega = W W' + diag(PSI)`.
#' @return An array `[nDraws x t x t]`.
#' @keywords internal
.read_fa_cov <- function(prefix, idx, t) {
  W   <- as.matrix(read.table(paste0(prefix, "W_", idx, ".dat")))
  PSI <- as.matrix(read.table(paste0(prefix, "PSI_", idx, ".dat")))
  n   <- nrow(W)
  nF  <- ncol(W) / t
  out <- array(0, dim = c(n, t, t))
  for (d in seq_len(n)) {
    Wd <- matrix(W[d, ], nrow = t, ncol = nF)   # column-major: trait varies fastest
    out[d, , ] <- tcrossprod(Wd) + diag(PSI[d, ], t)
  }
  out
}

#' Per-draw genetic covariance matrices for a multi-trait term
#' @keywords internal
.mt_gen_cov <- function(fit, key, chain = 1L) {
  t <- length(fit$response$traits)
  idx <- .eta_index(fit, key)
  if (.is_fa_term(fit, key)) {
    .read_fa_cov(fit$paths[chain], idx, t)
  } else {
    .read_cov_file(paste0(fit$paths[chain], "Omega_", idx, ".dat"), t)
  }
}

#' Per-draw residual covariance matrices (multi-trait)
#' @keywords internal
.mt_res_cov <- function(fit, chain = 1L) {
  t <- length(fit$response$traits)
  .read_cov_file(paste0(fit$paths[chain], "R.dat"), t)
}

#' Multi-trait variance components (pooled over chains)
#' @keywords internal
.varcomp_mt <- function(fit, prob) {
  traits <- fit$response$traits
  gkeys  <- .vm_keys(fit); if (!length(gkeys)) gkeys <- .random_keys(fit)
  gkey   <- gkeys[1]

  G <- do.call(abind_draws, lapply(seq_along(fit$paths), .mt_gen_cov, fit = fit, key = gkey))
  R <- do.call(abind_draws, lapply(seq_along(fit$paths), .mt_res_cov, fit = fit))

  rows <- list()
  for (i in seq_along(traits)) {
    rows[[length(rows) + 1]] <- cbind(
      data.frame(term = paste0("genetic:", traits[i])), .post_summary(G[, i, i], prob))
    rows[[length(rows) + 1]] <- cbind(
      data.frame(term = paste0("residual:", traits[i])), .post_summary(R[, i, i], prob))
  }
  if (length(traits) >= 2) {
    for (i in 1:(length(traits) - 1)) for (j in (i + 1):length(traits)) {
      rg <- G[, i, j] / sqrt(G[, i, i] * G[, j, j])
      rows[[length(rows) + 1]] <- cbind(
        data.frame(term = paste0("rg:", traits[i], "-", traits[j])), .post_summary(rg, prob))
    }
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

#' Multi-trait per-trait heritability (pooled over chains)
#' @keywords internal
.heritability_mt <- function(fit, prob) {
  traits <- fit$response$traits
  gkeys  <- .vm_keys(fit); if (!length(gkeys)) gkeys <- .random_keys(fit)
  gkey   <- gkeys[1]
  G <- do.call(abind_draws, lapply(seq_along(fit$paths), .mt_gen_cov, fit = fit, key = gkey))
  R <- do.call(abind_draws, lapply(seq_along(fit$paths), .mt_res_cov, fit = fit))
  draws <- sapply(seq_along(traits), function(i) G[, i, i] / (G[, i, i] + R[, i, i]))
  colnames(draws) <- traits
  summ <- do.call(rbind, lapply(traits, function(tr)
    cbind(data.frame(quantity = paste0("h2:", tr)), .post_summary(draws[, tr], prob))))
  structure(list(summary = summ, draws = draws), class = "breedRB_h2")
}

#' Bind a list of `[n x t x t]` draw arrays along the draw dimension
#' @keywords internal
abind_draws <- function(...) {
  parts <- list(...)
  if (length(parts) == 1L) return(parts[[1]])
  d <- dim(parts[[1]]); t <- d[2]
  n <- sum(vapply(parts, function(a) dim(a)[1], integer(1)))
  out <- array(0, dim = c(n, t, t)); off <- 0
  for (a in parts) { r <- dim(a)[1]; out[(off + 1):(off + r), , ] <- a; off <- off + r }
  out
}
