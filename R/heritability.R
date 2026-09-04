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
  # Whitespace-insensitive fallback: labels are deparsed with spaces (e.g.
  # "gen:leg(x, 4)"), but users routinely type them without ("gen:leg(x,4)").
  norm <- function(z) gsub("[[:space:]]+", "", z)
  nid  <- norm(id)
  hit  <- keys[norm(keys) == nid | norm(labels) == nid]
  if (length(hit)) return(hit[1])
  stop("Term '", id, "' not found. Available: ",
       paste(labels, collapse = ", "), call. = FALSE)
}

#' Per-coefficient heritability of a random-regression term
#'
#' For a `grouping x leg()` reaction-norm term, returns the heritability of each
#' reaction-norm coefficient separately (intercept and every Legendre degree),
#' \eqn{h^2_j = K_{jj} / (K_{jj} + \sigma^2_e)}, where \eqn{K_{jj}} is the
#' across-genotype variance of coefficient \eqn{j} (the realised diagonal of the
#' coefficient (co)variance matrix, computed per MCMC draw as in [rr_gradient()])
#' and \eqn{\sigma^2_e} the residual variance. This exposes the individual
#' per-coefficient heritabilities even though \pkg{BGLR} fits a single shared
#' variance component for the whole `leg()` interaction.
#' @keywords internal
.rr_coef_h2 <- function(fit, term, prob = 0.95) {
  cd   <- .rr_coef_draws(fit, term)
  key  <- .resolve_term(fit, term)
  rr   <- .rr_leg_parts(fit$meta[[key]])
  ikey <- .rr_intercept_key(fit, rr$gc)
  q    <- cd$q
  labs <- c(fit$meta[[ikey]]$label,
            if (q == 1L) fit$meta[[key]]$label
            else paste0(fit$meta[[key]]$label, ":deg", seq_len(q)))
  nD   <- nrow(cd$A)
  varE <- .varE_draws(fit)
  if (length(varE) != nD) {
    ve   <- varcomp(fit)
    varE <- rep(ve$mean[ve$term == "varE"], nD)      # fall back to posterior mean
  }
  # per-draw across-genotype variance of each coefficient (intercept, then slopes)
  Vg <- cbind(apply(cd$A, 1L, stats::var),
              vapply(cd$B, function(m) apply(m, 1L, stats::var), numeric(nD)))
  h2 <- Vg / (Vg + varE)
  colnames(h2) <- labs
  summ <- do.call(rbind, lapply(seq_along(labs), function(j)
    cbind(data.frame(quantity = paste0("h2(", labs[j], ")")),
          .post_summary(h2[, j], prob))))
  structure(list(summary = summ, draws = h2), class = "breedRB_h2")
}

#' Keys of the random terms that carry a genomic (`vm` or `mrk`) component
#' @keywords internal
.vm_keys <- function(fit) {
  names(fit$meta)[vapply(fit$meta, function(m)
    any(vapply(m$components, function(c) c$kind %in% c("vm", "mrk"), logical(1))),
    logical(1))]
}

#' Posterior distribution of heritability
#'
#' Computes narrow-sense (genomic) heritability for every MCMC draw, returning the
#' full posterior distribution rather than a single point estimate:
#' \deqn{h^2 = \sigma^2_{genetic} / (\sigma^2_{genetic} + \sum \sigma^2_{other} + \sigma^2_e).}
#'
#' For a **random regression**, passing a single `grouping x leg()` interaction
#' term as `genetic` (e.g. `heritability(fit, genetic = "gen:leg(x, 4)")`)
#' returns the heritability of **each reaction-norm coefficient** separately —
#' the intercept and every Legendre degree — as
#' \eqn{h^2_j = K_{jj} / (K_{jj} + \sigma^2_e)}, using the realised per-coefficient
#' variances \eqn{K_{jj}} (see [rr_gradient()] / [varcomp()]); this recovers the
#' individual coefficient heritabilities that the single shared \pkg{BGLR}
#' interaction component hides.
#'
#' @param fit A `breedRB_fit`.
#' @param genetic Character vector of genetic term identifiers (label or key)
#'   forming the numerator. Defaults to the `vm()` term(s). A single
#'   random-regression interaction term triggers the per-coefficient breakdown
#'   described above.
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
    if (!length(genetic)) {
      rk    <- .random_keys(fit)
      avail <- unname(vapply(fit$meta[rk], `[[`, character(1), "label"))
      stop("No vm() genetic term was found to use as the heritability numerator.\n",
           "Pass `genetic=` explicitly to pick the random effect, e.g. ",
           "heritability(fit, genetic = ", if (length(avail)) paste0("\"", avail[1], "\"") else "\"<term>\"", ").\n",
           "Available random terms: ",
           if (length(avail)) paste(avail, collapse = ", ") else "(none)",
           call. = FALSE)
    }
  }
  # A single random-regression interaction term -> per-coefficient heritability
  # (intercept + each Legendre degree), from the realised coefficient variances.
  if (length(genetic) == 1L) {
    gk <- .resolve_term(fit, genetic)
    if (!is.null(.rr_leg_parts(fit$meta[[gk]]))) return(.rr_coef_h2(fit, gk, prob))
  }

  gkeys <- vapply(genetic, .resolve_term, character(1), fit = fit)
  vc <- .read_varchains(fit)
  pooled <- do.call(rbind, vc)

  # Expand each logical term key to its variance-file column(s): a split
  # random-regression term contributes one variance component per degree.
  expand <- function(ks) unlist(lapply(ks, function(k) fit$meta[[k]]$eta_keys %||% k))

  rand_cols <- setdiff(colnames(pooled), "varE")
  dkeys <- if (is.null(denominator)) rand_cols else
    expand(vapply(denominator, .resolve_term, character(1), fit = fit))
  gcols <- expand(gkeys)

  num   <- rowSums(pooled[, gcols, drop = FALSE])
  denom <- rowSums(pooled[, union(dkeys, gcols), drop = FALSE]) + pooled[, "varE"]
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
