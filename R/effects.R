# ---------------------------------------------------------------------------
# Effect solutions: recover per-level solutions for any model term — random
# BLUPs (with or without a `vm()` genomic component) and fixed-effect estimates.
# ---------------------------------------------------------------------------

#' Posterior draws for a random term, back-mapped to interpretable columns
#' @keywords internal
.solution_random <- function(fit, key, meta) {
  comps <- meta$components
  is_vm_single <- length(comps) == 1L && identical(comps[[1]]$kind, "vm")
  has_vm       <- any(vapply(comps, function(c) identical(c$kind, "vm"), logical(1)))

  do.call(rbind, lapply(fit$paths, function(prefix) {
    b <- BGLR::readBinMat(paste0(prefix, "ETA_", key, "_b.bin"))   # [nSamples x p]
    if (is_vm_single) {
      PC  <- comps[[1]]$pc
      out <- b %*% t(PC)                                           # PC space -> levels
      colnames(out) <- rownames(PC)
      out
    } else {
      if (has_vm) {
        warning("Term '", meta$label, "' has a vm() component inside an interaction; ",
                "the returned coordinates are in the genomic PC basis, not per-genotype. ",
                "Use gebv() on the genomic main effect for breeding values.",
                call. = FALSE)
      }
      colnames(b) <- meta$coef_names
      b
    }
  }))
}

#' Posterior draws for a fixed term, read from BGLR's `_b.dat` trace
#' @keywords internal
.solution_fixed <- function(fit, key, meta) {
  nburn <- floor(fit$control$burnIn / fit$control$thin)
  draws <- do.call(rbind, lapply(fit$paths, function(prefix) {
    tr <- as.matrix(utils::read.table(paste0(prefix, "ETA_", key, "_b.dat"),
                                       header = TRUE))
    if (nrow(tr) > nburn) tr <- tr[(nburn + 1L):nrow(tr), , drop = FALSE]  # drop burn-in
    tr
  }))
  if (!is.null(meta$coef_names) && ncol(draws) == length(meta$coef_names)) {
    colnames(draws) <- meta$coef_names
  }
  draws
}

#' Posterior solutions (BLUPs / fixed effects) for a model term
#'
#' Extracts the posterior distribution of the effect solutions for any term in a
#' fitted model, pooled across chains:
#' * **Random terms** return BLUPs. A `vm()` genomic term is back-mapped through
#'   its principal-component rotation to per-genotype values (identical to
#'   [gebv()]); a plain random factor (fitted without a relationship matrix)
#'   returns one solution per level; random regression / covariate terms return
#'   one solution per basis coefficient.
#' * **Fixed terms** return the posterior of each estimable coefficient (treatment
#'   contrasts; the reference level is the implicit baseline at 0).
#'
#' @param fit A `breedRB_fit` (single-trait).
#' @param term Term identifier (the label as written in the formula, e.g. `"gen"`
#'   or `"env"`, or the internal key). If `NULL`, the first random term is used.
#' @param type Optional; one of `"random"` or `"fixed"`. When supplied it is
#'   validated against the term's actual role (a mismatch is an error), which is a
#'   convenient guard when a name could be ambiguous.
#' @param prob Central credible-interval mass (default 0.95).
#' @return A data frame with `effect` (level / coefficient name), `solution`
#'   (posterior mean), `sd`, `lower`, `upper`. Random-term rows are ordered by
#'   decreasing `solution`; fixed-term rows keep design order. The pooled draws
#'   matrix (`nDraws x nEffect`) is attached as attribute `"draws"`, and the
#'   resolved role as attribute `"type"`.
#' @examples
#' \donttest{
#' solution(fit, term = "gen", type = "random")   # random BLUPs (no G matrix needed)
#' solution(fit, term = "env", type = "fixed")    # fixed-effect estimates
#' }
#' @seealso [gebv()] for genomic breeding values from a `vm()` term.
#' @export
solution <- function(fit, term = NULL, type = NULL, prob = 0.95) {
  stopifnot(inherits(fit, "breedRB_fit"))
  if (isTRUE(fit$response$multitrait)) {
    stop("solution() currently supports single-trait fits.", call. = FALSE)
  }

  if (is.null(term)) {
    rk <- .random_keys(fit)
    if (!length(rk)) {
      stop("No random term to extract; pass `term=` (and optionally `type=`) explicitly.",
           call. = FALSE)
    }
    term <- fit$meta[[rk[1]]]$label
  }
  key  <- .resolve_term(fit, term)
  meta <- fit$meta[[key]]
  role <- meta$role

  if (!is.null(type)) {
    type <- match.arg(type, c("random", "fixed"))
    if (!identical(type, role)) {
      stop("Term '", term, "' is a ", role, " effect, but type = \"", type,
           "\" was requested.", call. = FALSE)
    }
  }

  draws <- if (identical(role, "random")) .solution_random(fit, key, meta)
           else                          .solution_fixed(fit, key, meta)

  a <- (1 - prob) / 2
  out <- data.frame(
    effect   = colnames(draws),
    solution = colMeans(draws),
    sd       = apply(draws, 2, stats::sd),
    lower    = apply(draws, 2, stats::quantile, probs = a),
    upper    = apply(draws, 2, stats::quantile, probs = 1 - a),
    row.names = NULL
  )
  if (identical(role, "random")) out <- out[order(-out$solution), ]
  attr(out, "draws") <- draws
  attr(out, "type")  <- role
  out
}

#' Posterior genomic estimated breeding values (GEBV)
#'
#' @description
#' **Deprecated.** Use [solution()], which extracts posterior solutions for any
#' term (genomic `vm()` BLUPs, plain random-factor BLUPs, and fixed effects).
#' `gebv(fit, term)` is equivalent to `solution(fit, term, type = "random")` with
#' the value column renamed `gebv`.
#'
#' It still works for a `vm()` genomic term for now, mapping the ridge
#' coefficients back through the PC rotation (`genetic value = PC %*% b`).
#'
#' @param fit A `breedRB_fit` (single-trait).
#' @param term Genetic term identifier (label or key). Defaults to the first
#'   `vm()` term. For a non-genomic random effect use [solution()] instead.
#' @param prob Central credible-interval mass (default 0.95).
#' @return A data frame with `ID`, posterior `gebv` (mean), `sd`, `lower`, `upper`,
#'   ordered by decreasing `gebv`, with the pooled draws matrix
#'   (`nDraws x nLevel`) attached as attribute `"draws"`.
#' @examples
#' \donttest{ head(gebv(fit)) }
#' @seealso [solution()], the general replacement.
#' @export
gebv <- function(fit, term = NULL, prob = 0.95) {
  .Deprecated("solution", package = "breedRBayes",
              msg = "gebv() is deprecated; use solution(fit, term, type = \"random\") instead.")
  stopifnot(inherits(fit, "breedRB_fit"))
  if (isTRUE(fit$response$multitrait)) {
    stop("gebv() currently supports single-trait fits; use the per-trait chains directly.",
         call. = FALSE)
  }
  if (is.null(term)) {
    vk <- .vm_keys(fit)
    if (!length(vk)) {
      stop("No vm() genomic term found. For a random effect fitted without a ",
           "relationship matrix, use solution(fit, term = ..., type = \"random\").",
           call. = FALSE)
    }
    term <- vk[1]
  }
  key  <- .resolve_term(fit, term)
  meta <- fit$meta[[key]]
  vmc  <- Filter(function(c) identical(c$kind, "vm"), meta$components)
  if (!length(vmc)) {
    stop("Term '", term, "' has no vm() genomic component. Use ",
         "solution(fit, term = \"", term, "\", type = \"random\") for its BLUPs instead.",
         call. = FALSE)
  }

  out <- solution(fit, term = key, type = "random", prob = prob)
  names(out)[names(out) == "effect"]   <- "ID"
  names(out)[names(out) == "solution"] <- "gebv"
  out
}
