# ---------------------------------------------------------------------------
# Effect solutions: recover per-level solutions for any model term — random
# BLUPs (with or without a `vm()` genomic component) and fixed-effect estimates.
# ---------------------------------------------------------------------------

#' Per-chain posterior draws for a random term, back-mapped to interpretable columns
#'
#' @return A list (one matrix per chain, `[nSamples x nEffect]`) of post-burn-in
#'   draws. Kept per-chain so the intercept can be aligned before pooling.
#' @keywords internal
.solution_random <- function(fit, key, meta) {
  comps    <- meta$components
  is_geno  <- function(c) c$kind %in% c("vm", "mrk")
  single_g <- length(comps) == 1L && is_geno(comps[[1]])
  has_geno <- any(vapply(comps, is_geno, logical(1)))

  lapply(fit$paths, function(prefix) {
    b <- BGLR::readBinMat(paste0(prefix, "ETA_", key, "_b.bin"))   # [nSamples x p]
    if (single_g) {
      bmap <- comps[[1]]$bmap                                      # PC (GBLUP) or markers (RRBLUP)
      out  <- b %*% t(bmap)                                        # -> per-genotype values
      colnames(out) <- rownames(bmap)
      out
    } else {
      if (has_geno) {
        warning("Term '", meta$label, "' has a genomic component inside an interaction; ",
                "the returned coordinates are in the fitted basis, not per-genotype. ",
                "Use gebv() on the genomic main effect for breeding values.",
                call. = FALSE)
      }
      colnames(b) <- meta$coef_names
      b
    }
  })
}

#' Per-chain posterior draws for a fixed term, read from BGLR's `_b.dat` trace
#'
#' @return A list (one matrix per chain) of post-burn-in draws.
#' @keywords internal
.solution_fixed <- function(fit, key, meta) {
  nburn <- floor(fit$control$burnIn / fit$control$thin)
  lapply(fit$paths, function(prefix) {
    tr <- as.matrix(utils::read.table(paste0(prefix, "ETA_", key, "_b.dat"),
                                       header = TRUE))
    if (nrow(tr) > nburn) tr <- tr[(nburn + 1L):nrow(tr), , drop = FALSE]  # drop burn-in
    if (!is.null(meta$coef_names) && ncol(tr) == length(meta$coef_names)) {
      colnames(tr) <- meta$coef_names
    }
    tr
  })
}

#' Per-chain post-burn-in draws of the intercept `mu`
#'
#' Reads BGLR's `mu.dat` trace (which includes the thinned burn-in) and aligns
#' each chain to a target number of draws by keeping the trailing rows.
#' @param n_per_chain Integer vector of target draw counts, one per chain.
#' @keywords internal
.solution_mu <- function(fit, n_per_chain) {
  Map(function(prefix, n) {
    f <- paste0(prefix, "mu.dat")
    if (!file.exists(f)) {
      stop("add_mu = TRUE but the intercept trace '", f, "' was not found; ",
           "the model appears to have been fitted without an intercept.",
           call. = FALSE)
    }
    mu <- scan(f, quiet = TRUE)
    utils::tail(mu, n)                                             # align to post-burn-in draws
  }, fit$paths, n_per_chain)
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
#' @param add_mu Logical (default `FALSE`). If `TRUE`, the model intercept `mu` is
#'   added to every draw, shifting the solutions onto the overall-mean scale
#'   (e.g. `BLUP + mu`). Most useful for random BLUPs; the addition is done
#'   per posterior draw so the credible intervals account for uncertainty in `mu`.
#' @return A data frame with `effect` (level / coefficient name), `solution`
#'   (posterior mean), `sd`, `lower`, `upper`. Random-term rows are ordered by
#'   decreasing `solution`; fixed-term rows keep design order. The pooled draws
#'   matrix (`nDraws x nEffect`) is attached as attribute `"draws"`, and the
#'   resolved role as attribute `"type"`.
#' @examples
#' \donttest{
#' solution(fit, term = "gen", type = "random")               # random BLUPs (no G matrix)
#' solution(fit, term = "gen", type = "random", add_mu = TRUE) # BLUP + mu
#' solution(fit, term = "env", type = "fixed")                 # fixed-effect estimates
#' }
#' @seealso [gebv()] for genomic breeding values from a `vm()` term.
#' @export
solution <- function(fit, term = NULL, type = NULL, prob = 0.95, add_mu = FALSE) {
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

  draws_list <- if (identical(role, "random")) .solution_random(fit, key, meta)
                else                          .solution_fixed(fit, key, meta)

  if (isTRUE(add_mu)) {
    mu_list    <- .solution_mu(fit, vapply(draws_list, nrow, integer(1)))
    draws_list <- Map(function(d, mu) d + mu, draws_list, mu_list)  # mu added per draw, per column
  }

  draws <- do.call(rbind, draws_list)

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
  vmc  <- Filter(function(c) c$kind %in% c("vm", "mrk"), meta$components)
  if (!length(vmc)) {
    stop("Term '", term, "' has no vm()/mrk() genomic component. Use ",
         "solution(fit, term = \"", term, "\", type = \"random\") for its BLUPs instead.",
         call. = FALSE)
  }

  out <- solution(fit, term = key, type = "random", prob = prob)
  names(out)[names(out) == "effect"]   <- "ID"
  names(out)[names(out) == "solution"] <- "gebv"
  out
}

#' Pooled posterior draws of per-marker effects on the centred-marker scale
#'
#' Returns an `[nDraws x nMarker]` matrix of allele-substitution effects \eqn{b}
#' such that \eqn{GEBV = M_c\, b}, pooled across chains. For a RR-BLUP fit the
#' effects are read directly (and rescaled); for a GBLUP fit they are back-solved
#' from the breeding-value draws. The per-chain draw counts are attached as
#' attribute `"n_per_chain"` so callers can align a per-draw intercept.
#' @keywords internal
.marker_draws <- function(fit, key, M = NULL) {
  meta   <- fit$meta[[key]]
  gc     <- Filter(function(c) c$kind %in% c("vm", "mrk"), meta$components)[[1]]
  method <- if (!is.null(gc$method)) gc$method else "GBLUP"

  if (identical(method, "RRBLUP")) {
    per   <- lapply(fit$paths, function(prefix)
      BGLR::readBinMat(paste0(prefix, "ETA_", key, "_b.bin")))
    draws <- do.call(rbind, per) / sqrt(gc$c_scale)   # undo the sqrt(c) design scaling
    colnames(draws) <- gc$markers
    attr(draws, "n_per_chain") <- vapply(per, nrow, integer(1))
    return(draws)
  }

  # GBLUP: back-solve b = Mc' (Mc Mc')^{-1} u per draw of the breeding values u.
  if (is.null(M)) M <- fit$relmat[[gc$relmat]]        # reuse the training markers held by the fit
  if (is.null(M)) {
    stop("Training marker matrix not found in the fit; pass `M` to back-solve SNP effects ",
         "from this GBLUP fit.", call. = FALSE)
  }
  M <- as.matrix(M)
  if (is.null(rownames(M))) stop("`M` must have genotype IDs as row names.", call. = FALSE)
  sol_list <- .solution_random(fit, key, meta)        # per-chain [nDraws x nGen]
  U   <- do.call(rbind, sol_list)
  ids <- colnames(U)
  miss <- setdiff(ids, rownames(M))
  if (length(miss)) {
    stop("`M` is missing ", length(miss), " genotype(s) present in the fit, e.g. '",
         miss[1], "'.", call. = FALSE)
  }
  Mc <- scale(M[ids, , drop = FALSE], center = TRUE, scale = FALSE)
  attr(Mc, "scaled:center") <- NULL
  MM   <- tcrossprod(Mc)
  Minv <- solve(MM + diag(1e-8 * mean(diag(MM)), nrow(MM)))   # tiny ridge for stability
  draws <- U %*% (Minv %*% Mc)
  colnames(draws) <- colnames(Mc)
  attr(draws, "n_per_chain") <- vapply(sol_list, nrow, integer(1))
  draws
}

#' Marker (SNP) effects from a genomic model
#'
#' Recovers the posterior distribution of per-marker allele-substitution effects
#' \eqn{b} on the centred-marker scale, such that the genomic breeding values
#' satisfy \eqn{GEBV = M_c\, b} (where \eqn{M_c} is the column-centred marker
#' matrix). This lets a **GBLUP** fit — which estimates breeding values rather
#' than marker effects — be back-transformed to marker effects:
#' \deqn{b = M_c^\top (M_c M_c^\top)^{-1}\, u}
#' applied to every posterior draw of the breeding values \eqn{u}. For a
#' **RR-BLUP** fit (`mrk()` chose RR-BLUP because genotypes outnumbered markers)
#' the marker effects are estimated directly and are simply rescaled and returned.
#'
#' The back-transform requires \eqn{M_c M_c^\top} to be invertible, which holds
#' when the number of markers is at least the number of genotypes — exactly the
#' regime in which `mrk()` selects GBLUP.
#'
#' @param fit A `breedRB_fit` (single-trait) with a `vm()` or `mrk()` genomic term.
#' @param M Marker matrix with genotypes in rows (row names matching the genotype
#'   IDs) and markers in columns. Only needed to override the training markers
#'   when back-solving a GBLUP fit; for a `mrk()` fit the training markers stored
#'   in the fit are used automatically, and a RR-BLUP fit ignores `M`.
#' @param term Genomic term identifier (label or key). Defaults to the first
#'   genomic term.
#' @param prob Central credible-interval mass (default 0.95).
#' @return A data frame with `marker`, `effect` (posterior mean), `sd`, `lower`,
#'   `upper`, in marker (design) order, with the pooled draws matrix
#'   (`nDraws x nMarker`) attached as attribute `"draws"` and the fitted method
#'   as attribute `"method"`.
#' @examples
#' \donttest{
#' snp <- solve_SNP(fit)      # marker effects (back-solved for GBLUP, direct for RR-BLUP)
#' head(snp[order(-abs(snp$effect)), ])
#' }
#' @seealso [predict.breedRB_fit()] to score new genotypes, [solution()].
#' @export
solve_SNP <- function(fit, M = NULL, term = NULL, prob = 0.95) {
  stopifnot(inherits(fit, "breedRB_fit"))
  if (isTRUE(fit$response$multitrait)) {
    stop("solve_SNP() currently supports single-trait fits.", call. = FALSE)
  }
  if (is.null(term)) {
    vk <- .vm_keys(fit)
    if (!length(vk)) {
      stop("No vm()/mrk() genomic term found in the fit.", call. = FALSE)
    }
    term <- vk[1]
  }
  key   <- .resolve_term(fit, term)
  meta  <- fit$meta[[key]]
  gcomp <- Filter(function(c) c$kind %in% c("vm", "mrk"), meta$components)
  if (!length(gcomp)) {
    stop("Term '", term, "' has no vm()/mrk() genomic component.", call. = FALSE)
  }
  gc     <- gcomp[[1]]
  method <- if (!is.null(gc$method)) gc$method else "GBLUP"

  draws <- .marker_draws(fit, key, M)

  a <- (1 - prob) / 2
  out <- data.frame(
    marker = colnames(draws),
    effect = colMeans(draws),
    sd     = apply(draws, 2, stats::sd),
    lower  = apply(draws, 2, stats::quantile, probs = a),
    upper  = apply(draws, 2, stats::quantile, probs = 1 - a),
    row.names = NULL
  )
  attr(out, "draws")  <- draws
  attr(out, "method") <- method
  out
}

#' Predict genomic values for new genotypes from their markers
#'
#' Scores genotypes that were **not** in the training data by combining the
#' posterior marker effects with the new genotypes' marker calls:
#' \deqn{\hat{g}_{new} = M_{c,new}\, b,}
#' where \eqn{M_{c,new}} is `M_new` centred by the **training** column means (so
#' the new genotypes sit on the same scale as the fitted breeding values) and
#' \eqn{b} are the per-draw marker effects (read directly for a RR-BLUP fit,
#' back-solved for a GBLUP fit — see [solve_SNP()]). The full posterior is
#' propagated, so every returned prediction carries a credible interval.
#'
#' This requires a marker-based `mrk()` term. A `vm()` fit holds only a
#' relationship matrix and cannot score genotypes outside it; refit the genomic
#' term as `mrk(gen, M)` to enable prediction.
#'
#' @param object A `breedRB_fit` (single-trait) with a `mrk()` genomic term.
#' @param M_new Marker matrix for the genotypes to predict: genotype IDs in the
#'   row names, markers in columns. Must contain every marker used in the fit
#'   (extra columns are ignored; column order need not match). May include
#'   training genotypes, in which case the prediction reproduces their fitted
#'   value.
#' @param term Genomic term identifier (label or key). Defaults to the first
#'   genomic term.
#' @param add_mu Logical (default `TRUE`). Add the model intercept `mu` to every
#'   draw, putting predictions on the overall-mean (phenotype) scale rather than
#'   the deviation scale. Done per posterior draw, so `mu`'s uncertainty enters
#'   the interval.
#' @param prob Central credible-interval mass (default 0.95).
#' @param ... Unused; for S3 compatibility.
#' @return A data frame with `ID`, `prediction` (posterior mean), `sd`, `lower`,
#'   `upper`, ordered by decreasing `prediction`, with the pooled draws matrix
#'   (`nDraws x nGenotype`) attached as attribute `"draws"`.
#' @examples
#' \donttest{
#' fit  <- bbglr(y ~ 1, random = ~ mrk(gen, M), data = dat, relmat = list(M = M))
#' pred <- predict(fit, M_new, add_mu = TRUE)   # score unobserved genotypes
#' head(pred)
#' }
#' @seealso [solve_SNP()], [solution()].
#' @export
predict.breedRB_fit <- function(object, M_new, term = NULL, add_mu = TRUE,
                                prob = 0.95, ...) {
  fit <- object
  stopifnot(inherits(fit, "breedRB_fit"))
  if (isTRUE(fit$response$multitrait)) {
    stop("predict() currently supports single-trait fits.", call. = FALSE)
  }
  if (is.null(term)) {
    vk <- .vm_keys(fit)
    if (!length(vk)) {
      stop("predict() needs a genomic mrk() term to score new genotypes; ",
           "none was found in the fit.", call. = FALSE)
    }
    term <- vk[1]
  }
  key   <- .resolve_term(fit, term)
  meta  <- fit$meta[[key]]
  gcomp <- Filter(function(c) c$kind %in% c("vm", "mrk"), meta$components)
  if (!length(gcomp)) {
    stop("Term '", term, "' has no genomic component.", call. = FALSE)
  }
  gc <- gcomp[[1]]
  if (is.null(gc$markers) || is.null(gc$center)) {
    stop("predict() requires a marker-based term (mrk(gen, M)); a vm() fit holds ",
         "only a relationship matrix and cannot score new genotypes. Refit the ",
         "genomic term as mrk(gen, M).", call. = FALSE)
  }

  M_new <- as.matrix(M_new)
  if (is.null(rownames(M_new))) {
    stop("`M_new` must have genotype IDs as row names.", call. = FALSE)
  }
  miss <- setdiff(gc$markers, colnames(M_new))
  if (length(miss)) {
    stop("`M_new` is missing ", length(miss), " marker(s) used in the fit, e.g. '",
         miss[1], "'.", call. = FALSE)
  }

  B  <- .marker_draws(fit, key)                          # [nDraws x nMarker], cols = markers
  Xn <- M_new[, gc$markers, drop = FALSE]                # align to training marker order
  Mc_new <- sweep(Xn, 2L, gc$center[gc$markers], "-")    # centre by TRAINING means
  preds  <- t(Mc_new %*% t(B))                           # [nDraws x nGenotype]

  if (isTRUE(add_mu)) {
    mu <- unlist(.solution_mu(fit, attr(B, "n_per_chain")), use.names = FALSE)
    preds <- preds + mu                                  # mu[draw] added to each genotype's draw
  }

  a <- (1 - prob) / 2
  out <- data.frame(
    ID         = colnames(preds),
    prediction = colMeans(preds),
    sd         = apply(preds, 2, stats::sd),
    lower      = apply(preds, 2, stats::quantile, probs = a),
    upper      = apply(preds, 2, stats::quantile, probs = 1 - a),
    row.names  = NULL
  )
  out <- out[order(-out$prediction), ]
  attr(out, "draws")  <- preds
  attr(out, "add_mu") <- isTRUE(add_mu)
  out
}
