# ---------------------------------------------------------------------------
# Effect solutions: recover per-level solutions for any model term — random
# BLUPs (with or without a `vm()` genomic component), genomic random-regression
# coefficients, and fixed-effect estimates.
# ---------------------------------------------------------------------------

#' Detect a genomic random-regression term (a `vm()`/`mrk()` x `leg()` interaction)
#'
#' @return `NULL` unless the term is an interaction of exactly one genomic
#'   component with one `leg()` basis; otherwise a list with the genomic component
#'   `gc`, the `leg()` order `q`, and the two components' positions `gpos`/`lpos`.
#' @keywords internal
.geno_leg_parts <- function(meta) {
  comps <- meta$components
  if (length(comps) != 2L) return(NULL)
  kinds <- vapply(comps, `[[`, character(1), "kind")
  gpos  <- which(kinds %in% c("vm", "mrk"))
  lpos  <- which(kinds == "leg")
  if (length(gpos) != 1L || length(lpos) != 1L) return(NULL)
  list(gc = comps[[gpos]], q = comps[[lpos]]$order, gpos = gpos, lpos = lpos)
}

#' Columns of a Khatri-Rao interaction block belonging to `leg()` degree `j`
#'
#' The earlier component in the interaction varies slowest (verified against
#' [.khatri_rao_rows()]): with the genomic component first the degree indexes the
#' fast axis, otherwise the slow axis.
#' @keywords internal
.rr_degree_cols <- function(gpos, lpos, ng, q, j) {
  if (gpos < lpos) seq(j, by = q, length.out = ng)                 # genomic slow, degree fast
  else             ((j - 1L) * ng + 1L):((j - 1L) * ng + ng)       # degree slow, genomic fast
}

#' Per-genotype effect names for a genomic random-regression degree
#' @keywords internal
.rr_effect_names <- function(ids, q, j) if (q == 1L) ids else paste0(ids, ":deg", j)

#' Per-chain posterior draws for a random term, back-mapped to interpretable columns
#'
#' Genomic terms are back-mapped through their principal-component rotation to
#' per-genotype values. A genomic random regression (`mrk(gen, M):leg(x, q)`) is
#' back-mapped **per `leg()` degree**, so each genotype gets one column per degree
#' (`deg1` = linear slope, `deg2` = quadratic, ...); with `q = 1` the columns are
#' just the genotypes. Other genomic interactions fall back to the fitted basis.
#'
#' @return A list (one matrix per chain, `[nSamples x nEffect]`) of post-burn-in
#'   draws. Kept per-chain so the intercept can be aligned before pooling.
#' @keywords internal
.solution_random <- function(fit, key, meta) {
  comps    <- meta$components
  is_geno  <- function(c) c$kind %in% c("vm", "mrk")
  single_g <- length(comps) == 1L && is_geno(comps[[1]])
  has_geno <- any(vapply(comps, is_geno, logical(1)))
  rr       <- .geno_leg_parts(meta)                                # non-NULL for a genomic RR term
  split    <- isTRUE(meta$rr_split)                                # per-degree split blocks
  eks      <- meta$eta_keys %||% key

  lapply(fit$paths, function(prefix) {
    if (split) {
      # Per-degree blocks: each Legendre degree is a separate BGLR term. Read its
      # own draws and (for a genomic grouping) back-map to per-genotype values.
      q <- meta$rr_order
      blocks <- lapply(seq_len(q), function(j) {
        bj <- BGLR::readBinMat(paste0(prefix, "ETA_", eks[j], "_b.bin"))  # [nSamples x p_j]
        if (!is.null(rr)) {
          bmap <- rr$gc$bmap; ids <- rownames(bmap)
          U <- bj %*% t(bmap)                                            # basis coeff -> per-genotype
          colnames(U) <- .rr_effect_names(ids, q, j)
          U
        } else {
          colnames(bj) <- meta$block_coef_names[[j]]                     # "level:degj"
          bj
        }
      })
      return(do.call(cbind, blocks))
    }

    b <- BGLR::readBinMat(paste0(prefix, "ETA_", key, "_b.bin"))   # [nSamples x p]
    if (single_g) {
      bmap <- comps[[1]]$bmap                                      # PC (GBLUP) or markers (RRBLUP)
      out  <- b %*% t(bmap)                                        # -> per-genotype values
      colnames(out) <- rownames(bmap)
      out
    } else if (!is.null(rr)) {
      # Genomic order-1 random regression (single block): back-map the slope.
      bmap <- rr$gc$bmap; ng <- ncol(bmap); q <- rr$q; ids <- rownames(bmap)
      blocks <- lapply(seq_len(q), function(j) {
        bj <- b[, .rr_degree_cols(rr$gpos, rr$lpos, ng, q, j), drop = FALSE]
        U  <- bj %*% t(bmap)
        colnames(U) <- .rr_effect_names(ids, q, j)
        U
      })
      do.call(cbind, blocks)
    } else {
      if (has_geno) {
        warning("Term '", meta$label, "' mixes a genomic component with a ",
                "non-leg() interaction; the returned coordinates are in the fitted ",
                "basis, not per-genotype.", call. = FALSE)
      }
      colnames(b) <- meta$coef_names
      b
    }
  })
}

#' Prior variance of each level/coefficient of a random term
#'
#' The denominator of a BLUP reliability: \eqn{\mathrm{Var}(u_i) = D_{ii}\,\sigma^2_u},
#' where \eqn{\sigma^2_u} is the term's variance component (posterior mean of its
#' `varB` chain) and \eqn{D_{ii}} is the prior relationship diagonal for level `i`.
#' For a genomic term (`u = PC\,b`, `PC PC' = G`) this is \eqn{G_{ii}\,\sigma^2_u}
#' with \eqn{G_{ii} = \sum_k \mathrm{PC}_{ik}^2}; for a plain random factor or a
#' random-regression basis (\eqn{D = I}) it is simply \eqn{\sigma^2_u} for every
#' level.
#' @return A named numeric vector of prior variances, one per effect column.
#' @keywords internal
.random_prior_var <- function(fit, key, meta) {
  vc      <- do.call(rbind, .read_varchains(fit))
  varB_of <- function(k) mean(vc[, k])                          # posterior-mean variance component
  comps   <- meta$components
  is_geno <- function(c) c$kind %in% c("vm", "mrk")
  rr      <- .geno_leg_parts(meta)
  eks     <- meta$eta_keys %||% key
  if (length(comps) == 1L && is_geno(comps[[1]])) {
    bmap  <- comps[[1]]$bmap                                     # PC / marker map, rows = genotypes
    gdiag <- rowSums(bmap^2)                                     # G_ii under the fitted (low-rank) basis
    stats::setNames(gdiag * varB_of(eks[1]), rownames(bmap))
  } else if (!is.null(rr)) {
    # Genomic random regression: Var(u_{g,j}) = G_ii * varB_j. With per-degree
    # split blocks each degree has its own variance component (eks[j]); an order-1
    # RR is a single block (eks[1]).
    ids   <- rownames(rr$gc$bmap); q <- rr$q
    gdiag <- rowSums(rr$gc$bmap^2)
    nm <- character(0); v <- numeric(0)
    for (j in seq_len(q)) {
      kj <- if (isTRUE(meta$rr_split)) eks[j] else eks[1]
      v  <- c(v, gdiag * varB_of(kj)); nm <- c(nm, .rr_effect_names(ids, q, j))
    }
    stats::setNames(v, nm)
  } else if (isTRUE(meta$rr_split)) {
    # Plain-factor random regression (D = I): each degree its own variance.
    bcn <- meta$block_coef_names
    v <- unlist(lapply(seq_along(bcn),
                       function(j) rep(varB_of(eks[j]), length(bcn[[j]]))))
    stats::setNames(v, unlist(bcn))
  } else {
    stats::setNames(rep(varB_of(eks[1]), length(meta$coef_names)), meta$coef_names)
  }
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
#' * **Random terms** return BLUPs. A `vm()`/`mrk()` genomic term is back-mapped
#'   through its principal-component rotation to per-genotype values (identical to
#'   [gebv()]); a plain random factor (fitted without a relationship matrix)
#'   returns one solution per level. A genomic **random regression**
#'   (`mrk(gen, M):leg(x, q)`) returns per-genotype reaction-norm coefficients,
#'   one column per `leg()` degree (`deg1` = linear slope, `deg2` = quadratic,
#'   ...); with `q = 1` the rows are simply the per-genotype slopes. The
#'   **intercept** comes from the genomic main-effect term (`mrk(gen, M)`), so a
#'   full reaction norm is read as
#'   `solution(fit, "mrk(gen, M)")` (intercept) plus
#'   `solution(fit, "mrk(gen, M):leg(x, 1)")` (slope). Non-genomic random
#'   regression / covariate terms return one solution per basis coefficient.
#' * **Fixed terms** return the posterior of each estimable coefficient. A
#'   single fixed factor uses cell-means coding, so **every level** gets a row
#'   (no reference level dropped); the estimates are deviations from the aliased
#'   intercept `mu`, and `add_mu = TRUE` puts them on the per-level mean scale.
#'   A factor **inside an interaction** keeps treatment contrasts (its reference
#'   level is the implicit baseline at 0) so the term stays identifiable.
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
#'   (posterior mean), `sd`, `lower`, `upper`. **Random terms** additionally carry
#'   `pev` — the prediction error variance, i.e. the posterior variance of the
#'   effect (computed on the deviation scale, before any `add_mu` shift) — and
#'   `reliability`, the BLUP reliability
#'   \eqn{r^2_i = 1 - \mathrm{PEV}_i / \mathrm{Var}(u_i)} bounded to `[0, 1]`, where
#'   the prior variance \eqn{\mathrm{Var}(u_i)} is \eqn{G_{ii}\sigma^2_u} for a
#'   genomic term and \eqn{\sigma^2_u} for a plain random factor / basis
#'   coefficient. Random-term rows are ordered by decreasing `solution`;
#'   fixed-term rows keep design order. The pooled draws matrix
#'   (`nDraws x nEffect`) is attached as attribute `"draws"`, and the resolved role
#'   as attribute `"type"`.
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

  raw <- do.call(rbind, draws_list)          # pre-mu draws — PEV is a property of the effect, not mu

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
  if (identical(role, "random")) {
    # Prediction error variance (posterior variance of the effect) and BLUP
    # reliability r^2 = 1 - PEV / Var(u_i), using the pre-mu draws for PEV.
    pev       <- apply(raw, 2, stats::var)
    prior_var <- .random_prior_var(fit, key, meta)[colnames(raw)]
    rel       <- pmin(pmax(1 - pev / prior_var, 0), 1)   # bounded to [0, 1] against MC noise
    out$pev         <- unname(pev)
    out$reliability <- unname(rel)
    out <- out[order(-out$solution), ]
  }
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

#' Back-solve per-marker effects from per-genotype breeding-value draws (GBLUP)
#'
#' Given a matrix of per-draw genotype breeding values `U` (`[nDraws x nGen]`,
#' columns named by genotype), recover the centred-marker allele-substitution
#' effects \eqn{b = M_c^\top (M_c M_c^\top)^{-1} u} for every draw. Shared by the
#' genomic main effect ([.marker_draws()]) and each random-regression coefficient
#' ([.rr_marker_draws()]).
#' @param gc The genomic component metadata (`vm`/`mrk`), carrying `markers`,
#'   `center` and the cached `gblup_rot` eigenpairs.
#' @param U Per-draw breeding values, `[nDraws x nGen]`, columns = genotype IDs.
#' @param M Training marker matrix (genotypes in rows). Required.
#' @keywords internal
.markers_from_U <- function(gc, U, M) {
  M <- as.matrix(M)
  if (is.null(rownames(M))) stop("`M` must have genotype IDs as row names.", call. = FALSE)
  ids  <- colnames(U)
  miss <- setdiff(ids, rownames(M))
  if (length(miss)) {
    stop("`M` is missing ", length(miss), " genotype(s) present in the fit, e.g. '",
         miss[1], "'.", call. = FALSE)
  }
  mk <- gc$markers
  if (is.null(mk)) {
    # vm() fit: no markers were retained by the model, centre the supplied M as-is
    Mc <- scale(M[ids, , drop = FALSE], center = TRUE, scale = FALSE)
    attr(Mc, "scaled:center") <- NULL
  } else {
    # mrk() fit: reproduce the exact cleaned + centred training markers used to
    # build G (retained markers only, centred by the training means, NAs imputed).
    miss_mk <- setdiff(mk, colnames(M))
    if (length(miss_mk)) {
      stop("`M` is missing ", length(miss_mk), " marker(s) used in the fit, e.g. '",
           miss_mk[1], "'.", call. = FALSE)
    }
    Mc <- sweep(M[ids, mk, drop = FALSE], 2L, gc$center[mk], "-")
    Mc[is.na(Mc)] <- 0                                # same mean-imputation as fitting
  }
  rot <- gc$gblup_rot
  if (!is.null(rot)) {
    # Reuse the cached fit-time eigendecomposition of Mc Mc' — no re-inversion.
    # b = Mc' (Mc Mc' + eps I)^{-1} u, and each u lies in span(V) (u = PC b_brr),
    # so only the retained eigenpairs contribute:
    #   (Mc Mc' + eps I)^{-1} u = V diag(1/(lambda + eps)) V' u.
    # (Directions in the left-null space of Mc drop out since V' Mc = 0 there.)
    V    <- rot$vectors[ids, , drop = FALSE]          # eigenvectors of Mc Mc', rows = genotypes
    eps  <- 1e-8 * rot$mean_diag_MM                   # == 1e-8 * mean(diag(Mc Mc')) used before
    VtMc <- crossprod(V, Mc)                          # r x p  (V' Mc)
    draws <- (U %*% V) %*% (VtMc / (rot$evals_MM + eps))   # row-scale V'Mc by 1/(lambda+eps)
  } else {
    # vm() fit (no cached rotation): fall back to a direct solve.
    MM   <- tcrossprod(Mc)
    Minv <- solve(MM + diag(1e-8 * mean(diag(MM)), nrow(MM)))   # tiny ridge for stability
    draws <- U %*% (Minv %*% Mc)
  }
  colnames(draws) <- colnames(Mc)
  draws
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
  sol_list <- .solution_random(fit, key, meta)        # per-chain [nDraws x nGen]
  U     <- do.call(rbind, sol_list)
  draws <- .markers_from_U(gc, U, M)
  attr(draws, "n_per_chain") <- vapply(sol_list, nrow, integer(1))
  draws
}

#' Pooled per-marker effects of a genomic random-regression term, per `leg()` degree
#'
#' The marker-effect counterpart of the random-regression back-map in
#' [.solution_random()]: splits the interaction coefficients by `leg()` degree and
#' returns, for each degree, the per-marker allele-substitution effects
#' (back-solved for GBLUP, read directly for RR-BLUP). Columns are the markers
#' (suffixed `:deg{j}` when the `leg()` order exceeds 1).
#' @keywords internal
.rr_marker_draws <- function(fit, key, M = NULL) {
  meta   <- fit$meta[[key]]
  rr     <- .geno_leg_parts(meta)
  gc     <- rr$gc; q <- rr$q; bmap <- gc$bmap; ng <- ncol(bmap)
  method <- if (!is.null(gc$method)) gc$method else "GBLUP"
  if (identical(method, "GBLUP")) {
    if (is.null(M)) M <- fit$relmat[[gc$relmat]]
    if (is.null(M)) {
      stop("Training marker matrix not found in the fit; pass `M` to back-solve ",
           "random-regression SNP effects from this GBLUP fit.", call. = FALSE)
    }
  }
  split <- isTRUE(meta$rr_split)
  eks   <- meta$eta_keys %||% key
  # Per-degree draws: separate blocks when split, else slice the single block.
  per_by_deg <- if (split) {
    lapply(seq_len(q), function(j) do.call(rbind, lapply(fit$paths, function(prefix)
      BGLR::readBinMat(paste0(prefix, "ETA_", eks[j], "_b.bin")))))
  } else {
    b <- do.call(rbind, lapply(fit$paths, function(prefix)
      BGLR::readBinMat(paste0(prefix, "ETA_", key, "_b.bin"))))   # [nDraws x (ng * q)]
    lapply(seq_len(q), function(j) b[, .rr_degree_cols(rr$gpos, rr$lpos, ng, q, j), drop = FALSE])
  }
  n_per_chain <- vapply(fit$paths, function(prefix)
    nrow(BGLR::readBinMat(paste0(prefix, "ETA_", eks[1], "_b.bin"))), integer(1))
  blocks <- lapply(seq_len(q), function(j) {
    bj <- per_by_deg[[j]]
    d  <- if (identical(method, "RRBLUP")) {
      out <- bj / sqrt(gc$c_scale); colnames(out) <- gc$markers; out
    } else {
      U <- bj %*% t(bmap); colnames(U) <- rownames(bmap)  # basis coeff -> per-genotype
      .markers_from_U(gc, U, M)                            # then back-solve markers
    }
    if (q > 1L) colnames(d) <- paste0(colnames(d), ":deg", j)
    d
  })
  draws <- do.call(cbind, blocks)
  attr(draws, "n_per_chain") <- n_per_chain
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
#'   genomic term. For a genomic **random-regression** term
#'   (`mrk(gen, M):leg(x, q)`) the marker effects are returned **per `leg()`
#'   degree**: one effect per marker per degree, with markers suffixed `:deg{j}`
#'   when `q > 1` (the intercept marker effects come from the genomic main-effect
#'   term, e.g. `solve_SNP(fit, term = "mrk(gen, M)")`).
#' @param prob Central credible-interval mass (default 0.95).
#' @return A data frame with `marker`, `effect` (posterior mean), `sd`, `lower`,
#'   `upper`, in marker (design) order, with the pooled draws matrix
#'   (`nDraws x nMarker`) attached as attribute `"draws"` and the fitted method
#'   as attribute `"method"`.
#' @examples
#' \donttest{
#' snp <- solve_SNP(fit)      # marker effects (back-solved for GBLUP, direct for RR-BLUP)
#' head(snp[order(-abs(snp$effect)), ])
#' # genomic random regression: intercept vs slope marker effects
#' int_snp   <- solve_SNP(fit, term = "mrk(gen, M)")            # intercept
#' slope_snp <- solve_SNP(fit, term = "mrk(gen, M):leg(x, 1)")  # slope (deg1)
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

  # A genomic random-regression term (mrk():leg()) is back-solved per leg() degree;
  # a plain genomic term yields one marker effect per marker.
  draws <- if (!is.null(.geno_leg_parts(meta))) .rr_marker_draws(fit, key, M)
           else                                 .marker_draws(fit, key, M)

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
#'   row names, markers in columns. Must contain every marker retained in the fit
#'   (extra columns are ignored; column order need not match). Missing calls are
#'   imputed with the training marker means. May include training genotypes, in
#'   which case the prediction reproduces their fitted value.
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
  Mc_new[is.na(Mc_new)] <- 0                             # impute missing new calls with the training mean
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
