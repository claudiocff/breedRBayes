# ---------------------------------------------------------------------------
# Design construction: turn parsed term specs into a BGLR ETA list, carrying
# metadata needed to reconstruct effects on the original scale afterwards.
# ---------------------------------------------------------------------------

#' Row-wise Khatri-Rao (face-splitting) product of two design matrices
#'
#' Builds the interaction design of two component matrices: for each row,
#' the Kronecker product of the two rows. Result has `ncol(A) * ncol(B)` columns.
#' @keywords internal
.khatri_rao_rows <- function(A, B) {
  p <- ncol(A); q <- ncol(B)
  out <- A[, rep(seq_len(p), each = q), drop = FALSE] *
         B[, rep(seq_len(q), times = p), drop = FALSE]
  colnames(out) <- as.vector(t(outer(colnames(A), colnames(B), paste, sep = ":")))
  out
}

#' Build the design + metadata for a single interaction component
#' @param exp_var Explained-variance target for the genomic PC rotation (GBLUP
#'   `vm()` / `mrk()`); `NA` keeps every non-null component.
#' @param solo `TRUE` when the component is the whole term (not part of an
#'   interaction). A solo fixed factor is given full "cell-means" coding (one
#'   column per level, no reference dropped) so every level appears in
#'   [solution()]; factors inside an interaction keep treatment contrasts to stay
#'   identifiable against the intercept and the main effects.
#' @keywords internal
.component_design <- function(comp, data, relmat, role, exp_var = NA_real_,
                              solo = TRUE) {
  switch(comp$type,
    factor = {
      if (identical(role, "fixed")) {
        f <- factor(data[[comp$var]])
        if (isTRUE(solo)) {
          # cell-means coding: one column per level (all levels shown). The intercept
          # mu is then aliased with the level means, so BGLR's flat prior splits an
          # overall mean into mu and each level's deviation; add_mu = TRUE in
          # solution() recovers the per-level means.
          X <- model.matrix(~ -1 + f)
          colnames(X) <- levels(f)
          list(X = X, meta = list(kind = "factor", var = comp$var,
                                  levels = levels(f), reference = NULL))
        } else {
          # treatment contrasts (drop reference level) so mu + factor is full rank
          X <- model.matrix(~ f)[, -1, drop = FALSE]
          colnames(X) <- levels(f)[-1]
          list(X = X, meta = list(kind = "factor", var = comp$var,
                                  levels = levels(f), reference = levels(f)[1]))
        }
      } else {
        X <- .incidence(data[[comp$var]], center = TRUE)
        list(X = X, meta = list(kind = "factor", var = comp$var,
                                levels = colnames(X)))
      }
    },
    vm = {
      K <- relmat[[comp$relmat]]
      if (is.null(K)) {
        stop("Relationship matrix '", comp$relmat, "' for vm(", comp$var,
             ") not found in `relmat`.", call. = FALSE)
      }
      lev <- intersect(colnames(K), unique(as.character(data[[comp$var]])))
      K <- K[lev, lev, drop = FALSE]
      rot <- .pc_rotation(K, exp_var = exp_var)
      X <- .expand_rows(rot$PC, data[[comp$var]], lev)   # incidence %*% PC, as a row-gather
      colnames(X) <- if (ncol(X) == length(lev)) lev else paste0("PC", seq_len(ncol(X)))
      list(X = X, meta = list(kind = "vm", var = comp$var, relmat = comp$relmat,
                              levels = lev, method = "GBLUP", bmap = rot$PC))
    },
    mrk = {
      # Marker term with automatic GBLUP / RR-BLUP selection.
      M0 <- relmat[[comp$relmat]]
      if (is.null(M0)) {
        stop("Marker matrix '", comp$relmat, "' for mrk(", comp$var,
             ") not found in `relmat`.", call. = FALSE)
      }
      lev <- intersect(rownames(M0), unique(as.character(data[[comp$var]])))
      M1  <- .clean_markers(M0[lev, , drop = FALSE], var_name = comp$var)  # impute NAs, drop constant markers
      Mc  <- scale(M1, center = TRUE, scale = FALSE)                       # centre markers: Z = M - 2p
      center <- attr(Mc, "scaled:center")         # training column means = 2p (for predicting new genos)
      attr(Mc, "scaled:center") <- NULL
      n <- length(lev); p <- ncol(Mc)
      # VanRaden (2008) method-1 scaling: 2 * sum(p_j q_j), allele freqs p_j from
      # the 0/1/2 dosage column means (center = 2p). Keeps the same denominator for
      # GBLUP and RR-BLUP so the two stay prediction-equivalent.
      pf <- center / 2
      cc <- 2 * sum(pf * (1 - pf))
      if (!is.finite(cc) || cc <= 0) {
        stop("mrk(", comp$var, "): VanRaden scaling 2*sum(pq) is zero/invalid — ",
             "markers appear monomorphic. Expecting a 0/1/2 dosage marker matrix.",
             call. = FALSE)
      }
      method <- if (identical(comp$method, "auto")) {
        if (p >= n) "GBLUP" else "RRBLUP"          # markers >= genotypes -> GBLUP
      } else comp$method
      gblup_rot <- NULL
      if (identical(method, "GBLUP")) {
        # VanRaden G = Mc Mc'/cc + ridge, fitted in its PC basis. Built straight
        # from Mc (avoids forming the n x n G for a low-rank panel when
        # RSpectra is available). A small ridge shifts every eigenvalue up so G is
        # positive-definite (duplicate genotypes / p = n / collinear markers). The
        # eigenpairs of Mc Mc' are cached for the marker back-solve (solve_SNP/predict).
        rot  <- .gblup_rotation(Mc, cc, exp_var = exp_var)
        X    <- .expand_rows(rot$PC, data[[comp$var]], lev)  # incidence %*% PC, as a row-gather
        colnames(X) <- if (ncol(X) == length(lev)) lev else paste0("PC", seq_len(ncol(X)))
        bmap <- rot$PC
        gblup_rot <- list(vectors = rot$vectors, evals_MM = rot$evals_MM,
                          mean_diag_MM = rot$mean_diag_MM)
      } else {
        Msc  <- Mc / sqrt(cc)                      # so var(Msc a) = G * sigma^2_a  (== GBLUP)
        X    <- .expand_rows(Msc, data[[comp$var]], lev)     # incidence %*% Msc, as a row-gather
        colnames(X) <- colnames(Mc)
        bmap <- Msc                                # per-genotype value = b %*% t(bmap)
      }
      message("mrk(", comp$var, "): ", n, " genotypes x ", p, " markers -> ",
              method, if (identical(comp$method, "auto")) " (auto)" else "",
              if (identical(method, "GBLUP") && ncol(bmap) < n)
                paste0(" (rank ", ncol(bmap), ")") else "", ".")
      list(X = X, meta = list(kind = "mrk", var = comp$var, relmat = comp$relmat,
                              levels = lev, method = method, bmap = bmap,
                              markers = colnames(Mc), center = center, c_scale = cc,
                              gblup_rot = gblup_rot))
    },
    leg = {
      x <- as.numeric(data[[comp$var]])
      rng <- range(x, na.rm = TRUE)
      B <- legendre_basis(.scale_unit(x, rng), order = comp$order, orthonormal = TRUE)
      B <- B[, -1, drop = FALSE]              # drop degree 0 (intercept handled separately)
      colnames(B) <- paste0("deg", seq_len(comp$order))
      list(X = B, meta = list(kind = "leg", var = comp$var,
                              order = comp$order, range = rng))
    },
    fa = {
      # Factor-analytic genetic covariance (multi-trait). The design is the plain
      # incidence of the grouping factor; the FA structure is carried on the term
      # via Cov = list(type = "FA", nF = k) and applied by BGLR::Multitrait.
      X <- .incidence(data[[comp$var]])
      list(X = X, meta = list(kind = "fa", var = comp$var, nfac = comp$nfac))
    },
    covariate = {
      X <- matrix(as.numeric(data[[comp$var]]), ncol = 1,
                  dimnames = list(NULL, comp$var))
      list(X = X, meta = list(kind = "covariate", var = comp$var))
    },
    stop("Cannot build design for component type '", comp$type, "'.", call. = FALSE)
  )
}

#' Build the ETA block(s) (design + BGLR model) for a whole term
#'
#' Most terms produce a single BGLR block. A **random regression** — a random
#' term carrying a `leg()` basis of order \eqn{q \ge 2} (e.g. `gen:leg(x, 2)`) —
#' is instead split into one block **per Legendre degree**: `grouping x deg1`,
#' `grouping x deg2`, ... Each block is a separate `BRR` term, so \pkg{BGLR}
#' estimates its **own** variance component. This is the natural diagonal random
#' regression (co)variance structure — a single shared component across degrees
#' (the old single-block build) forces the slope and higher-order coefficients to
#' the same variance, which is a misspecification. (The between-degree covariances
#' still require a multi-trait formulation and are not estimated here.) Fixed
#' `leg()` terms and order-1 random regressions are left as a single block
#' (identical to the unsplit build).
#'
#' @param split_rr Logical; split random `leg()` interactions of order >= 2 into
#'   per-degree blocks. `FALSE` reproduces the single-block behaviour (used for
#'   multi-trait fits).
#' @return A list with `meta` (one logical-term meta record, carrying
#'   `coef_names` and, for a split term, `rr_split`/`rr_order`/`block_coef_names`)
#'   and `blocks` — a list of `{suffix, degree, eta}` BGLR blocks.
#' @keywords internal
.term_eta <- function(term, data, relmat, exp_var = NA_real_, split_rr = TRUE) {
  comp_designs <- lapply(term$components, .component_design,
                         data = data, relmat = relmat, role = term$role,
                         exp_var = exp_var, solo = length(term$components) == 1L)
  model <- if (identical(term$role, "fixed")) "FIXED" else "BRR"
  fa_comp <- Filter(function(cd) identical(cd$meta$kind, "fa"), comp_designs)

  base_meta <- list(label = term$label, role = term$role, model = model,
                    components = lapply(comp_designs, `[[`, "meta"))

  leg_pos <- which(vapply(comp_designs,
                          function(cd) identical(cd$meta$kind, "leg"), logical(1)))
  q <- if (length(leg_pos) == 1L) comp_designs[[leg_pos]]$meta$order else 0L
  do_split <- isTRUE(split_rr) && identical(model, "BRR") &&
              length(leg_pos) == 1L && q >= 2L && !length(fa_comp)

  if (!do_split) {
    # single block: Khatri-Rao product of all components (original behaviour)
    X <- comp_designs[[1]]$X
    if (length(comp_designs) > 1L) {
      for (i in 2:length(comp_designs)) X <- .khatri_rao_rows(X, comp_designs[[i]]$X)
    }
    entry <- list(X = X, model = model)
    if (model == "BRR") entry$saveEffects <- TRUE
    if (length(fa_comp)) entry$Cov <- list(type = "FA", nF = fa_comp[[1]]$meta$nfac)
    base_meta$coef_names <- colnames(X)
    return(list(meta = base_meta,
                blocks = list(list(suffix = "", degree = NA_integer_, eta = entry))))
  }

  # split by Legendre degree: P = Khatri-Rao of the non-leg components (the
  # grouping design, e.g. the genomic PC block or a factor incidence); each
  # degree j then gets its own block  P (x) B[, j]  with its own variance.
  legB   <- comp_designs[[leg_pos]]$X                     # [n x q], colnames deg1..degq
  others <- comp_designs[-leg_pos]
  P <- if (length(others)) {
    Po <- others[[1]]$X
    if (length(others) > 1L) for (i in 2:length(others)) Po <- .khatri_rao_rows(Po, others[[i]]$X)
    Po
  } else NULL
  blocks <- lapply(seq_len(q), function(j) {
    Bj <- legB[, j, drop = FALSE]
    Xj <- if (is.null(P)) Bj else .khatri_rao_rows(P, Bj)
    list(suffix = paste0("deg", j), degree = j,
         eta = list(X = Xj, model = "BRR", saveEffects = TRUE))
  })
  base_meta$coef_names       <- unlist(lapply(blocks, function(b) colnames(b$eta$X)))
  base_meta$block_coef_names <- lapply(blocks, function(b) colnames(b$eta$X))
  base_meta$rr_split <- TRUE
  base_meta$rr_order <- q
  list(meta = base_meta, blocks = blocks)
}

#' Drop rows whose `vm()` factor levels are absent from the relationship matrix
#'
#' GBLUP requires every modelled level to appear in `K`. Rows referencing an
#' unknown level are removed (with a warning) so that all designs stay aligned
#' with the response.
#' @keywords internal
.filter_to_relmat <- function(parsed, data, relmat) {
  keep <- rep(TRUE, nrow(data))
  for (term in c(parsed$fixed, parsed$random)) {
    for (comp in term$components) {
      if (identical(comp$type, "vm")) {
        K <- relmat[[comp$relmat]]
        if (is.null(K)) next
        keep <- keep & as.character(data[[comp$var]]) %in% colnames(K)
      } else if (identical(comp$type, "mrk")) {
        M <- relmat[[comp$relmat]]
        if (is.null(M)) next
        keep <- keep & as.character(data[[comp$var]]) %in% rownames(M)  # genotypes are rows
      }
    }
  }
  if (!all(keep)) {
    warning(sum(!keep), " row(s) dropped: their vm()/mrk() levels are absent from the ",
            "supplied matrix.", call. = FALSE)
  }
  droplevels(data[keep, , drop = FALSE])
}

#' Assemble the full BGLR ETA list from parsed model terms
#'
#' @param parsed A `breedRB_terms` object from [parse_model()].
#' @param data Model data frame.
#' @param relmat Named list of relationship / covariance matrices for `vm()`.
#' @param exp_var Explained-variance target for the genomic PC rotation of
#'   GBLUP-style terms (see [bbglr_control()]); `NA` keeps every non-null
#'   component.
#' @return A list with `ETA` (the BGLR ETA list) and `meta` (per-term metadata,
#'   named by a safe key) plus `group` (residual grouping vector or `NULL`).
#' @keywords internal
build_eta <- function(parsed, data, relmat = list(), exp_var = NA_real_) {
  all_terms <- c(parsed$fixed, parsed$random)
  split_rr  <- !isTRUE(parsed$response$multitrait)     # per-degree split is single-trait only
  built <- lapply(all_terms, .term_eta, data = data, relmat = relmat,
                  exp_var = exp_var, split_rr = split_rr)

  # One base key per logical term; a split random-regression term expands into
  # several BGLR ETA entries (<key>_deg1, <key>_deg2, ...), recorded in the term
  # meta as `eta_keys` so downstream readers can find every block's files.
  base_keys <- make.unique(gsub("[^A-Za-z0-9]", "_",
                                vapply(built, function(b) b$meta$label, character(1))))
  ETA <- list(); meta <- list()
  for (i in seq_along(built)) {
    bk     <- base_keys[i]
    blocks <- built[[i]]$blocks
    m      <- built[[i]]$meta
    if (length(blocks) == 1L && identical(blocks[[1]]$suffix, "")) {
      ETA[[bk]]  <- blocks[[1]]$eta
      m$eta_keys <- bk
    } else {
      bkeys <- paste0(bk, "_", vapply(blocks, `[[`, character(1), "suffix"))
      for (j in seq_along(blocks)) ETA[[bkeys[j]]] <- blocks[[j]]$eta
      m$eta_keys <- bkeys
    }
    meta[[bk]] <- m
  }

  group <- NULL
  if (isTRUE(parsed$residual$hetero)) {
    group <- as.integer(factor(data[[parsed$residual$group]]))
  }
  list(ETA = ETA, meta = meta, group = group)
}
