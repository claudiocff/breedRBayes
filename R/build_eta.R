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
#' @keywords internal
.component_design <- function(comp, data, relmat, role) {
  switch(comp$type,
    factor = {
      if (identical(role, "fixed")) {
        # treatment contrasts (drop reference level) so mu + factor is full rank
        f <- factor(data[[comp$var]])
        X <- model.matrix(~ f)[, -1, drop = FALSE]
        colnames(X) <- levels(f)[-1]
        list(X = X, meta = list(kind = "factor", var = comp$var,
                                levels = levels(f), reference = levels(f)[1]))
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
      rot <- .pc_rotation(K)
      Z <- .incidence(factor(data[[comp$var]], levels = lev))
      X <- Z %*% rot$PC
      colnames(X) <- lev
      list(X = X, meta = list(kind = "vm", var = comp$var, relmat = comp$relmat,
                              levels = lev, pc = rot$PC))
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

#' Build one ETA entry (design + BGLR model) for a whole term
#' @keywords internal
.term_eta <- function(term, data, relmat) {
  comp_designs <- lapply(term$components, .component_design,
                         data = data, relmat = relmat, role = term$role)
  X <- comp_designs[[1]]$X
  if (length(comp_designs) > 1L) {
    for (i in 2:length(comp_designs)) X <- .khatri_rao_rows(X, comp_designs[[i]]$X)
  }
  model <- if (identical(term$role, "fixed")) "FIXED" else "BRR"
  entry <- list(X = X, model = model)
  if (model == "BRR") entry$saveEffects <- TRUE
  # Factor-analytic genetic covariance (multi-trait only)
  fa_comp <- Filter(function(cd) identical(cd$meta$kind, "fa"), comp_designs)
  if (length(fa_comp)) entry$Cov <- list(type = "FA", nF = fa_comp[[1]]$meta$nfac)
  list(eta = entry,
       meta = list(label = term$label, role = term$role, model = model,
                   components = lapply(comp_designs, `[[`, "meta")))
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
      }
    }
  }
  if (!all(keep)) {
    warning(sum(!keep), " row(s) dropped: their vm() levels are absent from the ",
            "relationship matrix.", call. = FALSE)
  }
  droplevels(data[keep, , drop = FALSE])
}

#' Assemble the full BGLR ETA list from parsed model terms
#'
#' @param parsed A `breedRB_terms` object from [parse_model()].
#' @param data Model data frame.
#' @param relmat Named list of relationship / covariance matrices for `vm()`.
#' @return A list with `ETA` (the BGLR ETA list) and `meta` (per-term metadata,
#'   named by a safe key) plus `group` (residual grouping vector or `NULL`).
#' @keywords internal
build_eta <- function(parsed, data, relmat = list()) {
  all_terms <- c(parsed$fixed, parsed$random)
  built <- lapply(all_terms, .term_eta, data = data, relmat = relmat)

  keys <- make.unique(gsub("[^A-Za-z0-9]", "_", vapply(built, function(b) b$meta$label,
                                                       character(1))))
  ETA  <- stats::setNames(lapply(built, `[[`, "eta"),  keys)
  meta <- stats::setNames(lapply(built, `[[`, "meta"), keys)

  group <- NULL
  if (isTRUE(parsed$residual$hetero)) {
    group <- as.integer(factor(data[[parsed$residual$group]]))
  }
  list(ETA = ETA, meta = meta, group = group)
}
