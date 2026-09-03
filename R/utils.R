#' Orthogonal / orthonormal Legendre polynomial basis on \[-1, 1]
#'
#' Constructs a Legendre polynomial basis evaluated at a numeric vector scaled to
#' \[-1, 1], up to a specified order, using the three-term recurrence. Optionally
#' returns the orthonormal basis.
#'
#' @param x Numeric vector (already scaled to \[-1, 1]) where the basis is evaluated.
#' @param order Integer. Maximum polynomial degree. Returns columns for degrees 0..order.
#' @param orthonormal Logical. If `TRUE`, applies orthonormal scaling to each degree.
#'
#' @return A numeric matrix `length(x) x (order + 1)`; column `j` is the degree
#'   `j - 1` polynomial.
#'
#' @examples
#' legendre_basis(seq(-1, 1, length.out = 5), order = 3)
#' @export
legendre_basis <- function(x, order, orthonormal = FALSE) {
  n <- length(x)
  B <- matrix(0, nrow = n, ncol = order + 1)
  B[, 1] <- 1                                   # P0(x) = 1
  if (order >= 1) B[, 2] <- x                   # P1(x) = x
  if (order >= 2) {
    for (k in 1:(order - 1)) {
      # P_{k+1}(x) = ((2k+1) x P_k(x) - k P_{k-1}(x)) / (k+1)
      B[, k + 2] <- ((2 * k + 1) * x * B[, k + 1] - k * B[, k]) / (k + 1)
    }
  }
  if (orthonormal) {
    for (d in 0:order) B[, d + 1] <- sqrt((2 * d + 1) / 2) * B[, d + 1]
  }
  B
}

#' Scale a numeric vector to the interval \[-1, 1]
#'
#' @param x Numeric vector.
#' @param rng Optional length-2 numeric giving the range to scale against
#'   (defaults to `range(x)`). Useful to reuse a training range at prediction time.
#' @return A numeric vector on \[-1, 1] (all zeros if the range has zero width).
#' @keywords internal
.scale_unit <- function(x, rng = range(x, na.rm = TRUE)) {
  if (diff(rng) > 0) 2 * (x - rng[1]) / diff(rng) - 1 else rep(0, length(x))
}

#' Build a factor incidence matrix with clean level names
#'
#' @param f A vector coercible to factor.
#' @param center Logical; if `TRUE` columns are centered (mean 0), as used for
#'   BGLR random main effects.
#' @return A numeric incidence matrix with one column per level; column names are
#'   the factor levels.
#' @keywords internal
.incidence <- function(f, center = FALSE) {
  f <- factor(f)
  Z <- model.matrix(~ -1 + f)
  colnames(Z) <- levels(f)
  if (center) Z <- scale(Z, center = TRUE, scale = FALSE)
  Z
}

#' Expand a level-indexed matrix to observation rows by gathering
#'
#' Equivalent to `.incidence(f) %*% M` (the 0/1 factor design of `f` times a
#' matrix `M` whose rows are indexed by `levels`), but computed as a pure
#' row-gather: no dense `N x nLevel` incidence matrix is built and no matrix
#' product is formed. Row `i` of the result is `M[level_of(f_i), ]`. This is
#' numerically identical to the incidence product (an incidence matrix only
#' selects rows) but avoids the `O(N x nLevel)` allocation and the
#' `O(N x nLevel x ncol(M))` multiply — the main memory cost under replication.
#'
#' @param M Matrix with one row per level (row order matching `levels`).
#' @param f Vector of observations, coercible to a factor with these `levels`.
#' @param levels Character vector of levels, in the same order as `rows(M)`.
#' @return A matrix with one row per element of `f`.
#' @keywords internal
.expand_rows <- function(M, f, levels) {
  idx <- as.integer(factor(f, levels = levels))
  M[idx, , drop = FALSE]
}

#' Clean a marker matrix before building a genomic model
#'
#' Mean-imputes missing values (per marker/column) and drops markers with
#' (near-)zero variance, so the downstream GBLUP / RR-BLUP construction never
#' sees `NA`s or constant columns. Markers that are entirely missing become
#' non-finite after imputation and are dropped alongside the constant ones.
#'
#' @param M Numeric marker matrix (genotypes in rows, markers in columns).
#' @param var_tol Markers with variance below this are treated as constant and
#'   removed.
#' @param var_name Factor name, used only for informative messages.
#' @return The cleaned marker matrix (possibly with fewer columns).
#' @keywords internal
.clean_markers <- function(M, var_tol = 1e-8, var_name = "gen") {
  M <- as.matrix(M)
  storage.mode(M) <- "double"

  na_idx <- is.na(M)
  n_na   <- sum(na_idx)
  if (n_na) {                                   # mean-impute each marker's missing calls
    cm <- colMeans(M, na.rm = TRUE)
    M[na_idx] <- cm[col(M)[na_idx]]             # all-NA columns -> NaN, dropped below
  }

  v    <- apply(M, 2L, stats::var)
  drop <- !is.finite(v) | v < var_tol           # constant / near-constant / all-NA markers
  n_drop <- sum(drop)
  if (n_drop) M <- M[, !drop, drop = FALSE]

  if (!ncol(M)) {
    stop("mrk(", var_name, "): no markers remain after imputing missing values and ",
         "dropping near-constant markers.", call. = FALSE)
  }
  if (n_na || n_drop) {
    message("mrk(", var_name, "): imputed ", n_na, " missing value(s) (marker means); ",
            "dropped ", n_drop, " near-constant marker(s).")
  }
  M
}

#' Number of leading components explaining a target fraction of variance
#'
#' Given eigenvalues in **descending** order, returns the smallest `k` such that
#' the top `k` account for at least `exp_var` of their total. Used to pick a
#' self-scaling low-rank GBLUP truncation from the genomic relationship's
#' spectrum.
#'
#' @param evals Numeric eigenvalues, largest first.
#' @param exp_var Target cumulative fraction in `(0, 1]`.
#' @return An integer rank `>= 1`.
#' @keywords internal
.rank_for_expvar <- function(evals, exp_var) {
  evals <- pmax(as.numeric(evals), 0)
  tot   <- sum(evals)
  if (!is.finite(tot) || tot <= 0) return(1L)
  k <- which(cumsum(evals) / tot >= exp_var)[1]
  if (is.na(k)) k <- length(evals)
  max(1L, as.integer(k))
}

#' Eigen (principal-component) rotation of a relationship matrix for GBLUP
#'
#' Reproduces the eigen-decomposition trick used to fit GBLUP through a Bayesian
#' ridge regression in \pkg{BGLR}: `K = PC PC'`, so a `BRR` model on the
#' incidence-times-`PC` design has genotype effects recoverable as
#' `tcrossprod(PC, B)`.
#'
#' Numerically-null directions (eigenvalue below `tol`) are dropped: they carry
#' no genetic variance, so an all-zero `PC` column would only enlarge the design
#' without changing the fit. If `exp_var` is supplied, only the leading
#' eigenvectors explaining that fraction of variance are kept (a low-rank /
#' principal-component GBLUP approximation).
#'
#' @param K Symmetric relationship / covariance matrix (row/col names = levels).
#' @param tol Numeric tolerance below which eigenvalues are treated as zero.
#' @param exp_var Optional explained-variance target in `(0, 1]`: keep the
#'   smallest set of leading eigenvectors accounting for at least this fraction of
#'   the total variance. `NA` (default) keeps every non-null direction (exact).
#' @return A list with `PC` (the rotation, `nLevel x rank`), `vectors` and
#'   `values` (the retained eigenpairs) and `levels` (row/column names of `K`).
#' @keywords internal
.pc_rotation <- function(K, tol = 1e-8, exp_var = NA_real_) {
  evd  <- eigen(K, symmetric = TRUE)
  vals <- evd$values
  vecs <- evd$vectors
  if (!is.na(exp_var) && exp_var < 1) {                # variance-based truncation
    k <- .rank_for_expvar(vals, exp_var)
    vals <- vals[seq_len(k)]; vecs <- vecs[, seq_len(k), drop = FALSE]
  }
  vals[vals < tol] <- 0
  drop <- vals <= 0                                    # exact: null directions carry no variance
  if (all(drop)) drop[which.max(vals)] <- FALSE        # always keep at least one column
  vals <- vals[!drop]; vecs <- vecs[, !drop, drop = FALSE]
  PC <- sweep(vecs, MARGIN = 2, STATS = sqrt(vals), FUN = "*")
  rownames(PC) <- rownames(K); rownames(vecs) <- rownames(K)
  list(PC = PC, vectors = vecs, values = vals, levels = rownames(K))
}

#' Leading eigenpairs of `Mc Mc'` explaining a variance target, without the tail
#'
#' Finds the smallest set of top principal components of the genomic relationship
#' that explain at least `exp_var` of its total variance, **without decomposing
#' the full spectrum**. The trick: the total variance is the trace of `Mc Mc'`,
#' i.e. `sum(rowSums(Mc^2))` — available in `O(np)` with no decomposition — so it
#' is the exact denominator up front. Only the top eigenpairs are then computed,
#' grown a block at a time with `RSpectra::svds` (which never forms `Mc Mc'`),
#' stopping as soon as the running cumulative share crosses the target.
#'
#' Cheap when the relationship has strong low-rank structure (a handful of PCs
#' dominate, as in real genomic panels): only a small `k << n` is ever computed.
#' If the target would require more than about half the spectrum (little
#' structure), the search gives up and returns `NULL` so the caller can fall back
#' to a single full eigendecomposition — which is cheaper than svds at that point.
#'
#' @param Mc Column-centred marker matrix (`n` genotypes x `p` markers).
#' @param exp_var Target cumulative variance fraction in `(0, 1)`.
#' @param total Total variance (`sum(rowSums(Mc^2))`); pass it if already known.
#' @param k0 Initial block size; grown by `growth` each round.
#' @param growth Multiplicative block-growth factor (`> 1`).
#' @param k_max Largest `k` to attempt before giving up (defaults to half the
#'   spectrum). Beyond this a full eigendecomposition is cheaper.
#' @return A list with `vectors` (top eigenvectors of `Mc Mc'`) and `evals` (their
#'   eigenvalues), truncated to the target; or `NULL` to signal "fall back to a
#'   full eigendecomposition".
#' @keywords internal
.svds_expvar <- function(Mc, exp_var, total = sum(rowSums(Mc^2)),
                         k0 = NULL, growth = 2, k_max = NULL) {
  if (!requireNamespace("RSpectra", quietly = TRUE)) return(NULL)
  n <- nrow(Mc)
  if (!is.finite(total) || total <= 0) return(NULL)
  if (is.null(k_max)) k_max <- max(1L, floor((n - 1L) / 2L))  # svds no longer pays past ~half
  if (is.null(k0))    k0    <- min(k_max, max(16L, as.integer(ceiling(0.05 * n))))
  k <- min(as.integer(k0), n - 1L, k_max)

  repeat {
    sv     <- RSpectra::svds(Mc, k)                   # top-k SVD of Mc; Mc Mc' never formed
    lambda <- sv$d^2                                  # eigenvalues of Mc Mc', descending
    frac   <- cumsum(lambda) / total
    if (frac[length(frac)] >= exp_var) {              # target reached within these top-k
      kstar <- which(frac >= exp_var)[1]
      return(list(vectors = sv$u[, seq_len(kstar), drop = FALSE],
                  evals   = lambda[seq_len(kstar)]))
    }
    if (k >= k_max || k >= n - 1L) return(NULL)       # too many needed -> caller does full eigen
    k <- min(as.integer(ceiling(k * growth)), k_max, n - 1L)
  }
}

#' GBLUP principal-component rotation built directly from a marker matrix
#'
#' Computes the same rotation as `.pc_rotation()` applied to the VanRaden genomic
#' relationship `G = Mc Mc'/cc + ridge`, but works from the centred markers `Mc`
#' so it can (a) avoid ever forming the `n x n` matrix `G` and (b) return the
#' eigenpairs of `Mc Mc'` for later reuse by the marker back-solve.
#'
#' When an `exp_var` target is set and \pkg{RSpectra} is available, the required
#' leading components are grown incrementally ([.svds_expvar()]), so neither `G`
#' nor `Mc Mc'` is materialised for a low-rank panel. Otherwise `Mc Mc'` is formed
#' and eigendecomposed. The tiny diagonal ridge (`1e-6 * mean(diag(G))`) that
#' keeps `G` positive-definite is applied through the eigenvalues.
#'
#' @param Mc Column-centred marker matrix (`n` genotypes x `p` markers), row
#'   names = genotype IDs.
#' @param cc VanRaden scaling `2 * sum(p_j q_j)`.
#' @param exp_var Optional explained-variance target in `(0, 1]`: keep the
#'   smallest set of leading components accounting for at least this fraction of
#'   the total genomic variance. When \pkg{RSpectra} is available this is found by
#'   growing the top eigenpairs incrementally ([.svds_expvar()]) — the total
#'   variance is the trace of `Mc Mc'`, so only the leading components are computed
#'   and the `n x n` matrix is never formed for a low-rank panel. It falls back to
#'   a single full eigendecomposition of `Mc Mc'` when \pkg{RSpectra} is absent,
#'   `n` is small, or the target needs most of the spectrum. `NA` (default) keeps
#'   every non-null component.
#' @param tol Eigenvalue tolerance for dropping null directions.
#' @return A list with `PC` (rotation `n x r`), `vectors` (eigenvectors of
#'   `Mc Mc'`), `evals_MM` (its eigenvalues) and `mean_diag_MM`
#'   (`mean(diag(Mc Mc'))`, needed to reproduce the back-solve ridge).
#' @keywords internal
.gblup_rotation <- function(Mc, cc, exp_var = NA_real_, tol = 1e-8) {
  n <- nrow(Mc)
  mean_diag_MM <- mean(rowSums(Mc^2))                  # = mean(diag(Mc Mc')), cheap
  total_var    <- mean_diag_MM * n                     # = sum(rowSums(Mc^2)) = sum of all eigvals
  ridge <- 1e-6 * mean_diag_MM / cc                    # == 1e-6 * mean(diag(G))
  want_expvar <- !is.na(exp_var) && exp_var < 1
  V <- NULL; lambda <- NULL

  if (want_expvar) {
    # Reach exp_var by growing the top eigenpairs, using the (free) trace as the
    # denominator -> never touches the tail / forms Mc Mc' for a low-rank panel.
    hit <- .svds_expvar(Mc, exp_var, total = total_var)
    if (!is.null(hit)) { V <- hit$vectors; lambda <- hit$evals }
  }

  if (is.null(V)) {                                    # full-spectrum fallback
    MM <- tcrossprod(Mc)                               # n x n (the only path that forms it)
    ev <- eigen(MM, symmetric = TRUE)
    V  <- ev$vectors
    lambda <- pmax(ev$values, 0)
    if (want_expvar) {
      k <- .rank_for_expvar(lambda, exp_var)           # variance-based truncation
      V <- V[, seq_len(k), drop = FALSE]; lambda <- lambda[seq_len(k)]
    }
  }
  gamma <- lambda / cc + ridge                         # eigenvalues of G = Mc Mc'/cc + ridge I
  keep  <- gamma > tol
  if (!any(keep)) keep[which.max(gamma)] <- TRUE
  V <- V[, keep, drop = FALSE]; lambda <- lambda[keep]; gamma <- gamma[keep]
  rownames(V) <- rownames(Mc)
  PC <- sweep(V, MARGIN = 2, STATS = sqrt(gamma), FUN = "*")
  rownames(PC) <- rownames(Mc)
  list(PC = PC, vectors = V, evals_MM = lambda, mean_diag_MM = mean_diag_MM)
}
