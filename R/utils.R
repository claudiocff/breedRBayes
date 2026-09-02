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

#' Eigen (principal-component) rotation of a relationship matrix for GBLUP
#'
#' Reproduces the eigen-decomposition trick used to fit GBLUP through a Bayesian
#' ridge regression in \pkg{BGLR}: `K = PC PC'`, so a `BRR` model on the
#' incidence-times-`PC` design has genotype effects recoverable as
#' `tcrossprod(PC, B)`.
#'
#' @param K Symmetric relationship / covariance matrix (row/col names = levels).
#' @param tol Numeric tolerance below which eigenvalues are floored to 0.
#' @return A list with `PC` (the rotation, `nLevel x nLevel`) and `levels`
#'   (row/column names of `K`).
#' @keywords internal
.pc_rotation <- function(K, tol = 1e-8) {
  evd <- eigen(K, symmetric = TRUE)
  vals <- evd$values
  vals[vals < tol] <- 0
  PC <- sweep(evd$vectors, MARGIN = 2, STATS = sqrt(vals), FUN = "*")
  rownames(PC) <- rownames(K)
  list(PC = PC, levels = rownames(K))
}
