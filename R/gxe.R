# ---------------------------------------------------------------------------
# gxe(): summarise a random-regression (reaction-norm) fit into per-genotype
# GxE descriptors -- broad adaptability (intercept / overall level),
# responsiveness (slope across the gradient) and stability (CVge, the
# coefficient of variation of the genotype's curve across the gradient).
# ---------------------------------------------------------------------------

#' Genotype-by-environment summary of a random-regression fit
#'
#' Reduces each genotype's fitted reaction norm \eqn{v_i(x)} across the
#' environmental gradient to three interpretable genotype-by-environment
#' descriptors:
#'
#' * **Broad adaptability** — the genotype intercept, i.e. its overall
#'   performance level across the gradient (mean of its curve). Higher = better
#'   general performance.
#' * **Responsiveness** — the slope of the genotype's curve along the gradient
#'   (the linear reaction-norm slope; for a `leg()` order \eqn{> 1} the average
#'   slope from a straight-line fit to the curve). Larger magnitude = more
#'   sensitive to the environment; near zero = flat / broadly adapted.
#' * **Stability (`cv_ge`)** — the coefficient of variation of the genotype's
#'   curve across the gradient, \eqn{100 \times \mathrm{sd}(v_i(x)) /
#'   |\mathrm{mean}(v_i(x))|}. **Lower `cv_ge` means a more stable genotype**
#'   (flatter curve relative to its level).
#'
#' The curves are the posterior-mean per-genotype reaction norms from
#' [reaction_norm()]. `cv_ge` is only meaningful on the phenotype scale, so
#' `add_fixed_reg = TRUE` (the default) is recommended.
#'
#' @param fit A `breedRB_fit` (single-trait) containing a random regression.
#' @param term The interaction term identifier, e.g. `"gen:leg(x)"` or
#'   `"mrk(gen, M):leg(x, 1)"`; its matching intercept term is detected
#'   automatically (as in [reaction_norm()]).
#' @param n_grid Integer; number of gradient points used to evaluate each curve
#'   (default 100).
#' @param add_fixed_reg Logical (default `TRUE`); evaluate the curves on the
#'   phenotype scale (`mu` + fixed population regression added) so `cv_ge` has a
#'   meaningful mean in the denominator.
#' @param higher Logical (default `TRUE`); if `TRUE` a higher value is better
#'   (adaptability rank 1 = highest); set `FALSE` to rank the smallest as best.
#'
#' @return A data frame of class `breedRB_gxe`, one row per genotype ordered by
#'   broad adaptability, with columns `id`, `adaptability` (overall level),
#'   `responsiveness` (gradient slope), `cv_ge` (coefficient of variation, \%),
#'   `adaptability_rank` and `stability_rank` (rank 1 = most stable, i.e.
#'   smallest `cv_ge`).
#' @examples
#' \donttest{
#' fit <- bbglr(yield ~ leg(x), random = ~ gen + gen:leg(x) + env:rep, data = dat)
#' gxe(fit, term = "gen:leg(x)")
#' }
#' @seealso [reaction_norm()] for the curves and [rr_gradient()] for the
#'   across-gradient (co)variance / heritability / selection analytics.
#' @export
gxe <- function(fit, term, n_grid = 100L, add_fixed_reg = TRUE, higher = TRUE) {
  stopifnot(inherits(fit, "breedRB_fit"))
  if (isTRUE(fit$response$multitrait)) {
    stop("gxe() currently supports single-trait fits.", call. = FALSE)
  }

  rn <- reaction_norm(fit, term = term, add_fixed_reg = add_fixed_reg,
                      leg_basis = TRUE, n_grid = n_grid, plot = FALSE)
  if (!isTRUE(add_fixed_reg)) {
    warning("cv_ge is computed on the genotype-deviation scale (add_fixed_reg = ",
            "FALSE); its denominator (curve mean) is near zero, so the ",
            "coefficient of variation may be unstable. Use add_fixed_reg = TRUE ",
            "for a phenotype-scale stability measure.", call. = FALSE)
  }

  sp <- split(rn, rn$id)
  g  <- data.frame(
    id             = names(sp),
    adaptability   = vapply(sp, function(d) mean(d$value), numeric(1)),
    responsiveness = vapply(sp, function(d) {
      stats::cov(d$value, d$gradient) / stats::var(d$gradient)
    }, numeric(1)),
    cv_ge          = vapply(sp, function(d) {
      m <- mean(d$value)
      100 * stats::sd(d$value) / abs(m)
    }, numeric(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )

  g$adaptability_rank <- rank(if (isTRUE(higher)) -g$adaptability else g$adaptability,
                              ties.method = "min")
  g$stability_rank    <- rank(g$cv_ge, ties.method = "min")   # 1 = most stable

  g <- g[order(g$adaptability_rank), , drop = FALSE]
  rownames(g) <- NULL
  attr(g, "term")   <- fit$meta[[.resolve_term(fit, term)]]$label
  attr(g, "higher") <- higher
  class(g) <- c("breedRB_gxe", "data.frame")
  g
}

#' @export
print.breedRB_gxe <- function(x, ...) {
  cat("<breedRB_gxe>  term:", attr(x, "term"),
      sprintf(" (%s is better)\n", if (isTRUE(attr(x, "higher"))) "higher" else "lower"))
  cat("  adaptability = overall level | responsiveness = gradient slope",
      "| cv_ge = stability (lower = more stable)\n")
  print(as.data.frame(x), row.names = FALSE, ...)
  invisible(x)
}
