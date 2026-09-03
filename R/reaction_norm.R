# ---------------------------------------------------------------------------
# reaction_norm(): evaluate per-genotype reaction-norm curves across the
# environmental gradient for a genomic random-regression fit.
# ---------------------------------------------------------------------------

# aes() columns of the long data frame built inside reaction_norm()
utils::globalVariables(c("gradient", "value", "id"))

#' Find the genomic main-effect (intercept) term matching a random-regression term
#'
#' Given the genomic component `gc` of an interaction term, locates the random
#' term whose single component is the same genomic effect (same variable and
#' matrix) — the reaction-norm intercept.
#' @return The matching term key, or `NULL` if none is found.
#' @keywords internal
.rr_intercept_key <- function(fit, gc) {
  for (k in .random_keys(fit)) {
    comps <- fit$meta[[k]]$components
    if (length(comps) != 1L) next
    c1 <- comps[[1]]
    if (c1$kind %in% c("vm", "mrk") &&
        identical(c1$var, gc$var) && identical(c1$relmat, gc$relmat)) {
      return(k)
    }
  }
  NULL
}

#' Fixed-effect population regression curve evaluated on a Legendre grid
#'
#' Returns the fitted fixed `leg()` regression for gradient variable `legvar`
#' evaluated at `grid` (on the standardized \[-1, 1] domain), i.e. the population
#' reaction norm relative to `mu`. Zero if the model has no such fixed term.
#' @keywords internal
.fixed_leg_curve <- function(fit, legvar, grid) {
  for (k in names(fit$meta)) {
    m <- fit$meta[[k]]
    if (!identical(m$role, "fixed")) next
    comps <- m$components
    if (length(comps) == 1L && identical(comps[[1]]$kind, "leg") &&
        identical(comps[[1]]$var, legvar)) {
      qf  <- comps[[1]]$order
      sol <- solution(fit, term = k, type = "fixed")
      cf  <- stats::setNames(sol$solution, sol$effect)[paste0("deg", seq_len(qf))]
      Bf  <- legendre_basis(grid, order = qf, orthonormal = TRUE)[, -1, drop = FALSE]
      return(as.numeric(Bf %*% cf))
    }
  }
  rep(0, length(grid))
}

#' Posterior-mean intercept `mu`, dropping the thinned burn-in
#' @keywords internal
.mu_mean <- function(fit) {
  nburn <- floor(fit$control$burnIn / fit$control$thin)
  mean(unlist(lapply(fit$paths, function(p) {
    f <- paste0(p, "mu.dat")
    if (!file.exists(f)) return(numeric(0))
    mu <- scan(f, quiet = TRUE)
    if (length(mu) > nburn) mu[(nburn + 1L):length(mu)] else mu
  })), na.rm = TRUE)
}

#' Genotype reaction-norm curves across the environmental gradient
#'
#' For a genomic random-regression (reaction-norm) fit
#' (`random = ~ mrk(gen, M) + mrk(gen, M):leg(x, q)`), evaluates each genotype's
#' fitted curve over the environmental gradient,
#' \deqn{v_i(x) = \hat a_i + \sum_{j=1}^{q} \hat b_{i,j}\, L_j(x),}
#' where the intercept \eqn{\hat a_i} comes from the genomic main-effect term and
#' the reaction-norm coefficients \eqn{\hat b_{i,j}} from the interaction `term`
#' (`deg1` = linear slope, `deg2` = quadratic, ...), and \eqn{L_j} is the same
#' orthonormal Legendre basis used in fitting, on the standardized domain
#' \eqn{[-1, 1]}.
#'
#' @param fit A `breedRB_fit` (single-trait) containing a genomic random
#'   regression.
#' @param term The interaction term identifier, e.g. `"mrk(gen, M):leg(x, 1)"`
#'   (label as written in the formula, or the internal key). Its matching
#'   intercept term (`mrk(gen, M)`) is detected automatically.
#' @param type Optional; validated against the term's role (must be `"random"`
#'   for a reaction norm).
#' @param add_fixed_reg Logical (default `TRUE`). Add the fixed population
#'   regression (`mu` + the fixed `leg()` curve on the same gradient) so the
#'   curves sit on the phenotype scale; otherwise curves are genotype deviations
#'   around zero.
#' @param plot Logical (default `TRUE`). Draw and print a \pkg{ggplot2} figure of
#'   value versus gradient, one line per genotype. The plot is also returned as
#'   attribute `"plot"`.
#' @param leg_basis Logical (default `TRUE`). If `TRUE` the gradient axis is the
#'   standardized Legendre domain \eqn{[-1, 1]}; if `FALSE` it is back-transformed
#'   to the original covariate scale using the range stored at fitting.
#' @param n_grid Integer; number of gradient points (default 100).
#'
#' @return A long-format data frame with `id` (genotype), `gradient` (the
#'   evaluation point, on the chosen scale) and `value` (fitted reaction-norm
#'   value). The gradient grid is attached as attribute `"grid"` and, when
#'   `plot = TRUE`, the \pkg{ggplot2} object as attribute `"plot"`.
#' @examples
#' \donttest{
#' fit <- bbglr(y ~ 1 + leg(x, 1),
#'              random = ~ mrk(gen, M) + mrk(gen, M):leg(x, 1),
#'              data = dat, relmat = list(M = M))
#' rn <- reaction_norm(fit, term = "mrk(gen, M):leg(x, 1)", type = "random")
#' }
#' @seealso [solution()] for the underlying intercept and slope coefficients.
#' @export
reaction_norm <- function(fit, term, type = NULL, add_fixed_reg = TRUE,
                          plot = TRUE, leg_basis = TRUE, n_grid = 100L) {
  stopifnot(inherits(fit, "breedRB_fit"))
  if (isTRUE(fit$response$multitrait)) {
    stop("reaction_norm() currently supports single-trait fits.", call. = FALSE)
  }

  key  <- .resolve_term(fit, term)
  meta <- fit$meta[[key]]
  rr   <- .geno_leg_parts(meta)
  if (is.null(rr)) {
    stop("Term '", term, "' is not a genomic random-regression interaction ",
         "(expected a mrk()/vm() x leg() term such as \"mrk(gen, M):leg(x, 1)\").",
         call. = FALSE)
  }
  if (!is.null(type)) {
    type <- match.arg(type, c("random", "fixed"))
    if (!identical(type, meta$role)) {
      stop("Term '", term, "' is a ", meta$role, " effect, but type = \"", type,
           "\" was requested.", call. = FALSE)
    }
  }

  gc     <- rr$gc; q <- rr$q
  legvar <- meta$components[[rr$lpos]]$var
  rng    <- meta$components[[rr$lpos]]$range

  # Intercept (genomic main effect) and per-degree slopes, both per genotype.
  ikey <- .rr_intercept_key(fit, gc)
  if (is.null(ikey)) {
    stop("No genomic main-effect (intercept) term matching '", term,
         "' was found; a reaction norm needs both mrk(gen, M) and ",
         "mrk(gen, M):leg(x, q) in the model.", call. = FALSE)
  }
  ii <- solution(fit, term = ikey, type = "random")
  ss <- solution(fit, term = key,  type = "random")

  ids <- ii$effect
  int <- stats::setNames(ii$solution, ii$effect)[ids]

  # Reshape the per-degree slope solutions into an [nGen x q] coefficient matrix.
  S <- matrix(0, nrow = length(ids), ncol = q, dimnames = list(ids, NULL))
  slp <- stats::setNames(ss$solution, ss$effect)
  for (j in seq_len(q)) {
    nm <- .rr_effect_names(ids, q, j)
    S[, j] <- slp[nm]
  }

  # Legendre grid on the standardized [-1, 1] domain (as used in fitting).
  gstd <- seq(-1, 1, length.out = n_grid)
  B    <- legendre_basis(gstd, order = q, orthonormal = TRUE)[, -1, drop = FALSE]

  # Per-genotype curves: intercept + basis %*% slopes  ([nGen x nGrid]).
  curves <- outer(int, rep(1, n_grid)) + S %*% t(B)

  if (isTRUE(add_fixed_reg)) {
    pop <- .mu_mean(fit) + .fixed_leg_curve(fit, legvar, gstd)
    curves <- sweep(curves, 2L, pop, "+")
  }

  # Gradient axis: standardized domain, or back-transformed to the covariate scale.
  gx <- if (isTRUE(leg_basis)) gstd
        else rng[1] + (gstd + 1) / 2 * diff(rng)

  out <- data.frame(
    id       = rep(ids, times = n_grid),
    gradient = rep(gx,  each  = length(ids)),
    value    = as.numeric(curves),
    row.names = NULL
  )
  attr(out, "grid") <- gx

  if (isTRUE(plot)) {
    xlab <- if (isTRUE(leg_basis)) paste0("Legendre gradient  [-1, 1]  (", legvar, ")")
            else legvar
    p <- ggplot2::ggplot(out, ggplot2::aes(x = gradient, y = value, group = id)) +
      ggplot2::geom_line(alpha = 0.4, linewidth = 0.3, colour = "steelblue") +
      ggplot2::labs(x = xlab,
                    y = if (isTRUE(add_fixed_reg)) "Reaction norm (phenotype scale)"
                        else "Reaction norm (genotype deviation)",
                    title = paste0("Reaction norms: ", meta$label)) +
      ggplot2::theme_minimal()
    print(p)
    attr(out, "plot") <- p
  }
  out
}
