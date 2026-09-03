# ---------------------------------------------------------------------------
# rr_gradient(): reaction-norm analytics across the environmental gradient for a
# random-regression fit -- genetic correlation surface (phi K phi'), heritability
# over the gradient, per-genotype reliability, and top-fraction selection
# probabilities, all with heat-map / ribbon visualisations.
# ---------------------------------------------------------------------------

# aes() columns of the long data frames built inside rr_gradient()
utils::globalVariables(c("x1", "x2", "gcor", "h2", "lower", "upper",
                         "reliability", "prob"))

#' Per-draw residual variance, burn-in-dropped to align with the effect draws
#'
#' The variance-component traces on disk (`varE.dat`) still carry the thinned
#' burn-in, whereas the effect draws (`readBinMat`) do not. Dropping
#' `floor(burnIn/thin)` samples per chain aligns `varE` draw-for-draw with the
#' pooled effect draws (same chain order, same post-burn-in samples).
#' @keywords internal
.varE_draws <- function(fit) {
  nburn <- floor(fit$control$burnIn / fit$control$thin)
  unlist(lapply(fit$paths, function(p) {
    v <- scan(paste0(p, "varE.dat"), quiet = TRUE)
    if (length(v) > nburn) v[(nburn + 1L):length(v)] else v
  }), use.names = FALSE)
}

#' Per-draw reaction-norm coefficients (intercept + per-degree slopes)
#'
#' Pulls the pooled posterior draws of the grouping main effect (the intercept
#' \eqn{a_i}) and of the `leg()` interaction (`term`, the per-degree slopes
#' \eqn{b_{i,j}}) and reshapes them into per-genotype coefficient draw matrices,
#' all aligned to a common genotype order. Works for both genomic
#' (`mrk(gen, M):leg(x, q)`) and plain-factor (`gen:leg(x)`) random regressions.
#' @return A list with `ids`, `q`, `legvar`, `range`, `genomic`, the intercept
#'   draws `A` (`[nDraws x nGeno]`) and a length-`q` list `B` of per-degree slope
#'   draws (each `[nDraws x nGeno]`, columns in the `ids` order).
#' @keywords internal
.rr_coef_draws <- function(fit, term) {
  key  <- .resolve_term(fit, term)
  meta <- fit$meta[[key]]
  rr   <- .rr_leg_parts(meta)
  if (is.null(rr)) {
    stop("Term '", term, "' is not a random-regression interaction ",
         "(expected a grouping x leg() term such as \"gen:leg(x)\" or ",
         "\"mrk(gen, M):leg(x, 1)\").", call. = FALSE)
  }
  ikey <- .rr_intercept_key(fit, rr$gc)
  if (is.null(ikey)) {
    stop("No main-effect (intercept) term matching '", term, "' was found; a ",
         "reaction norm needs both the main effect (e.g. ", rr$gc$var, ") and ",
         "its leg() interaction in the model.", call. = FALSE)
  }

  A <- attr(solution(fit, term = ikey, type = "random"), "draws")   # [nDraws x nGeno]
  ids <- colnames(A)

  S    <- attr(solution(fit, term = key, type = "random"), "draws")  # [nDraws x (nGeno*q)]
  nm   <- colnames(S)
  hasd <- grepl(":deg[0-9]+$", nm)
  sid  <- ifelse(hasd, sub(":deg[0-9]+$", "", nm), nm)
  sdeg <- ifelse(hasd, as.integer(sub(".*:deg([0-9]+)$", "\\1", nm)), 1L)

  B <- lapply(seq_len(rr$q), function(j) {
    mat  <- matrix(0, nrow = nrow(S), ncol = length(ids), dimnames = list(NULL, ids))
    take <- which(sdeg == j)
    ord  <- match(ids, sid[take])
    ok   <- !is.na(ord)
    mat[, ok] <- S[, take[ord[ok]], drop = FALSE]
    mat
  })

  list(ids = ids, q = rr$q, legvar = rr$legvar, range = rr$range,
       genomic = rr$genomic, A = A[, ids, drop = FALSE], B = B)
}

#' Posterior summaries of the reaction-norm coefficient (co)variances
#'
#' For every random-regression term in the fit (a grouping x `leg()` interaction
#' with a matching intercept term), computes the posterior distribution of the
#' entries of the across-genotype coefficient (co)variance matrix \eqn{K}: the
#' per-coefficient **variances** `var(...)` (intercept and each Legendre degree —
#' the diagonal) and the intercept--slope / slope--slope **covariances**
#' `cov(...)` (the off-diagonal). Each entry is computed per MCMC draw as `cov()`
#' of the per-genotype coefficients in that draw (the same realised (co)variance
#' used by [rr_gradient()]) and then summarised. This exposes the individual
#' per-degree variances even though \pkg{BGLR} fits a **single shared** variance
#' component for the whole `leg()` interaction term, and gives a genuine
#' intercept--slope covariance even though the terms have independent priors.
#' @return A data frame with the same columns as [varcomp()]
#'   (`term`, `mean`, `median`, `sd`, `lower`, `upper`) — the `var(...)` rows
#'   first, then the `cov(...)` rows — or `NULL` if the model has no random
#'   regression.
#' @keywords internal
.rr_cov_summaries <- function(fit, prob = 0.95) {
  rows <- list()
  for (key in .random_keys(fit)) {
    rr <- .rr_leg_parts(fit$meta[[key]])
    if (is.null(rr)) next
    ikey <- .rr_intercept_key(fit, rr$gc)
    if (is.null(ikey)) next
    cd <- tryCatch(.rr_coef_draws(fit, key), error = function(e) NULL)
    if (is.null(cd)) next
    q <- cd$q; G <- length(cd$ids); nD <- nrow(cd$A); p1 <- q + 1L
    coef_lab <- c(fit$meta[[ikey]]$label,
                  if (q == 1L) fit$meta[[key]]$label
                  else paste0(fit$meta[[key]]$label, ":deg", seq_len(q)))
    pairs  <- if (p1 >= 2L) utils::combn(p1, 2L) else matrix(integer(0), 2L, 0L)
    diag_d <- matrix(0, nD, p1)                       # per-draw variances (diagonal of K)
    off_d  <- matrix(0, nD, ncol(pairs))              # per-draw covariances (off-diagonal)
    for (t in seq_len(nD)) {
      Ct <- cbind(cd$A[t, ], vapply(cd$B, function(m) m[t, ], numeric(G)))
      Kt <- stats::cov(Ct)
      diag_d[t, ] <- diag(Kt)
      if (ncol(pairs)) off_d[t, ] <- Kt[cbind(pairs[1L, ], pairs[2L, ])]
    }
    for (i in seq_len(p1)) {
      rows[[length(rows) + 1L]] <-
        cbind(data.frame(term = paste0("var(", coef_lab[i], ")")),
              .post_summary(diag_d[, i], prob))
    }
    for (p in seq_len(ncol(pairs))) {
      nm <- paste0("cov(", coef_lab[pairs[1L, p]], ", ", coef_lab[pairs[2L, p]], ")")
      rows[[length(rows) + 1L]] <- cbind(data.frame(term = nm),
                                         .post_summary(off_d[, p], prob))
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

#' Reaction-norm analytics across the environmental gradient
#'
#' For a random-regression (reaction-norm) fit, summarises the fitted
#' reaction-norm coefficients across the environmental gradient. Writing each
#' genotype's curve as \eqn{v_i(x) = \phi(x)^\top c_i} with
#' \eqn{\phi(x) = (1, L_1(x), \dots, L_q(x))} (intercept column = 1, then the
#' orthonormal Legendre slopes) and coefficient vector
#' \eqn{c_i = (a_i, b_{i,1}, \dots, b_{i,q})}, it returns:
#'
#' * **Genetic covariance / correlation surface.** With
#'   \eqn{K = \mathrm{Var}(c_i)} the (co)variance matrix of the reaction-norm
#'   coefficients across genotypes (intercept--slope covariance included), the
#'   genetic covariance between two gradient points is
#'   \eqn{G(x_1, x_2) = \phi(x_1)^\top K\, \phi(x_2)}, i.e. the matrix
#'   \eqn{\Phi K \Phi^\top}; `cov2cor()` of it is the across-gradient genetic
#'   **correlation** surface (how a genotype ranking at one environment carries
#'   over to another).
#' * **Heritability over the gradient.** \eqn{h^2(x) = \sigma^2_g(x) /
#'   (\sigma^2_g(x) + \sigma^2_e)} with \eqn{\sigma^2_g(x) = \phi(x)^\top K \phi(x)}
#'   the diagonal of the covariance surface and \eqn{\sigma^2_e} the residual
#'   variance, computed per MCMC draw for a full posterior (mean + credible band).
#' * **Reliability over the gradient.** Per genotype,
#'   \eqn{r^2_i(x) = 1 - \mathrm{PEV}_i(x)/\sigma^2_g(x)}, where
#'   \eqn{\mathrm{PEV}_i(x)} is the posterior variance of \eqn{v_i(x)}.
#' * **Selection probability over the gradient.** For each gradient point, the
#'   posterior probability that each genotype ranks in the top `threshold`
#'   fraction *at that point* (each MCMC draw ranked independently, as in [pr()]).
#'
#' `K` is estimated from the posterior draws of the per-genotype coefficients
#' (the realised across-genotype coefficient (co)variance, averaged over draws),
#' so it carries an intercept--slope covariance even though `~ gen + gen:leg(x)`
#' fits the two terms with independent priors.
#'
#' @param fit A `breedRB_fit` (single-trait) containing a random regression.
#' @param term The interaction term identifier, e.g. `"gen:leg(x)"` or
#'   `"mrk(gen, M):leg(x, 1)"`. Its matching intercept term is detected
#'   automatically.
#' @param n_grid Integer; number of gradient points (default 20). The
#'   correlation surface is `n_grid x n_grid`; keep it modest for large panels.
#' @param threshold Selected fraction in `(0, 1)` for the selection
#'   probabilities; `0.20` = top 20%.
#' @param higher Logical (default `TRUE`); if `FALSE`, "top" means the smallest
#'   values (used for both the selection probability and genotype ordering).
#' @param leg_basis Logical (default `TRUE`). Gradient axis on the standardized
#'   Legendre domain \eqn{[-1, 1]}; if `FALSE`, back-transformed to the covariate
#'   scale stored at fitting.
#' @param top_n Optional integer; restrict the per-genotype reliability and
#'   probability heat maps to the `top_n` genotypes (by mean value across the
#'   gradient) so the figures stay legible. `NULL` (default) keeps all genotypes.
#' @param prob Central credible-interval mass for the heritability band
#'   (default 0.95).
#' @param plot Logical (default `TRUE`); draw and print the four \pkg{ggplot2}
#'   figures. They are always returned in the `plots` element.
#'
#' @return An object of class `breedRB_rrgradient` (a list) with: `grid` (gradient
#'   points on the chosen scale), `K` (posterior-mean coefficient (co)variance),
#'   `gcov` / `gcor` (`n_grid x n_grid` genetic covariance / correlation
#'   surfaces), `h2` (data frame `gradient`, `mean`, `median`, `lower`, `upper`),
#'   `reliability` and `prob_top` (`nGeno x n_grid` matrices), and `plots` (named
#'   list of ggplots: `cor`, `h2`, `reliability`, `prob`).
#' @examples
#' \donttest{
#' fit <- bbglr(yield ~ leg(x), random = ~ gen + gen:leg(x) + env:rep, data = dat)
#' g <- rr_gradient(fit, term = "gen:leg(x)", threshold = 0.20)
#' g$gcor          # across-gradient genetic correlation matrix
#' g$plots$prob    # who is likely to be top-20% where along the gradient
#' }
#' @seealso [reaction_norm()] for the per-genotype curves, [pr()] for the
#'   single-point selection probabilities, [heritability()].
#' @export
rr_gradient <- function(fit, term, n_grid = 20L, threshold = 0.20, higher = TRUE,
                        leg_basis = TRUE, top_n = NULL, prob = 0.95, plot = TRUE) {
  stopifnot(inherits(fit, "breedRB_fit"))
  if (isTRUE(fit$response$multitrait)) {
    stop("rr_gradient() currently supports single-trait fits.", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      threshold <= 0 || threshold >= 1) {
    stop("`threshold` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }

  cd  <- .rr_coef_draws(fit, term)
  ids <- cd$ids; q <- cd$q; G <- length(ids)
  A   <- cd$A; B <- cd$B
  nD  <- nrow(A)

  # Legendre design across the gradient: intercept column = 1, then L_1..L_q.
  gstd <- seq(-1, 1, length.out = n_grid)
  Bful <- legendre_basis(gstd, order = q, orthonormal = TRUE)[, -1, drop = FALSE]
  Phi  <- cbind(1, Bful)                                    # n_grid x (q + 1)
  gx   <- if (isTRUE(leg_basis)) gstd
          else cd$range[1] + (gstd + 1) / 2 * diff(cd$range)

  # --- K (coefficient (co)variance) and h2(x), per MCMC draw ------------------
  varE <- .varE_draws(fit)
  use_varE_draws <- length(varE) == nD
  if (!use_varE_draws) {
    varE_mean <- varcomp(fit)$mean[varcomp(fit)$term == "varE"]
    if (!length(varE_mean)) varE_mean <- NA_real_
  }

  Ksum   <- matrix(0, q + 1L, q + 1L)
  vgD    <- matrix(0, nD, n_grid)                           # per-draw genetic variance sigma^2_g(x)
  for (t in seq_len(nD)) {
    Ct <- cbind(A[t, ], vapply(B, function(m) m[t, ], numeric(G)))   # G x (q + 1)
    Kt <- stats::cov(Ct)
    Ksum <- Ksum + Kt
    vgD[t, ] <- rowSums((Phi %*% Kt) * Phi)                 # diag(Phi Kt Phi')
  }
  K    <- Ksum / nD
  gcov <- Phi %*% K %*% t(Phi)                              # posterior-mean covariance surface
  gcor <- stats::cov2cor(gcov)
  dimnames(gcov) <- dimnames(gcor) <- list(round(gx, 3), round(gx, 3))

  h2D <- if (use_varE_draws) vgD / (vgD + varE)             # [nD x n_grid]
         else                vgD / (vgD + varE_mean)
  h2  <- data.frame(gradient = gx,
                    mean   = colMeans(h2D),
                    median = apply(h2D, 2L, stats::median),
                    lower  = apply(h2D, 2L, stats::quantile, probs = (1 - prob) / 2),
                    upper  = apply(h2D, 2L, stats::quantile, probs = 1 - (1 - prob) / 2),
                    row.names = NULL)

  # --- reliability(x) and top-fraction probability(x), per gradient point -----
  vg   <- diag(gcov)                                        # population genetic var at each x
  sgn  <- if (isTRUE(higher)) -1 else 1
  k    <- max(1L, round(threshold * G))
  Rel  <- matrix(NA_real_, G, n_grid, dimnames = list(ids, round(gx, 3)))
  Prob <- matrix(NA_real_, G, n_grid, dimnames = list(ids, round(gx, 3)))
  meanGV <- matrix(NA_real_, G, n_grid, dimnames = list(ids, NULL))
  for (g in seq_len(n_grid)) {
    GV <- A                                                 # [nD x G]; phi[,1] = 1
    for (j in seq_len(q)) GV <- GV + Phi[g, j + 1L] * B[[j]]
    meanGV[, g] <- colMeans(GV)
    pev         <- apply(GV, 2L, stats::var)
    Rel[, g]    <- pmin(pmax(1 - pev / vg[g], 0), 1)
    flag        <- apply(GV, 1L, function(r) rank(sgn * r, ties.method = "min") <= k)
    Prob[, g]   <- rowMeans(flag)                           # flag is [G x nD]
  }

  # Genotype ordering for the heat maps: by mean value across the gradient.
  ord <- order(rowMeans(meanGV), decreasing = isTRUE(higher))
  sel <- if (is.null(top_n)) ord else ord[seq_len(min(top_n, G))]

  out <- structure(
    list(grid = gx, gstd = gstd, K = K, gcov = gcov, gcor = gcor, h2 = h2,
         reliability = Rel[ord, , drop = FALSE],
         prob_top    = Prob[ord, , drop = FALSE],
         ids = ids, q = q, legvar = cd$legvar, genomic = cd$genomic,
         threshold = threshold, higher = higher, leg_basis = leg_basis,
         term = fit$meta[[.resolve_term(fit, term)]]$label),
    class = "breedRB_rrgradient")

  out$plots <- .rr_gradient_plots(out, sel_ids = ids[sel])
  if (isTRUE(plot)) for (p in out$plots) print(p)
  out
}

#' Build the four ggplot figures for an rr_gradient() result
#' @keywords internal
.rr_gradient_plots <- function(x, sel_ids) {
  xlab <- if (isTRUE(x$leg_basis))
    paste0("Legendre gradient  [-1, 1]  (", x$legvar, ")") else x$legvar

  # 1. Genetic correlation surface (heat map).
  cor_df <- expand.grid(x1 = x$grid, x2 = x$grid)
  cor_df$gcor <- as.numeric(x$gcor)
  p_cor <- ggplot2::ggplot(cor_df, ggplot2::aes(x1, x2, fill = gcor)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradient2(midpoint = 0, limits = c(-1, 1),
                                  low = "#b2182b", mid = "white", high = "#2166ac",
                                  name = "genetic\ncorrelation") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = xlab, y = xlab,
                  title = paste0("Across-gradient genetic correlation: ", x$term)) +
    ggplot2::theme_minimal()

  # 2. Heritability over the gradient (posterior mean + band).
  p_h2 <- ggplot2::ggplot(x$h2, ggplot2::aes(x = gradient, y = mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = "steelblue", alpha = 0.25) +
    ggplot2::geom_line(colour = "steelblue", linewidth = 0.8) +
    ggplot2::ylim(0, NA) +
    ggplot2::labs(x = xlab, y = expression(h^2),
                  title = paste0("Heritability across the gradient: ", x$term)) +
    ggplot2::theme_minimal()

  # 3. Reliability heat map (genotype x gradient).
  rel <- x$reliability[rownames(x$reliability) %in% sel_ids, , drop = FALSE]
  rel_df <- expand.grid(id = factor(rownames(rel), levels = rev(rownames(rel))),
                        gradient = x$grid)
  rel_df$reliability <- as.numeric(rel)
  p_rel <- ggplot2::ggplot(rel_df, ggplot2::aes(gradient, id, fill = reliability)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), name = "reliability") +
    ggplot2::labs(x = xlab, y = "genotype",
                  title = paste0("Reliability across the gradient: ", x$term)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = if (nrow(rel) > 40)
      ggplot2::element_blank() else ggplot2::element_text(size = 6))

  # 4. Top-fraction selection probability heat map.
  pt <- x$prob_top[rownames(x$prob_top) %in% sel_ids, , drop = FALSE]
  pt_df <- expand.grid(id = factor(rownames(pt), levels = rev(rownames(pt))),
                       gradient = x$grid)
  pt_df$prob <- as.numeric(pt)
  p_prob <- ggplot2::ggplot(pt_df, ggplot2::aes(gradient, id, fill = prob)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "magma",
                                  name = sprintf("P(top %g%%)", 100 * x$threshold)) +
    ggplot2::labs(x = xlab, y = "genotype",
                  title = sprintf("P(top %g%%) across the gradient: %s",
                                  100 * x$threshold, x$term)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = if (nrow(pt) > 40)
      ggplot2::element_blank() else ggplot2::element_text(size = 6))

  list(cor = p_cor, h2 = p_h2, reliability = p_rel, prob = p_prob)
}

#' @export
print.breedRB_rrgradient <- function(x, ...) {
  cat("<breedRB_rrgradient>  term:", x$term,
      if (x$genomic) "(genomic)" else "(factor)", "\n")
  cat(sprintf("  Gradient points: %d on %s scale\n", length(x$grid),
              if (x$leg_basis) "Legendre [-1,1]" else "covariate"))
  cat(sprintf("  Genotypes: %d   leg() order q = %d\n", length(x$ids), x$q))
  cat("  Coefficient (co)variance K (intercept, slope1, ...):\n")
  print(round(x$K, 4))
  hr <- range(x$gcor[upper.tri(x$gcor)])
  cat(sprintf("  Across-gradient genetic correlation: %.3f to %.3f\n", hr[1], hr[2]))
  cat(sprintf("  Heritability over gradient: %.3f to %.3f (posterior mean)\n",
              min(x$h2$mean), max(x$h2$mean)))
  cat(sprintf("  Selection probabilities: P(top %g%%), %s is better\n",
              100 * x$threshold, if (x$higher) "higher" else "lower"))
  cat("  $plots: cor, h2, reliability, prob\n")
  invisible(x)
}
