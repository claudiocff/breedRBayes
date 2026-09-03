# ---------------------------------------------------------------------------
# MCMC convergence diagnostics, built on `coda`.
# ---------------------------------------------------------------------------

#' Markov-chain convergence diagnostics for a fitted model
#'
#' Summarises convergence of the variance-component chains: effective sample
#' size, Geweke's z, and (when the fit has >1 chain) the Gelman-Rubin
#' potential-scale-reduction factor (R-hat).
#'
#' @param fit A `breedRB_fit`.
#' @param what Passed to [as_mcmc()] (default `"varcomp"`; add `"mu"` to include
#'   the intercept).
#' @param plot Logical. If `TRUE`, also render `ggplot2` trace plots of the
#'   monitored chains (via [plot_trace()]) and attach the plot object to the
#'   returned data frame as the `"plot"` attribute. Default `FALSE`.
#' @return A data frame with one row per monitored parameter: `param`, `n_eff`,
#'   `geweke_z`, and `Rhat` (`NA` for single-chain fits). When `plot = TRUE`,
#'   the corresponding `ggplot` is attached as `attr(., "plot")`.
#' @examples
#' \donttest{
#' mcmc_diag(fit)
#' d <- mcmc_diag(fit, plot = TRUE)   # prints trace plots
#' attr(d, "plot")                    # the ggplot object
#' }
#' @export
mcmc_diag <- function(fit, what = "varcomp", plot = FALSE) {
  stopifnot(inherits(fit, "breedRB_fit"))
  mc <- as_mcmc(fit, what = what)
  params <- coda::varnames(mc)

  n_eff <- coda::effectiveSize(mc)

  gew <- tryCatch(coda::geweke.diag(mc[[1]])$z, error = function(e) rep(NA_real_, length(params)))

  rhat <- rep(NA_real_, length(params))
  if (length(mc) > 1L) {
    gd <- tryCatch(coda::gelman.diag(mc, multivariate = FALSE, autoburnin = FALSE),
                   error = function(e) NULL)
    if (!is.null(gd)) rhat <- gd$psrf[, 1]
  }

  out <- data.frame(param = params,
                    n_eff = as.numeric(n_eff[params]),
                    geweke_z = as.numeric(gew[params]),
                    Rhat = as.numeric(rhat),
                    row.names = NULL)

  if (isTRUE(plot)) {
    p <- plot_trace(fit, what = what)
    print(p)
    attr(out, "plot") <- p
  }
  out
}

#' Trace plots of the variance-component chains
#'
#' @param fit A `breedRB_fit`.
#' @param what Passed to [as_mcmc()].
#' @return A `ggplot` object (iteration on x, value on y, one panel per parameter,
#'   coloured by chain).
#' @export
plot_trace <- function(fit, what = "varcomp") {
  mc <- as_mcmc(fit, what = what)
  df <- do.call(rbind, lapply(seq_along(mc), function(i) {
    m <- as.matrix(mc[[i]])
    data.frame(iter = rep(seq_len(nrow(m)), times = ncol(m)),
               value = as.vector(m),
               param = rep(colnames(m), each = nrow(m)),
               chain = factor(i))
  }))
  ggplot2::ggplot(df, ggplot2::aes(.data$iter, .data$value, colour = .data$chain)) +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::facet_wrap(~ param, scales = "free_y") +
    ggplot2::labs(x = "Iteration", y = "Value", colour = "Chain") +
    ggplot2::theme_bw()
}

#' Posterior density plot of variance components or heritability
#'
#' @param x A `breedRB_fit` (variance components) or a `breedRB_h2` object.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
plot_posterior <- function(x, ...) UseMethod("plot_posterior")

#' @rdname plot_posterior
#' @export
plot_posterior.breedRB_fit <- function(x, ...) {
  d <- attr(varcomp(x, draws = TRUE), "draws")
  df <- data.frame(value = as.vector(d),
                   param = rep(colnames(d), each = nrow(d)))
  ggplot2::ggplot(df, ggplot2::aes(.data$value)) +
    ggplot2::geom_density(fill = "#1F77B4", alpha = 0.5) +
    ggplot2::facet_wrap(~ param, scales = "free") +
    ggplot2::labs(x = "Posterior value", y = "Density") +
    ggplot2::theme_bw()
}

#' @rdname plot_posterior
#' @export
plot_posterior.breedRB_h2 <- function(x, ...) {
  ggplot2::ggplot(data.frame(h2 = x$draws), ggplot2::aes(.data$h2)) +
    ggplot2::geom_density(fill = "#2CA02C", alpha = 0.5) +
    ggplot2::geom_vline(xintercept = x$summary$median, linetype = 2) +
    ggplot2::labs(x = expression(h^2), y = "Density") +
    ggplot2::theme_bw()
}
