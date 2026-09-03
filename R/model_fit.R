# ---------------------------------------------------------------------------
# Goodness-of-fit helpers: fitted values / residuals from BGLR's posterior-mean
# fit, and a one-call model_fit() summary (R^2, RMSE, DIC, per-term reliability).
# ---------------------------------------------------------------------------

#' Pooled posterior-mean fitted values, aligned to the modelled data
#' @keywords internal
.pooled_yhat <- function(fit) {
  yh <- lapply(fit$chains, `[[`, "yHat")
  if (isTRUE(fit$response$multitrait)) {
    Reduce(`+`, yh) / length(yh)                       # [n x nTrait]
  } else {
    rowMeans(do.call(cbind, yh))                       # length n
  }
}

#' Observed response(s) used in the fit
#' @keywords internal
.observed <- function(fit) {
  tr <- fit$response$traits
  if (isTRUE(fit$response$multitrait)) as.matrix(fit$data[, tr, drop = FALSE])
  else fit$data[[tr]]
}

#' Fitted values from a `breedRB_fit`
#'
#' Posterior-mean fitted values from the BGLR engine (`yHat`), pooled across
#' chains. These include **every** model term — the fixed part, all random terms
#' and, for a random regression, the fitted reaction norm at each observation's
#' covariate value. Rows are aligned to the data actually modelled
#' (`fit$data`), so missing responses (which BGLR imputes) still receive a fitted
#' value.
#'
#' @param object A `breedRB_fit`.
#' @param ... Unused; for S3 compatibility.
#' @return A numeric vector of fitted values (single-trait), or an
#'   `n x nTrait` matrix (multi-trait). The per-observation posterior SD
#'   (`SD.yHat`, averaged across chains) is attached as attribute `"sd"`.
#' @seealso [residuals.breedRB_fit()], [model_fit()].
#' @export
fitted.breedRB_fit <- function(object, ...) {
  yhat <- .pooled_yhat(object)
  sdl  <- lapply(object$chains, function(c) c$SD.yHat)
  if (!any(vapply(sdl, is.null, logical(1)))) {
    attr(yhat, "sd") <- if (isTRUE(object$response$multitrait))
      Reduce(`+`, sdl) / length(sdl) else rowMeans(do.call(cbind, sdl))
  }
  yhat
}

#' Residuals from a `breedRB_fit`
#'
#' Observed response minus the posterior-mean fitted value ([fitted.breedRB_fit()]).
#' Observations with a missing response (imputed by BGLR) return `NA`.
#'
#' @param object A `breedRB_fit`.
#' @param ... Unused; for S3 compatibility.
#' @return A numeric vector (single-trait) or `n x nTrait` matrix (multi-trait)
#'   of residuals.
#' @seealso [fitted.breedRB_fit()], [model_fit()].
#' @export
residuals.breedRB_fit <- function(object, ...) {
  obs  <- .observed(object)
  yhat <- .pooled_yhat(object)
  res  <- obs - yhat
  res[!is.finite(obs)] <- NA_real_          # response was missing -> residual undefined
  res
}

#' Goodness-of-fit summary for a fitted model
#'
#' Bundles the diagnostics most useful for judging how well a model — in
#' particular a random-regression / reaction-norm model — fits the data:
#' the observed-vs-fitted \eqn{R^2} and RMSE (on the observations with a
#' non-missing response), the deviance information criterion (`DIC`) and
#' effective number of parameters (`pD`) for model comparison, the posterior-mean
#' residual variance, and, for every random term, the mean and median BLUP
#' reliability (from [solution()]). A random-regression interaction is reported
#' **per Legendre degree** (`...:deg1`, `...:deg2`, ...) when its order exceeds 1,
#' so the drop in reliability from the linear slope to the higher-order
#' curvature coefficients is visible. Lower `DIC` indicates a better
#' complexity-penalised fit; comparing the `DIC` of the full random regression
#' against a slope-free model shows whether the reaction-norm (`gen:leg(x)`) term
#' is warranted.
#'
#' @param fit A `breedRB_fit` (single-trait).
#' @return An object of class `breedRB_modelfit` (a list with `n`, `r2`, `rmse`,
#'   `dic`, `pD`, `varE`, and a `reliability` data frame of per-term
#'   `mean`/`median`/`n_effects` — one row per Legendre degree for a
#'   random-regression interaction of order > 1), with a `print` method. The
#'   fitted values and
#'   residuals are attached as attributes `"fitted"` / `"residuals"`.
#' @examples
#' \donttest{
#' fit_rr <- bbglr(yield ~ leg(x), random = ~ gen + gen:leg(x) + env:rep,
#'                 data = dat)
#' model_fit(fit_rr)
#' # justify the random regression: compare DIC to a slope-free model
#' fit0 <- bbglr(yield ~ leg(x), random = ~ gen + env:rep, data = dat)
#' model_fit(fit_rr)$dic < model_fit(fit0)$dic
#' }
#' @seealso [fitted.breedRB_fit()], [residuals.breedRB_fit()], [mcmc_diag()],
#'   [heritability()], [varcomp()].
#' @export
model_fit <- function(fit) {
  stopifnot(inherits(fit, "breedRB_fit"))
  if (isTRUE(fit$response$multitrait)) {
    stop("model_fit() currently supports single-trait fits; inspect ",
         "fitted()/residuals() and DIC per chain for multi-trait models.",
         call. = FALSE)
  }

  obs  <- .observed(fit)
  yhat <- .pooled_yhat(fit)
  ok   <- is.finite(obs)
  res  <- obs - yhat

  r2   <- if (sum(ok) > 1L) stats::cor(obs[ok], yhat[ok])^2 else NA_real_
  rmse <- sqrt(mean(res[ok]^2))
  dic  <- mean(vapply(fit$chains, function(c) c$fit$DIC, numeric(1)))
  pD   <- mean(vapply(fit$chains, function(c) c$fit$pD,  numeric(1)))

  vc   <- varcomp(fit)
  varE <- vc$mean[vc$term == "varE"]
  if (!length(varE)) varE <- NA_real_

  # Per random-term BLUP reliability (skip terms with no reliability column). A
  # random-regression interaction is broken out per Legendre degree (deg1 =
  # linear slope, deg2 = quadratic, ...), since higher degrees are estimated far
  # less reliably than the intercept/linear coefficients.
  rk  <- .random_keys(fit)
  rel <- lapply(rk, function(k) {
    s <- tryCatch(solution(fit, term = k, type = "random"), error = function(e) NULL)
    if (is.null(s) || is.null(s$reliability)) return(NULL)
    lab <- fit$meta[[k]]$label
    deg  <- ifelse(grepl(":deg[0-9]+$", s$effect),
                   as.integer(sub(".*:deg([0-9]+)$", "\\1", s$effect)), NA_integer_)
    degs <- sort(unique(deg[!is.na(deg)]))
    if (length(degs) > 1L) {                      # RR interaction (q > 1): one row per degree
      do.call(rbind, lapply(degs, function(j) {
        ix <- which(deg == j)
        data.frame(term = paste0(lab, ":deg", j),
                   mean_reliability   = mean(s$reliability[ix]),
                   median_reliability = stats::median(s$reliability[ix]),
                   n_effects = length(ix), row.names = NULL)
      }))
    } else {
      data.frame(term = lab,
                 mean_reliability   = mean(s$reliability),
                 median_reliability = stats::median(s$reliability),
                 n_effects = nrow(s), row.names = NULL)
    }
  })
  rel <- do.call(rbind, Filter(Negate(is.null), rel))

  structure(
    list(n = sum(ok), r2 = r2, rmse = rmse, dic = dic, pD = pD,
         varE = as.numeric(varE), reliability = rel),
    fitted    = yhat,
    residuals = { res[!ok] <- NA_real_; res },
    class     = "breedRB_modelfit"
  )
}

#' @export
print.breedRB_modelfit <- function(x, ...) {
  cat("<breedRB_modelfit>\n")
  cat(sprintf("  Observations (non-missing): %d\n", x$n))
  cat(sprintf("  Observed-vs-fitted R2: %.4f   RMSE: %.4g\n", x$r2, x$rmse))
  cat(sprintf("  DIC: %.1f   pD: %.1f   residual var: %.4g\n", x$dic, x$pD, x$varE))
  if (!is.null(x$reliability) && nrow(x$reliability)) {
    cat("  Random-term BLUP reliability:\n")
    r <- x$reliability
    for (i in seq_len(nrow(r))) {
      cat(sprintf("    %-24s mean %.3f  median %.3f  (%d effects)\n",
                  r$term[i], r$mean_reliability[i], r$median_reliability[i],
                  r$n_effects[i]))
    }
  }
  invisible(x)
}
