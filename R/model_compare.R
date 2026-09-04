# ---------------------------------------------------------------------------
# Model comparison for (nested) random-regression models: WAIC / PSIS-LOO, an
# approximate deviance "pseudo-LRT", and DIC, via an anova() method.
#
# A classical likelihood-ratio test is not strictly valid here: (1) the fit is
# Bayesian (BGLR Gibbs) with no maximised likelihood, and (2) adding a Legendre
# degree tests variance components on the boundary of the parameter space, so the
# null distribution is a chi-bar-squared mixture rather than a plain chi-square.
# The defensible analogues are the predictive criteria (WAIC / LOO, with the SE of
# their difference) and DIC. The deviance pseudo-LRT is reported only as a rough
# heuristic and is flagged as such.
# ---------------------------------------------------------------------------

#' Pointwise log-likelihood matrix (draws x observations) for a Gaussian fit
#'
#' Reconstructs the per-draw linear predictor and evaluates the Gaussian
#' pointwise log-likelihood `log p(y_i | theta^(s))` for every retained MCMC
#' draw `s`, over the observations with a non-missing response. This is the input
#' WAIC and PSIS-LOO need.
#'
#' Deterministic term designs (fixed effects, plain random factors, Legendre
#' bases, covariates, and their interactions) are rebuilt directly and multiplied
#' by BGLR's saved coefficient draws. `vm()` / `mrk()` genomic terms are **not**
#' rebuilt — their principal-component rotation is data-dependent and, under the
#' default `RSpectra` truncation, not sign-reproducible — so their per-observation
#' contribution is taken from [solution()], which back-maps through the rotation
#' stored on the fit. Draws are pooled across chains in `fit$paths` order,
#' matching how [solution()] pools, so all pieces align draw-by-draw.
#'
#' Heterogeneous residuals (`dsum(~units | g)`) are supported: each observation
#' is scored against its own residual-variance group's per-draw variance.
#'
#' @param fit A single-trait `breedRB_fit` (homogeneous or heterogeneous
#'   residual).
#' @return A numeric matrix `[nDraws x nObs]` of pointwise log-likelihoods, with
#'   attributes `"chain_id"` (chain of each draw row), `"S_per_chain"` and
#'   `"obs_index"` (which rows of `fit$data` the columns correspond to).
#' @keywords internal
.pointwise_loglik <- function(fit) {
  if (isTRUE(fit$response$multitrait)) {
    stop("WAIC/LOO currently support single-trait fits.", call. = FALSE)
  }
  hetero <- !is.null(fit$group)                    # dsum(~units | g): per-group varE

  y  <- fit$data[[fit$response$traits]]
  ok <- is.finite(y)
  paths  <- fit$paths
  C      <- length(paths)
  nburn  <- floor(fit$control$burnIn / fit$control$thin)
  ev     <- fit$control$exp_var_rank
  exp_var <- if (is.null(ev) || is.na(ev)) NA_real_ else ev

  # Deterministic designs (genomic terms rebuilt here are ignored below).
  built <- suppressMessages(build_eta(fit$parsed, fit$data, fit$relmat, exp_var = exp_var))
  keys  <- names(fit$meta)

  is_genomic <- vapply(keys, function(k) {
    any(vapply(fit$meta[[k]]$components,
               function(cm) isTRUE(cm$kind %in% c("vm", "mrk")), logical(1)))
  }, logical(1))

  # --- per-draw intercept and residual variance, pooled in fit$paths order ----
  # varE.dat holds one column per residual group (a single column when the
  # residual is homogeneous); keep it as a matrix so hetero fits are handled too.
  S_c <- integer(C); mu <- numeric(0); vE_list <- vector("list", C)
  for (i in seq_len(C)) {
    mu_all <- scan(paste0(paths[i], "mu.dat"), quiet = TRUE)
    vE_all <- as.matrix(utils::read.table(paste0(paths[i], "varE.dat")))  # [rows x G]
    s <- length(mu_all) - nburn
    if (s <= 0) stop("No post-burn-in draws found for chain ", i, ".", call. = FALSE)
    S_c[i] <- s
    mu <- c(mu, utils::tail(mu_all, s))
    vE_list[[i]] <- vE_all[(nrow(vE_all) - s + 1L):nrow(vE_all), , drop = FALSE]
  }
  if (length(unique(S_c)) != 1L) {
    stop("Chains have unequal numbers of retained draws; cannot pool.", call. = FALSE)
  }
  varE <- do.call(rbind, vE_list)                  # [S x G]  (G = 1 if homogeneous)
  S <- sum(S_c); n <- nrow(fit$data)

  # --- sum every term's per-observation contribution [n x S] ------------------
  contrib <- matrix(0, n, S)
  for (k in keys) {
    if (is_genomic[[k]]) {
      contrib <- contrib + .genomic_obs_contrib(fit, k)
      next
    }
    Xk <- built$ETA[[k]]$X
    isFixed <- identical(fit$meta[[k]]$model, "FIXED")
    Blist <- lapply(paths, function(prefix) {
      if (isFixed) {
        B <- as.matrix(utils::read.table(paste0(prefix, "ETA_", k, "_b.dat"),
                                         header = TRUE))
        if (nrow(B) > nburn) B <- B[(nburn + 1L):nrow(B), , drop = FALSE]
        B
      } else {
        BGLR::readBinMat(paste0(prefix, "ETA_", k, "_b.bin"))       # [S_c x p]
      }
    })
    Bpool <- do.call(rbind, Blist)                                   # [S x p]
    if (nrow(Bpool) != S) {
      stop("Draw count mismatch reconstructing term '", k, "'.", call. = FALSE)
    }
    contrib <- contrib + tcrossprod(Xk, Bpool)                       # X %*% t(B) = [n x S]
  }

  yhat <- sweep(contrib, 2L, mu, "+")                                # add mu per draw

  # --- Gaussian pointwise log-likelihood on the observed rows -----------------
  yo <- y[ok]; YH <- yhat[ok, , drop = FALSE]; n_ok <- sum(ok)
  resid <- YH - yo                                                   # (i,s) - yo_i (column recycle)
  # Per-observation residual-variance trace [n_ok x S]: each observed row uses
  # its residual group's varE column (all rows share column 1 when homogeneous).
  if (hetero) {
    grp_ok   <- fit$group[which(ok)]                                 # group index per obs row
    varE_obs <- t(varE[, grp_ok, drop = FALSE])                      # [n_ok x S]
  } else {
    varE_obs <- matrix(varE[, 1L], n_ok, S, byrow = TRUE)            # [n_ok x S]
  }
  ll  <- -0.5 * log(2 * pi) - 0.5 * log(varE_obs) - resid^2 / (2 * varE_obs)
  out <- t(ll)                                                       # [S x n_ok]
  attr(out, "chain_id")    <- rep(seq_len(C), times = S_c)
  attr(out, "S_per_chain") <- S_c[1]
  attr(out, "obs_index")   <- which(ok)
  out
}

#' Per-observation contribution of a genomic (`vm()`/`mrk()`) term, all draws
#'
#' Uses [solution()] (which back-maps through the stored PC rotation, so it is
#' sign-safe) to obtain per-genotype BLUP draws (a main effect) or per-genotype,
#' per-degree reaction-norm coefficient draws (a genomic random regression), then
#' maps them to observations — directly for a main effect, or via the same
#' orthonormal Legendre basis used at fitting for a random regression.
#' @keywords internal
.genomic_obs_contrib <- function(fit, key) {
  meta  <- fit$meta[[key]]
  comps <- meta$components
  gi    <- which(vapply(comps, function(cm) isTRUE(cm$kind %in% c("vm", "mrk")), logical(1)))
  li    <- which(vapply(comps, function(cm) identical(cm$kind, "leg"), logical(1)))
  gen_var <- comps[[gi[1]]]$var
  gof     <- as.character(fit$data[[gen_var]])

  s <- solution(fit, term = key, type = "random")
  D <- attr(s, "draws")                                             # [S x nEff]
  # NB: solution() sorts its data-frame rows, but the `draws` columns stay in the
  # term's natural level order — always map through colnames(D), never s$effect.
  eff <- colnames(D)

  if (!length(li)) {                                                # genomic main effect
    idx <- match(gof, eff)
    Cn  <- t(D[, idx, drop = FALSE])                                # [n x S]
    Cn[is.na(idx), ] <- 0
    return(Cn)
  }

  legm <- comps[[li[1]]]
  q    <- legm$order
  xs   <- .scale_unit(as.numeric(fit$data[[legm$var]]), rng = legm$range)
  L    <- legendre_basis(xs, order = q, orthonormal = TRUE)[, -1, drop = FALSE]  # [n x q]
  deg  <- ifelse(grepl(":deg[0-9]+$", eff),
                 as.integer(sub(".*:deg([0-9]+)$", "\\1", eff)),
                 if (q == 1L) 1L else NA_integer_)
  id   <- sub(":deg[0-9]+$", "", eff)

  n <- nrow(fit$data); S <- nrow(D); out <- matrix(0, n, S)
  for (j in seq_len(q)) {
    cj  <- which(deg == j)
    col <- cj[match(gof, id[cj])]                                   # D column per obs at degree j
    Dj  <- t(D[, col, drop = FALSE]); Dj[is.na(col), ] <- 0         # [n x S]
    out <- out + Dj * L[, j]                                        # row i scaled by L[i, j]
  }
  out
}

#' WAIC from a pointwise log-likelihood matrix (fallback when `loo` is absent)
#' @keywords internal
.waic_manual <- function(ll) {
  S <- nrow(ll)
  lpd    <- apply(ll, 2L, function(col) { m <- max(col); m + log(mean(exp(col - m))) })
  p_waic <- apply(ll, 2L, stats::var)
  elpd_i <- lpd - p_waic
  list(elpd_waic = sum(elpd_i), p_waic = sum(p_waic),
       waic = -2 * sum(elpd_i), se_waic = sqrt(length(elpd_i) * stats::var(2 * elpd_i)))
}

#' WAIC and PSIS-LOO for one fit, preferring the `loo` package
#' @keywords internal
.criteria_one <- function(fit, want_loo = TRUE) {
  ll  <- .pointwise_loglik(fit)
  S_c <- attr(ll, "S_per_chain"); C <- length(unique(attr(ll, "chain_id")))
  n   <- ncol(ll)
  res <- list(n = n, waic = NA_real_, p_waic = NA_real_, se_waic = NA_real_,
              looic = NA_real_, elpd_loo = NA_real_, p_loo = NA_real_,
              se_looic = NA_real_, max_pareto_k = NA_real_, loo_obj = NULL)

  if (requireNamespace("loo", quietly = TRUE)) {
    arr <- array(0, dim = c(S_c, C, n))
    for (ci in seq_len(C)) arr[, ci, ] <- ll[attr(ll, "chain_id") == ci, ]
    w <- loo::waic(arr)
    res$waic    <- w$estimates["waic", "Estimate"]
    res$p_waic  <- w$estimates["p_waic", "Estimate"]
    res$se_waic <- w$estimates["waic", "SE"]
    if (want_loo) {
      lo <- loo::loo(arr)
      res$looic       <- lo$estimates["looic", "Estimate"]
      res$se_looic    <- lo$estimates["looic", "SE"]
      res$elpd_loo    <- lo$estimates["elpd_loo", "Estimate"]
      res$p_loo       <- lo$estimates["p_loo", "Estimate"]
      res$max_pareto_k <- max(loo::pareto_k_values(lo), na.rm = TRUE)
      res$loo_obj     <- lo
    }
  } else {                                                          # manual WAIC only
    w <- .waic_manual(ll)
    res$waic <- w$waic; res$p_waic <- w$p_waic; res$se_waic <- w$se_waic
  }
  res
}

#' Compare fitted models: WAIC / PSIS-LOO, DIC, and an approximate deviance test
#'
#' Compares two or more models fitted with [bbglr()] — typically a ladder of
#' random-regression reaction norms of increasing Legendre degree — using the
#' criteria appropriate for a Bayesian fit. It reports, per model, the deviance
#' information criterion (`DIC`) and effective number of parameters (`pD`), the
#' widely-applicable information criterion (`WAIC`), and, when the \pkg{loo}
#' package is installed, PSIS-LOO (`LOOIC`) with the standard error of each
#' pairwise difference (via [loo::loo_compare()]) — the most defensible
#' "is-the-extra-degree-worth-it" comparison here. Lower `DIC`/`WAIC`/`LOOIC` is
#' better.
#'
#' It also prints an **approximate deviance ratio** between consecutive models (in
#' the order supplied): the drop in posterior-mean deviance `Dbar = DIC - pD`
#' against a chi-square with degrees of freedom equal to the increase in `pD`.
#' This is a rough heuristic only. A true likelihood-ratio test does not apply to
#' these Bayesian fits, and even under REML the added random-regression
#' (co)variances are tested on the boundary of the parameter space, so the naive
#' chi-square reference distribution is wrong (a chi-bar-squared mixture is
#' needed). Prefer the WAIC/LOO comparison and its difference SE; treat the
#' deviance ratio and its "p-value" as indicative, not a formal test.
#'
#' Pass the models in order of increasing complexity (degree 0, 1, 2, ...) so the
#' consecutive deviance ratios line up with adding one Legendre degree at a time.
#'
#' @param object A `breedRB_fit`.
#' @param ... Further `breedRB_fit` objects to compare against `object`.
#' @param loo Logical; compute PSIS-LOO in addition to WAIC (requires the
#'   \pkg{loo} package). Default `TRUE`.
#' @return An object of class `breedRB_anova` (a list with `table`, the deviance
#'   `lrt` data frame, and — when available — the `loo_compare` matrix), with a
#'   `print` method.
#' @examples
#' \donttest{
#' m0 <- bbglr(y ~ leg(x), random = ~ mrk(gen, M), data = dat, relmat = list(M = M))
#' m1 <- bbglr(y ~ leg(x), random = ~ mrk(gen, M) + mrk(gen, M):leg(x, 1),
#'             data = dat, relmat = list(M = M))
#' m2 <- bbglr(y ~ leg(x), random = ~ mrk(gen, M) + mrk(gen, M):leg(x, 2),
#'             data = dat, relmat = list(M = M))
#' anova(m0, m1, m2)
#' }
#' @seealso [model_fit()], [heritability()], [varcomp()].
#' @export
anova.breedRB_fit <- function(object, ..., loo = TRUE) {
  fits <- c(list(object), list(...))
  fits <- Filter(function(f) inherits(f, "breedRB_fit"), fits)
  if (length(fits) < 2L) {
    stop("anova(): provide at least two fitted models to compare.", call. = FALSE)
  }
  k <- length(fits)
  model <- paste0("model", seq_len(k))
  terms_lab <- vapply(fits, function(f) {
    r   <- f$call$random
    lab <- if (is.null(r)) "~1" else paste(deparse(r), collapse = "")
    # tag the residual structure so hetero-vs-homo comparisons are distinguishable
    if (isTRUE(f$parsed$residual$hetero)) {
      lab <- paste0(lab, "  [R: het by ", f$parsed$residual$group, "]")
    }
    lab
  }, character(1))

  dic <- vapply(fits, function(f) mean(vapply(f$chains, function(c) c$fit$DIC, numeric(1))), numeric(1))
  pD  <- vapply(fits, function(f) mean(vapply(f$chains, function(c) c$fit$pD,  numeric(1))), numeric(1))
  dbar <- dic - pD                                        # BGLR: DIC = Dbar + pD

  # WAIC / LOO (degrade gracefully if a fit can't be reconstructed)
  crit <- lapply(fits, function(f) tryCatch(.criteria_one(f, want_loo = loo),
                                            error = function(e) {
    warning("anova(): WAIC/LOO unavailable for a model (", conditionMessage(e), ").",
            call. = FALSE)
    NULL
  }))
  getc <- function(nm) vapply(crit, function(x) if (is.null(x)) NA_real_ else x[[nm]], numeric(1))
  n    <- vapply(fits, function(f) sum(is.finite(.observed(f))), integer(1))
  waic <- getc("waic"); looic <- getc("looic")
  # gap to the best model; all-NA columns (e.g. WAIC/LOO not computed) stay NA
  # instead of tripping min(na.rm=TRUE) on an empty set.
  dbest <- function(z) if (all(is.na(z))) z else z - min(z, na.rm = TRUE)

  tab <- data.frame(
    model = model, terms = terms_lab, n = n,
    pD = round(pD, 2), DIC = round(dic, 1), dDIC = round(dbest(dic), 1),
    WAIC = round(waic, 1), dWAIC = round(dbest(waic), 1),
    LOOIC = round(looic, 1), dLOOIC = round(dbest(looic), 1),
    max_pareto_k = round(getc("max_pareto_k"), 2),
    row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE)

  # Approximate deviance ratio between consecutive models (heuristic).
  lrt <- NULL
  if (k >= 2L) {
    comp <- character(k - 1L); dev <- df <- pval <- numeric(k - 1L)
    for (i in seq_len(k - 1L)) {
      comp[i] <- paste0(model[i], " -> ", model[i + 1L])
      dev[i]  <- dbar[i] - dbar[i + 1L]                   # deviance reduction
      df[i]   <- pD[i + 1L] - pD[i]                       # extra effective params
      pval[i] <- if (is.finite(df[i]) && df[i] > 0 && dev[i] > 0)
                   stats::pchisq(dev[i], df[i], lower.tail = FALSE) else NA_real_
    }
    lrt <- data.frame(comparison = comp,
                      deviance_reduction = round(dev, 2),
                      df = round(df, 2),
                      p_value_approx = signif(pval, 3),
                      row.names = NULL, stringsAsFactors = FALSE)
  }

  # Exact LOO difference + SE (the recommended comparison) when loo is available.
  loo_cmp <- NULL
  loo_objs <- lapply(crit, function(x) if (is.null(x)) NULL else x$loo_obj)
  if (loo && requireNamespace("loo", quietly = TRUE) &&
      all(!vapply(loo_objs, is.null, logical(1)))) {
    names(loo_objs) <- model
    loo_cmp <- loo::loo_compare(loo_objs)
  }

  structure(list(table = tab, lrt = lrt, loo_compare = loo_cmp,
                 has_loo = !is.null(loo_cmp)),
            class = "breedRB_anova")
}

#' @export
print.breedRB_anova <- function(x, ...) {
  cat("<breedRB_anova>  model comparison\n\n")
  print(x$table, row.names = FALSE)
  if (!is.null(x$loo_compare)) {
    cat("\nPSIS-LOO comparison (elpd_diff, se_diff; best model first):\n")
    print(round(x$loo_compare[, c("elpd_diff", "se_diff"), drop = FALSE], 2))
  }
  if (!is.null(x$lrt)) {
    cat("\nApproximate deviance ratio between consecutive models",
        "(HEURISTIC, not a\nformal test - boundary problem; prefer WAIC/LOO):\n")
    print(x$lrt, row.names = FALSE)
  }
  invisible(x)
}
