#' Adaptive environment grid and Legendre basis for profile evaluation
#'
#' Builds an adaptive grid on the environment index (l) and the corresponding
#' Legendre polynomial basis used to evaluate genotype-specific profiles.
#' The grid resolution is tied to the response scale, and the basis is
#' orthonormalized for degrees 1..D while forcing degree 0 to be 1 to match
#' the separate genotype intercept used in the Finlay–Wilkinson random regression.
#'
#' Conceptual background:
#' - Finlay–Wilkinson random regression with Legendre basis (intercept separated):
#'   \eqn{y_{ge} = \mu + l_e + \alpha_g + \sum_{d=1}^{D} \gamma_{g,d}\,\tilde P_d(x_e) + \varepsilon_{ge}}
#'   where \eqn{x_e} is the scaled environment index in [-1,1], \eqn{\tilde P_d} are orthonormal Legendre polynomials.
#' - Orthonormalization of Legendre polynomials on [-1,1]:
#'   \eqn{\tilde P_n(x) = \sqrt{\frac{2n+1}{2}}\,P_n(x), \quad \int_{-1}^{1}\tilde P_n(x)^2\,dx = 1}.
#'   For compatibility with a separate intercept, the degree-0 column is set to 1.
#'
#' @param model A `bglr_met` object returned by `FW_bglr`. Must contain:
#'   `data`, `trait`, `range` (the environment-index range from Step 1),
#'   and `order` (maximum Legendre degree used in Step 2).
#' @param min_points Integer. Minimum number of points in the grid (default 50).
#' @param step_fraction Numeric in (0,1). Fraction of the response scale used to
#'   set the maximum step size in the environment grid (default 0.01).
#' @param use_scale Character. Which response scale to use to set grid resolution:
#'   one of `c("range","sd","IQR")`. Default `"range"`.
#'
#' @return A list with:
#' - `grid`: numeric vector of environment index values (length L).
#' - `basis`: matrix of Legendre basis values with columns named `deg0..degD`,
#'   where `deg0` is constant 1 and `deg1..degD` are orthonormalized Legendre values.
#'
#' @details
#' The environment index grid is constructed over `model$range` with a step size
#' bounded by `resp_scale * step_fraction`, where `resp_scale` is computed from
#' the response variable using the chosen scale (`range`, `sd`, or `IQR`).
#' The grid values are then scaled to [-1,1] and fed to `legendre_basis(...)`.
#' The first column (degree 0) is replaced by 1 to ensure that genotype intercepts
#' multiply a constant column in predictions.
#'
#' @examples
#' \donttest{
#' # Assuming you have fitted FW_bglr and obtained `fit_fw`:
#' # gb <- .make_gradient_basis(fit_fw, min_points = 100, step_fraction = 0.02, use_scale = "sd")
#' # str(gb$grid)
#' # str(gb$basis)
#' }
#'
#' @seealso FW_bglr, legendre_basis
#' @keywords internal
.make_gradient_basis <- function(model,
                                 min_points = 50,
                                 step_fraction = 0.01,
                                 use_scale = c("range","sd","IQR")) {
  use_scale <- match.arg(use_scale)
  
  # Response scale used to set the resolution of the grid
  y <- model$data[[model$trait]]
  resp_scale <- switch(use_scale,
                       range = diff(range(y, na.rm = TRUE)),
                       sd    = stats::sd(y, na.rm = TRUE),
                       IQR   = IQR(y, na.rm = TRUE))
  max_step <- max(1e-12, resp_scale * step_fraction)
  
  # Use the l-range from Step 1 saved in the model
  rng_env <- model$range
  L <- max(min_points, ceiling(diff(rng_env) / max_step) + 1)
  grid_l <- seq(rng_env[1], rng_env[2], length.out = L)
  
  # Scale grid_l to [-1, 1] for Legendre polynomials
  if (diff(rng_env) > 0) {
    x_grid <- 2 * (grid_l - rng_env[1]) / diff(rng_env) - 1
  } else {
    x_grid <- rep(0, length(grid_l))
  }
  
  # Orthonormal Legendre basis; set degree 0 column to 1 to match separate intercept
  P_grid <- legendre_basis(x_grid, order = model$order, orthonormal = TRUE)
  P_grid[, 1] <- 1
  colnames(P_grid) <- paste0("deg", 0:model$order)
  
  # Return both the l-grid and the corresponding basis (deg0..degD)
  list(grid = grid_l, basis = P_grid)
}


#' Extract posterior summaries, Legendre random regression coefficients,
#' variance components, diagnostics, and predictive profiles from FW_bglr
#'
#' Summarizes posterior chains, variance components, diagnostic traces, and
#' genotype-level Finlay–Wilkinson random regression parameters (intercept and
#' Legendre coefficients) from a model fitted by `FW_bglr`. Also constructs
#' environment-by-genotype predictions and per-genotype R^2 against observed means.
#'
#' Conceptual model with Legendre basis and separate intercept:
#' \eqn{y_{ge} = \mu + l_e + \alpha_g + \sum_{d=1}^{D} \gamma_{g,d}\,\tilde P_d(x_e) + \varepsilon_{ge}}
#' where \eqn{x_e} is the environment index scaled to [-1,1], \eqn{\tilde P_d} are
#' orthonormal Legendre polynomials for degrees \eqn{d=1,\dots,D}, and \eqn{\alpha_g}
#' is the genotype-specific intercept. The function collects posterior chains for
#' \eqn{\alpha_g} and \eqn{\gamma_{g,d}}, variance components, and diagnostics.
#'
#' @param model An object of class `bglr_met` returned by `FW_bglr`, including its
#'   `path` (saveAt), posterior effect matrices (`INT`, `COEF`, `SLOPE`, `GE`, `GGE`),
#'   and metadata (`data`, `trait`, `env`, `order`, `check`, `X_gradient`, `range`).
#' @param gen Character. Column name for genotype IDs, matching the column used in `FW_bglr`.
#'
#' @return An object of class `extr_bglr` (list) containing:
#' - `post`: list of posterior chains:
#'   - `chain_int` [nGen x nIter]: genotype intercepts.
#'   - `chain_coef`: list of degree-specific chains (`slope1..slopeD`), each [nGen x nIter].
#'   - `chain_slope`: first-degree chain for backward compatibility (may be NULL if `order==0`).
#'   - `chain_ge`: centered genotype-by-environment predictions across iterations.
#'   - `chain_gge`: uncentered genotype + GE predictions across iterations.
#' - `variances`: data frame with posterior means and SDs of variance components:
#'   intercept, environment main effect (`l`), each slope_d, and error (from Step 2).
#' - `diagnostic_trace`: data frame with traces for `mu`, `intercept`, `l`,
#'   each `slope{d}`, and iteration index.
#' - `effective_sizes`: numeric vector of effective sample sizes for diagnostic traces
#'   computed via `coda::effectiveSize`.
#' - `gradient`: matrix with two columns (intercept=1, gradient values) spanning the
#'   range of the environment effect (`l`) for legacy compatibility.
#' - `gradient_basis`: Legendre basis matrix (deg0..degD) on the adaptive grid returned by
#'   `.make_gradient_basis`, suitable for evaluating profile curves.
#' - `G_param`: data frame with genotype-level posterior means for `intercept` and
#'   `slope{d}` (degrees 1..D), plus a `check` flag.
#' - `N_obs`: data frame with counts of observations per genotype (`ID`, `Obs`).
#' - `check`: vector of check genotype IDs propagated from the model.
#' - `ID_R2`: data frame with per-genotype R^2 (squared correlation across environments
#'   between predicted and observed means).
#'
#' @details
#' Files read for diagnostics:
#' - `mu.dat`, `ETA_int_varB.dat`, `ETA_l_varB.dat`, and `ETA_slope{d}_varB.dat` (for each degree).
#' - Ensure `FW_bglr` was run with `saveAt` pointing to a writable directory that ends
#'   with a path separator (e.g., `"/"`), so these files exist at extraction time.
#'
#' Predictions and R^2:
#' - Environment-by-genotype predictions use the model's environmental basis `X_gradient`
#'   (columns `deg0..degD`, with `deg0=1`) and posterior mean coefficients
#'   (`intercept` and `slope{d}`).
#' - Observed means by genotype×environment are computed from `model$data`, and R^2 is
#'   the squared correlation between predictions and observed means across environments per genotype.
#'
#' @examples
#' \donttest{
#' # Fit FW_bglr first (see FW_bglr examples), then:
#' # ext <- extr_outs_bglr(model = fit_fw, gen = "gid")
#' # names(ext)
#' # head(ext$G_param)
#' # ext$variances
#' # ext$effective_sizes
#' # str(ext$gradient_basis)
#' }
#'
#' @seealso FW_bglr, .make_gradient_basis, legendre_basis
#' @export
#' @importFrom coda effectiveSize
#' @importFrom dplyr rename select group_by summarise filter reframe left_join
#' @importFrom tidyr drop_na
#' @importFrom reshape2 melt
#' @importFrom utils scan
extr_outs_bglr <- function(model, gen) {
  
  check <- model$check
  
  # Counts of observations per genotype
  N <- data.frame(table(model$data[[gen]])) |>
    dplyr::rename(ID = Var1, Obs = Freq)
  
  # Posterior chains from FW_bglr
  post <- list(
    chain_int   = model$INT,     # [nGen x nIter]
    chain_coef  = model$COEF,    # list: degree d -> [nGen x nIter]
    chain_slope = model$SLOPE,   # keeps slope1 for backward compatibility (may be NULL if order == 0)
    chain_ge    = model$GE,      # centered G×E
    chain_gge   = model$GGE      # uncentered G + GE
  )
  
  # Engine (bbglr) term keys, resolved at fit time by FW_bglr().
  keys <- model$keys

  # Variances: intercept (Step 2), l (Step-1 environment index), the fused
  # Legendre:genomic slope term (Step 2), and residual error (Step 2).
  variances_list <- list(
    data.frame(effect = "intercept",
               var = model$fit2$ETA[[keys$int]]$varB,
               sd  = model$fit2$ETA[[keys$int]]$SD.varB),
    data.frame(effect = "l",
               var = model$l_var$varB,
               sd  = model$l_var$SD.varB),
    data.frame(effect = "slope",
               var = model$fit2$ETA[[keys$slope]]$varB,
               sd  = model$fit2$ETA[[keys$slope]]$SD.varB),
    data.frame(effect = "error",
               var = model$fit2$varE,
               sd  = model$fit2$SD.varE)
  )
  variances <- do.call(rbind, variances_list)

  # Diagnostic variance-component traces read from the BGLR output files.
  # Step-2 traces live under `path`; the environment-index (l) trace under `path1`.
  path  <- model$path
  path1 <- if (!is.null(model$path1)) model$path1 else model$path
  mu_trace    <- scan(paste0(path,  "mu.dat"), quiet = TRUE)
  int_trace   <- scan(paste0(path,  "ETA_", keys$int,   "_varB.dat"), quiet = TRUE)
  slope_trace <- scan(paste0(path,  "ETA_", keys$slope, "_varB.dat"), quiet = TRUE)
  l_trace     <- scan(paste0(path1, "ETA_", keys$l,     "_varB.dat"), quiet = TRUE)
  n_iter <- length(mu_trace)
  # Align the (possibly separately-run) Step-1 l trace to the Step-2 length.
  if (length(l_trace) != n_iter) length(l_trace) <- n_iter

  # Assemble diagnostic data frame (mu, intercept, l, slope, and iteration)
  diagnostic <- data.frame(
    mu = mu_trace,
    intercept = int_trace,
    l = l_trace,
    slope1 = slope_trace,
    iter = seq_len(n_iter)
  )
  
  # Effective sample sizes for diagnostics (excluding the iteration index)
  effective_sizes <- coda::effectiveSize(diagnostic |> dplyr::select(-iter))
  
  # Gradient grid and Legendre basis for profiles (keep legacy gradient and add full basis)
  gb <- .make_gradient_basis(model, min_points = 50, step_fraction = 0.01, use_scale = "range")
  gradient_l     <- gb$grid
  gradient_basis <- gb$basis
  gradient       <- cbind(intercept = 1, gradient = gradient_l)  # legacy compatibility (intercept + l only)
  
  # Genotype-level parameters: posterior means for intercept and each slope_d
  G_int <- rowMeans(model$INT)
  G_coef_means <- list()
  if (model$order >= 1) {
    for (d in seq_len(model$order)) {
      G_coef_means[[paste0("slope", d)]] <- rowMeans(model$COEF[[d]])
    }
  }
  G_param <- data.frame(ID = names(G_int),
                        intercept = G_int,
                        check = names(G_int) %in% check)
  if (length(G_coef_means) > 0) {
    for (nm in names(G_coef_means)) {
      G_param[[nm]] <- G_coef_means[[nm]]
    }
    # Backward compatibility: slope = slope1
    if ("slope1" %in% names(G_param)) {
      G_param[["slope"]] <- G_param[["slope1"]]
    }
  } else {
    # Case order == 0 (no slopes)
    G_param[["slope"]] <- NA_real_
  }
  
  # Predictions by environment × genotype using the model's environmental basis
  # X_env: rows = environments, columns = deg0..degD (deg0 = 1)
  X_env <- model$X_gradient
  # Beta_means: columns deg0..degD (deg0 is intercept mean; deg1..degD are slope means)
  Beta_means <- cbind(
    deg0 = G_int,
    do.call(cbind, G_coef_means)
  )
  # Name columns to match deg0..degD (ordering is already deg0, slope1, ..., slopeD)
  col_order <- paste0("deg", 0:model$order)
  colnames(Beta_means) <- paste0("deg", 0:model$order)
  
  # Predictions matrix: [nEnv x nGen]
  Yhat_env_gen <- X_env %*% t(Beta_means)
  
  # Observed averages by environment × genotype for comparison
  obs_means <- model$data |>
    dplyr::group_by(.data[[gen]], .data[[model$env]]) |>
    dplyr::summarise(value = mean(.data[[model$trait]], na.rm = TRUE), .groups = "drop") |>
    dplyr::rename(ID = !!gen, env = !!model$env)
  
  # Stack predictions and join with observed means
  pred_data <- reshape2::melt(Yhat_env_gen) |>
    dplyr::rename(env = Var1, ID = Var2, yhat = value) |>
    dplyr::left_join(obs_means, by = c("ID", "env"))
  
  
  # Return structured list with chains, variances, diagnostics, gradients, parameters, and R^2
  structure(
    list(
      post = post,
      variances = variances,
      diagnostic_trace = diagnostic,
      effective_sizes = effective_sizes,
      gradient = gradient,               # legacy compatibility (intercept + l)
      gradient_basis = gradient_basis,   # new basis (deg0..degD) for profile evaluation
      G_param = G_param,
      N_obs = N,
      check = check
    ),
    class = "extr_bglr"
  )
}
