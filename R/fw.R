#' Bayesian Multi-Environment Trial Analysis with BGLR via Random Regression (Finlay–Wilkinson)
#'
#' Fits the Finlay–Wilkinson model in two steps using the BGLR package, with genotype
#' effects through a genomic relationship matrix, environment main effects, and a random-regression
#' slope per genotype based on the environment index estimated in Step 1.
#'
#' Conceptually, the model can be written as:
#' \eqn{y_{ge} = \mu + l_e + \alpha_g + \beta_g \, l_e + \varepsilon_{ge}}
#' where y is the phenotype for genotype g in environment e, l_e is the environment index,
#' alpha_g is the genotype-specific intercept, beta_g is the genotype-specific slope
#' (sensitivity to environment), and epsilon is the residual.
#'
#' @param data A data frame containing the multi-environment trial data.
#' @param gen Character. Column name for genotype IDs.
#' @param env Character. Column name for environment IDs.
#' @param trial Character (optional). Column name for experiment/trial factor (if multiple trials within environments).
#' @param block Character (optional). Column name for replication/block factor.
#' @param check Character (optional). Column name with logical values (TRUE/FALSE) marking check genotypes.
#' @param trait Character. Column name for the response variable (phenotype).
#' @param res.het Logical. Heterogeneous residuals? (Present for compatibility; not used directly). Default = FALSE.
#' @param nIter Integer. Number of MCMC iterations for BGLR. Default = 2000.
#' @param burnIn Integer. Burn-in period. Default = 500.
#' @param verbose Logical. Print progress from BGLR? Default = TRUE.
#' @param saveAt Character. Path/prefix to save binary files from BGLR (used later by `readBinMat`). Default = "" (current working directory).
#' @param kinship.matrix Matrix (required). Genomic relationship matrix with row/column names matching genotype IDs in `data[[gen]]`.
#' @param seed Integer. Random seed for reproducibility. Default = 123.
#'
#' @details
#' Workflow:
#' - Step 1: Fit genotype (G) and environment (E) main effects via BGLR; obtain the environment projections to build an environmental index.
#' - Step 2: Fit genotype-specific intercept and slope (G + G×E) using the index derived from Step 1 through a random-regression formulation.
#'
#' Implementation notes:
#' - The genomic relationship matrix is aligned internally to the genotype incidence matrix (Zg) order, ensuring consistency between model inputs.
#' - The environment design matrix is centered to represent the environmental index used in the slope computation.
#' - Posterior samples for intercept and slope are read from the binary files saved by BGLR (`ETA_int_b.bin`, `ETA_slope_b.bin`), which depend on `saveAt`.
#' - BLAS threading is controlled via RhpcBLASctl for reproducibility/performance.
#'
#' Checks performed by the code:
#' - Verifies that specified columns (`gen`, `env`, `trait`, and optionally `trial`, `block`, `check`) exist in `data`.
#' - Filters `data` to genotypes present in `rownames(kinship.matrix)`.
#' - Aligns `kinship.matrix` to match the genotype incidence matrix (`Zg`) ordering/content.
#'
#' @return
#' An object of class `bglr_met` containing:
#' - `path`: the `saveAt` path used to write/read BGLR binaries.
#' - `GGE`: wide matrix (by genotype:environment rows) with raw predictions (G + GE) across iterations.
#' - `GE`: wide matrix (by genotype:environment rows) with centered predictions (mu + XB) across iterations.
#' - `INT`: matrix of genotype intercepts by iterations.
#' - `SLOPE`: matrix of genotype slopes by iterations.
#' - `fit1`: BGLR fit object for Step 1 (G + E).
#' - `fit2`: BGLR fit object for Step 2 (G + GE).
#' - `data`: curated/aligned data used in the model.
#' - `trait`, `gen`, `env`, `trial`, `block`: column names used in the analysis.
#' - `check`: vector of check genotype IDs (if `check` column was provided).
#'
#' @examples
#' \donttest{
#' # Example using EnvRtype data
#' library(EnvRtype)
#' data("maizeG")     # genomic relationship matrix
#' data("maizeYield") # phenotypic data
#'
#' data <- maizeYield
#' kinship.matrix <- maizeG
#'
#' # Select 5 random genotypes as checks (as in your setup)
#' set.seed(123)
#' check_ids <- as.character(sample(unique(data$gid), 5))
#' data$check <- data$gid %in% check_ids
#'
#' # Identify a candidate trait column (exclude ID-like columns)
#' candidate_traits <- setdiff(names(data), c("gid","env","trial","block","check"))
#' trait_col <- candidate_traits[1] # use the first available trait column
#'
#' # Create a temporary output directory for BGLR binary files
#' out_dir <- tempfile(pattern = "bglr_fw_")
#' dir.create(out_dir)
#'
#' # Run the Finlay–Wilkinson two-step model
#' fit_fw <- FW_bglr(
#'   data = data,
#'   gen = "gid",
#'   env = "env",
#'   trait = trait_col,
#'   check = "check",
#'   kinship.matrix = kinship.matrix,
#'   nIter = 2000,
#'   burnIn = 500,
#'   saveAt = paste0(out_dir, "/"),
#'   verbose = TRUE
#' )
#'
#' # Inspect results
#' names(fit_fw)
#' dim(fit_fw$INT)   # intercepts by genotype x iterations
#' dim(fit_fw$SLOPE) # slopes by genotype x iterations
#' dim(fit_fw$GE)    # genotype:environment rows x iterations
#' dim(fit_fw$GGE)   # same, uncentered
#'
#' # Check genotypes
#' fit_fw$check
#' }
#'
#' @references
#' Finlay, K. W., & Wilkinson, G. N. (1963). The analysis of adaptation in a plant-breeding programme.
#' Australian Journal of Agricultural Research, 14(6), 742–754.
#'
#' @export
#' @importFrom BGLR BGLR Multitrait readBinMat
#' @importFrom tidyr expand_grid pivot_wider
#' @importFrom dplyr group_by mutate ungroup left_join select
#' @importFrom stats model.matrix
#' @importFrom reshape2 melt
#' @importFrom tibble column_to_rownames
FW_bglr <- function(data, # data frame containing the data
                    gen, # column name for genotypes
                    env, # column name for environments
                    trial = NULL, # optional, column name for trials
                    block = NULL, # optional, column name for replications
                    check = NULL, # optional, column containing TRUE or FALSE for checks
                    trait, # column name for the response variable
                    res.het = FALSE, # boolean for heterogeneous residuals
                    nIter = 2000, # number of iterations for Bayesian sampling
                    burnIn = 500, # number of burn-in iterations
                    verbose = TRUE, # boolean to print progress
                    saveAt = "", # path to save results
                    kinship.matrix = NULL, # genomic relationship matrix (required for some models)
                    W = NULL, #matrix containing environmental covariates (env x ec)
                    order = 1, #polynomial order
                    seed = 123) {
  # ------------------------------------------------------------------------
  # Refactored Finlay-Wilkinson: the two model fits are now delegated to the
  # general `bbglr()` engine (which centralises the PC-rotation GBLUP trick,
  # Legendre bases, incidence construction and multi-chain machinery). This
  # function only orchestrates the two FW steps and reconstructs the
  # value-added `bglr_met` object consumed by extr_outs_bglr()/prob_sup_bglr().
  #
  # Step 1  y ~ 1,               random = ~ vm(gen, G) + env      -> env index
  # Step 2  y ~ 1 + leg(idx, o), random = ~ vm(gen, G)            (intercept)
  #                                       + leg(idx, o):vm(gen, G) (slope)
  #
  # The specialised environmental-covariate (`W`) index remains a bespoke
  # single BGLR fit, as bbglr() has no `W` special; only Step 1 branches on it.
  # NOTE: for order > 1 the slope degrees share a single Legendre:genomic BRR
  # term (one variance) rather than the legacy per-degree terms.
  # ------------------------------------------------------------------------
  set.seed(seed)
  RhpcBLASctl::blas_set_num_threads(1)
  data <- droplevels(as.data.frame(data))
  stopifnot(
    gen %in% colnames(data),
    env %in% colnames(data),
    trait %in% colnames(data),
    is.null(trial) || trial %in% colnames(data),
    is.null(block) || block %in% colnames(data),
    is.null(check) || check %in% colnames(data)
  )

  check_names <- if (!is.null(check)) unique(data[data[[check]] == TRUE, gen]) else character(0)

  keep <- c(gen, env, trial, block, check, trait)
  data <- data[, unique(keep[!vapply(keep, is.null, logical(1))]), drop = FALSE]
  data <- data[order(data[[env]], data[[gen]]), ]
  if (!is.null(kinship.matrix)) data <- data[data[[gen]] %in% rownames(kinship.matrix), ]
  data <- droplevels(data)

  gblup <- !is.null(kinship.matrix)
  relmat <- if (gblup) list(G = kinship.matrix) else list()

  # Nested trial/block become combined-level factors so bbglr's centered
  # incidence reproduces the legacy scale(model.matrix(~ -1 + factor(paste(...)))).
  nest_terms <- character(0)
  if (!is.null(trial)) {
    data$.fw_trial <- interaction(data[[env]], data[[trial]], drop = TRUE, sep = ":")
    nest_terms <- c(nest_terms, ".fw_trial")
  }
  if (!is.null(block)) {
    if (!is.null(trial)) {
      data$.fw_block <- interaction(data[[env]], data[[trial]], data[[block]], drop = TRUE, sep = ":")
    } else {
      data$.fw_block <- interaction(data[[env]], data[[block]], drop = TRUE, sep = ":")
    }
    nest_terms <- c(nest_terms, ".fw_block")
  }

  bt <- function(x) paste0("`", x, "`")
  gen_term <- if (gblup) sprintf("vm(%s, G)", bt(gen)) else bt(gen)
  nest_rhs <- if (length(nest_terms)) paste("+", paste(bt(nest_terms), collapse = " + ")) else ""

  # ---- Step 1: G + E to estimate the environment index -------------------
  cat(if (gblup) "kinship matrix provided: Running GBLUP\n" else "kinship matrix not provided: Running BLUP\n")
  cat("Running First Step (G + E)\n")

  if (!is.null(W)) {
    # Bespoke environmental-covariate path (PCA index), fitted directly with BGLR.
    cat("W provided: Environmental gradient predicted from environmental covariates\n")
    Xl <- .incidence(factor(data[[env]]))
    W  <- W[match(colnames(Xl), rownames(W)), , drop = FALSE]
    stopifnot(all(rownames(W) == colnames(Xl)))
    W_std   <- scale(W, center = TRUE, scale = TRUE)
    mu_W    <- attr(W_std, "scaled:center")
    sigma_W <- attr(W_std, "scaled:scale")
    pr   <- stats::prcomp(W_std, center = FALSE, scale. = FALSE)
    expl <- cumsum(pr$sdev^2) / sum(pr$sdev^2)
    k    <- which(expl >= 0.95)[1]; if (is.na(k) || k < 1) k <- min(3, ncol(W_std))
    S     <- pr$x[, 1:k, drop = FALSE]
    S_std <- scale(S, center = TRUE, scale = TRUE)
    cS <- attr(S_std, "scaled:center"); sS <- attr(S_std, "scaled:scale")
    step1_dir <- tempfile("fw_step1_"); dir.create(step1_dir)
    fit1 <- bbglr(stats::as.formula(sprintf("%s ~ 1", bt(trait))),
                  random = stats::as.formula(sprintf("~ %s%s", gen_term, nest_rhs)),
                  data = data, relmat = relmat, nIter = nIter, burnIn = burnIn,
                  nChains = 1, seed = seed, saveAt = step1_dir, verbose = verbose)
    # The environment index comes from a separate BRR on the standardised PCA scores.
    Wfit <- BGLR::BGLR(y = data[[trait]],
                       ETA = list(l = list(X = Xl %*% S_std, model = "BRR", saveEffects = TRUE)),
                       nIter = nIter, burnIn = burnIn,
                       saveAt = file.path(step1_dir, "wpath_"), verbose = FALSE)
    beta_w <- Wfit$ETA[["l"]]$b
    lHat1  <- as.vector(S_std %*% beta_w); names(lHat1) <- rownames(W)
    l_var  <- list(varB = Wfit$ETA[["l"]]$varB, SD.varB = Wfit$ETA[["l"]]$SD.varB)
    l_path <- file.path(step1_dir, "wpath_"); l_key <- "l"; fit1_obj <- Wfit
  } else {
    cat("W not provided: Environmental gradient based on average performance\n")
    step1_dir <- tempfile("fw_step1_"); dir.create(step1_dir)
    fit1 <- bbglr(stats::as.formula(sprintf("%s ~ 1", bt(trait))),
                  random = stats::as.formula(sprintf("~ %s + %s%s", gen_term, bt(env), nest_rhs)),
                  data = data, relmat = relmat, nIter = nIter, burnIn = burnIn,
                  nChains = 1, seed = seed, saveAt = step1_dir, verbose = verbose)
    envkey <- .fw_find_key(fit1, function(m) m$role == "random" &&
                             length(m$components) == 1L &&
                             identical(m$components[[1]]$kind, "factor") &&
                             identical(m$components[[1]]$var, env))
    lHat1  <- fit1$chains[[1]]$ETA[[envkey]]$b
    if (is.null(names(lHat1))) names(lHat1) <- fit1$meta[[envkey]]$components[[1]]$levels
    beta_w <- NULL; W <- NULL
    l_var  <- list(varB = fit1$chains[[1]]$ETA[[envkey]]$varB,
                   SD.varB = fit1$chains[[1]]$ETA[[envkey]]$SD.varB)
    l_path <- fit1$paths[1]; l_key <- envkey; fit1_obj <- fit1$chains[[1]]
  }

  # Per-observation environment index for the random-regression covariate.
  data$.fw_index <- as.numeric(lHat1[as.character(data[[env]])])

  # ---- Step 2: intercept + Legendre random regression on the index -------
  cat("Running Second Step (G + GE)\n")
  step2_dir <- tempfile("fw_step2_"); dir.create(step2_dir)
  fixed2  <- stats::as.formula(sprintf("%s ~ 1 + leg(.fw_index, %d)", bt(trait), order))
  random2 <- stats::as.formula(sprintf("~ %s + leg(.fw_index, %d):%s%s",
                                        gen_term, order, gen_term, nest_rhs))
  fit2 <- bbglr(fixed2, random = random2, data = data, relmat = relmat,
                nIter = nIter, burnIn = burnIn, nChains = 1, seed = seed,
                saveAt = step2_dir, verbose = verbose)

  # Locate the engine keys for the intercept, slope and fixed-basis terms.
  intkey   <- .fw_find_key(fit2, function(m) m$role == "random" &&
                             !.fw_has(m, "leg") && (.fw_has(m, "vm") || .fw_has_factor(m, gen)))
  slopekey <- .fw_find_key(fit2, function(m) m$role == "random" &&
                             .fw_has(m, "leg") && (.fw_has(m, "vm") || .fw_has_factor(m, gen)))
  envfixkey <- .fw_find_key(fit2, function(m) m$role == "fixed" && .fw_has(m, "leg"))

  # Back-mapping matrix: PC rotation for GBLUP, identity for BLUP.
  int_comp <- .fw_pick(fit2$meta[[intkey]]$components, if (gblup) "vm" else "factor")
  if (gblup) {
    PC <- int_comp$pc; gid <- rownames(PC)
  } else {
    gid <- int_comp$levels; PC <- diag(length(gid)); rownames(PC) <- gid
  }
  nGen <- length(gid)

  prefix   <- fit2$paths[1]

  # Intercept: PC-back-mapped ridge coefficients plus the population mean.
  # The saved effect samples (.bin) contain only post-burn-in draws, whereas the
  # trace files (mu.dat, *_b.dat) also include thinned burn-in samples; align by
  # taking the tail of the traces to the number of effect draws.
  BInt   <- BGLR::readBinMat(paste0(prefix, "ETA_", intkey, "_b.bin"))   # [nDraws x nGen]
  nDraws <- nrow(BInt)
  mu_chain <- utils::tail(scan(paste0(prefix, "mu.dat"), quiet = TRUE), nDraws)
  # BGLR writes a header row of column names for FIXED-effect traces.
  fix_coef_chain <- as.matrix(utils::read.table(paste0(prefix, "ETA_", envfixkey, "_b.dat"),
                                                header = TRUE))
  fix_coef_chain <- fix_coef_chain[seq.int(nrow(fix_coef_chain) - nDraws + 1L, nrow(fix_coef_chain)), , drop = FALSE]
  colnames(fix_coef_chain) <- paste0("deg", seq_len(order))

  INT  <- tcrossprod(PC, BInt) + rep(mu_chain, each = nGen)            # [nGen x nDraws]
  rownames(INT) <- gid

  # Slope degrees from the fused Legendre:genomic interaction term. The
  # Khatri-Rao column layout follows the component order recorded in the
  # metadata: when the genomic factor is first, columns run gen-major
  # (gen1[deg1..degO], gen2[...], ...); otherwise degree-major.
  BSl <- BGLR::readBinMat(paste0(prefix, "ETA_", slopekey, "_b.bin"))  # [nDraws x order*nGen]
  slope_comps <- fit2$meta[[slopekey]]$components
  gen_first <- identical(slope_comps[[1]]$kind, "vm") || identical(slope_comps[[1]]$kind, "factor")
  COEF <- vector("list", order)
  for (d in seq_len(order)) {
    cols <- if (gen_first) seq.int(d, by = order, length.out = nGen) else ((d - 1L) * nGen + 1L):(d * nGen)
    Bd   <- BSl[, cols, drop = FALSE] + fix_coef_chain[, d]            # add population degree-d coef
    COEF[[d]] <- tcrossprod(PC, Bd)
    rownames(COEF[[d]]) <- gid
  }

  # Environment-level Legendre gradient basis (deg0 forced to the constant 1).
  rng_env <- range(lHat1, na.rm = TRUE)
  x_env <- if (diff(rng_env) > 0) 2 * (lHat1 - rng_env[1]) / diff(rng_env) - 1 else rep(0, length(lHat1))
  X_gradient <- legendre_basis(x_env, order, orthonormal = TRUE)
  rownames(X_gradient) <- names(lHat1)
  colnames(X_gradient) <- paste0("deg", 0:order)
  X_gradient[, "deg0"] <- 1

  # G + GE predictions (centered GE and uncentered GGE), wide by "gen:env".
  GE <- list(); GGE <- list()
  for (i in seq_len(ncol(INT))) {
    Beta_i <- cbind(INT[, i], do.call(cbind, lapply(COEF, function(M) M[, i])))
    rownames(Beta_i) <- gid
    GE_i <- scale(X_gradient %*% t(Beta_i), center = TRUE, scale = FALSE)
    GE_i <- scale(t(GE_i), center = TRUE, scale = FALSE)
    GE_i <- reshape2::melt(GE_i) |>
      dplyr::mutate(comb = paste(Var1, Var2, sep = ":")) |>
      dplyr::select(-c(Var1, Var2)) |>
      dplyr::mutate(iter = i) |>
      tidyr::pivot_wider(names_from = iter, values_from = value) |>
      tibble::column_to_rownames("comb")
    GE[[i]] <- GE_i
    GGE_i <- X_gradient %*% t(Beta_i)
    GGE_i <- reshape2::melt(GGE_i) |>
      dplyr::mutate(comb = paste(Var1, Var2, sep = ":")) |>
      dplyr::select(-c(Var1, Var2)) |>
      dplyr::mutate(iter = i) |>
      tidyr::pivot_wider(names_from = iter, values_from = value) |>
      tibble::column_to_rownames("comb")
    GGE[[i]] <- GGE_i
  }
  GE  <- do.call("cbind", GE)
  GGE <- do.call("cbind", GGE)

  structure(
    list(
      path  = prefix,       # Step-2 chain prefix (int/slope/mu/varB traces)
      path1 = l_path,       # Step-1 prefix (environment-index variance trace)
      keys  = list(int = intkey, slope = slopekey, envfix = envfixkey, l = l_key),
      GGE = GGE, GE = GE, INT = INT,
      SLOPE = if (order >= 1) COEF[[1]] else NULL,
      COEF = COEF,
      fit1 = fit1_obj, fit2 = fit2$chains[[1]],
      l_var = l_var,
      data = data, trait = trait, gen = gen, env = env, trial = trial, block = block,
      check = check_names, beta_w = beta_w, obs = "obs",
      W = W, order = order, X_gradient = X_gradient, range = rng_env,
      pca = if (!is.null(beta_w)) list(
        rotation = pr$rotation[, 1:k, drop = FALSE], k = k,
        mu_W = mu_W, sigma_W = sigma_W, score_center = cS, score_scale = sS
      ) else NULL
    ),
    class = "bglr_met"
  )
}

#' Find the first ETA key of a `breedRB_fit` whose metadata satisfies `pred`
#' @keywords internal
.fw_find_key <- function(fit, pred) {
  hit <- Filter(function(k) isTRUE(pred(fit$meta[[k]])), names(fit$meta))
  if (!length(hit)) stop("FW_bglr: could not locate a required model term via bbglr metadata.",
                         call. = FALSE)
  hit[[1]]
}

#' Does a term's metadata contain a component of the given kind?
#' @keywords internal
.fw_has <- function(m, kind) any(vapply(m$components, function(c) identical(c$kind, kind), logical(1)))

#' Does a term contain a plain factor component on variable `var`?
#' @keywords internal
.fw_has_factor <- function(m, var) any(vapply(m$components,
  function(c) identical(c$kind, "factor") && identical(c$var, var), logical(1)))

#' Pick the first component of a term with the given kind
#' @keywords internal
.fw_pick <- function(components, kind) {
  hit <- Filter(function(c) identical(c$kind, kind), components)
  if (!length(hit)) stop("FW_bglr: expected a '", kind, "' component but found none.", call. = FALSE)
  hit[[1]]
}

#' Predict Environment Gradient (Index) for New Environments via PCA Projection
#'
#' Uses the PCA mapping and coefficients learned in `FW_bglr` (when environmental covariates `W`
#' were provided) to predict the environment gradient/index \eqn{l_e} for a new set of environments.
#'
#' Procedure:
#' 1) Standardize the new `Wp` covariate matrix with the original training means and sds.
#' 2) Project standardized `Wp` onto the retained PCA loadings (scores).
#' 3) Standardize the new scores using the training score center/scale.
#' 4) Multiply by the learned regression coefficients \eqn{\beta} from Step 1 to obtain \eqn{l_e}.
#'
#' @param model A fitted object from `FW_bglr` with `pca` and `beta_w` components populated
#'   (i.e., `W` was provided during training).
#' @param Wp Numeric matrix of environmental covariates for new environments (rows = environments,
#'   columns = same covariates used in training). Row names should be the environment IDs.
#'
#' @return
#' A named numeric vector of predicted environment gradient values \eqn{l_e}, with names taken from
#' `rownames(Wp)`.
#'
#' @examples
#' # Suppose fit_fw was obtained with W provided:
#' # new_W is a matrix with the same columns as fit_fw$W and rownames as environment IDs
#' # l_new <- predict_env_gradient_pca(fit_fw, new_W)
#'
#' @export
predict_env_gradient_pca <- function(model, Wp) {
  # Ensure the model has PCA info (only available when W was provided during training)
  stopifnot(!is.null(model$pca))
  rot <- model$pca$rotation
  k   <- model$pca$k
  muW <- model$pca$mu_W
  sdW <- model$pca$sigma_W
  cS  <- model$pca$score_center
  sS  <- model$pca$score_scale
  beta <- model$beta_w
  
  
  # Standardize using training means (muW) and sds (sdW) to ensure compatibility
  Wp_std <- scale(Wp, center = muW, scale = sdW)
  
  
  # Project standardized covariates onto retained PCA loadings
  S_new <- as.matrix(Wp_std) %*% rot   # [env_new x k]
  
  
  # Standardize new scores using the same center/scale used during training PCA scores
  S_new_std <- scale(S_new, center = cS, scale = sS)
  
  
  # Environment gradient prediction as linear combination of standardized scores
  l_pred <- as.vector(S_new_std %*% beta)
  names(l_pred) <- rownames(Wp)
  l_pred
}