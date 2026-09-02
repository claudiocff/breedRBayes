#' Posterior probability summaries and rankings from extr_outs_bglr output
#'
#' Computes posterior-based probabilities and summaries for genotype performance,
#' stability, environment-specific performance, and Legendre random regression coefficients,
#' using the output produced by `extr_outs_bglr` (which itself summarizes a `FW_bglr` fit).
#'
#' Under the Finlay–Wilkinson formulation with Legendre basis and separate intercept:
#' #$# y_{ge} = \mu + l_e + \alpha_g + \sum_{d=1}^{D} \gamma_{g,d}\,\tilde P_d(x_e) + \varepsilon_{ge} #$#
#' this function provides probabilities and summaries for #$# \alpha_g #$# (genotypic intercepts),
#' the degree-specific Legendre coefficients #$# \gamma_{g,d} #$#, and their environment-specific predictions.
#'
#' Specifically, it computes:
#' - Across-environment performance probability: fraction of posterior samples in which a genotype’s intercept
#'   ranks in the top tail defined by `int` (or bottom tail if `increase = FALSE`).
#' - Pairwise performance probability: #$# P\left(\alpha_i > \alpha_j\right) #$# (or #$# < #$# if `increase = FALSE`).
#' - Across-environment stability probability: fraction of posterior samples in which the variance of the
#'   genotype-specific centered G×E predictions is below the `int`-quantile (i.e., more stable).
#' - Pairwise stability probability: #$# P\left(\mathrm{Var}_{e}(\text{GE}_i) < \mathrm{Var}_{e}(\text{GE}_j)\right) #$#.
#' - Joint probability: product of performance and stability probabilities (heuristic).
#' - Within-environment performance probability: fraction of posterior samples for each environment in which the
#'   genotype’s environment-specific prediction ranks in the top tail defined by `int` (or bottom tail if `increase = FALSE`).
#' - Posterior summaries (medians, SDs, 2.5% and 97.5% quantiles) for intercepts and all Legendre coefficients
#'   across degrees (returned in `hpd_coef`).
#' - Probability of being better than the best check (based on posterior draws).
#' - An approximate t-test-style comparison versus the best check using posterior SDs and observed counts.
#'
#' Note: Probabilities related to slope being greater than 1 (or less than 1) are not computed in this version.
#'
#' @param extr An object of class `extr_bglr` created by `extr_outs_bglr`, containing posterior chains,
#'   observation counts, gradient matrices (`gradient`, `gradient_basis`), genotype parameters, and check IDs.
#' @param int Numeric in (0,1). Intensity parameter defining the tail probability used for "top/bottom"
#'   performance thresholds (e.g., 0.2 means top 20%). Default = 0.2.
#' @param increase Logical. If TRUE, higher values of the response are better (e.g., yield);
#'   if FALSE, lower values are better (e.g., disease scores). Default = TRUE.
#' @param verbose Logical. Currently unused; reserved for future messages. Default = FALSE.
#'
#' @details
#' Inputs expected from `extr_outs_bglr`:
#' - `post$chain_int` (matrix): posterior samples of genotype intercepts (# genotypes × iterations).
#' - `post$chain_coef` (list): posterior samples per Legendre degree d = 1..D (each # genotypes × iterations).
#' - `post$chain_ge` (matrix): centered genotype-by-environment posterior predictions (combos × iterations).
#' - `post$chain_gge` (matrix): uncentered genotype-by-environment posterior predictions (combos × iterations).
#' - `gradient` (matrix): two columns [1, l], with an intercept of 1 and an environment index l (legacy).
#' - `gradient_basis` (matrix): Legendre basis on an adaptive grid (columns `deg0..degD`, with `deg0=1`).
#' - `G_param` (data frame): genotype-level posterior means for intercept and slopes `slope{d}`, plus a `check` flag.
#' - `N_obs` (data frame): counts of observations per genotype.
#' - `check` (character): vector of check genotype IDs.
#'
#' Interpretation notes:
#' - Performance probabilities are computed using global posterior thresholds (across all genotypes and samples),
#'   via empirical quantiles defined by `int`.
#' - Stability probabilities use the posterior distribution of the variance across environments for each genotype
#'   (computed per sample), comparing to the lower tail defined by `int`.
#' - The "best check" is defined by the highest posterior median among check genotypes; probabilities of being
#'   better than the best check are computed as the fraction of posterior samples in which each genotype’s
#'   intercept exceeds the best-check intercept sample.
#' - The t-test-style comparison uses posterior SDs (as SE proxies) and observed counts to construct a
#'   simple statistic; treat this as heuristic rather than a strict frequentist test.
#'
#' @return
#' An object of class `probsup_bglr` (a list) with:
#' - `across`:
#'   - `g_hpd`: posterior summaries for genotype intercepts (medians, quantiles, SD, probability vs best check, t-statistic, p-value, check flag, N).
#'   - `perfo`: across-environment performance probabilities (top tail by `int` or bottom if `increase = FALSE`).
#'   - `pair_perfo`: matrix of pairwise performance probabilities (#$# P(\alpha_i > \alpha_j) #$# or #$# < #$#).
#'   - `stabi`: across-environment stability probabilities (low-variance tail).
#'   - `pair_stabi`: matrix of pairwise stability probabilities (#$# P(\mathrm{Var}_e(\text{GE}_i) < \mathrm{Var}_e(\text{GE}_j)) #$#).
#'   - `best_probs`: probabilities of being better than the best check.
#'   - `joint_prob`: product of performance and stability probabilities.
#'   - `reg`: data frame with regression profiles over the environment gradient for each genotype (using `gradient_basis`).
#'   - `bayesian_stats`: consolidated table of main Bayesian summaries (intercept, degree-1 slope, top-prob, best-check-prob).
#'   - `hpd_coef`: long-format table with posterior summaries for all degrees (`slope1..slopeD`).
#' - `within`:
#'   - `perfo`: environment-specific performance probabilities for each genotype.
#'   - `hpd_gge`: posterior summaries for G×E predictions per environment-genotype combo.
#'   - `hpd_slope`: posterior summaries for degree-1 slopes (median, SD, 2.5% and 97.5% quantiles) for backward compatibility.
#'
#' @examples
#' \donttest{
#' # Assuming you've run FW_bglr and extracted with extr_outs_bglr:
#' # probs <- prob_sup_bglr(extr = ext, int = 0.2, increase = TRUE)
#' # names(probs)
#' # head(probs$across$g_hpd)
#' # probs$across$perfo[1:5, ]
#' # head(probs$across$hpd_coef)
#' }
#'
#' @seealso FW_bglr, extr_outs_bglr
#' @export
#' @importFrom dplyr arrange left_join mutate rename rowwise select group_by summarise filter reframe
#' @importFrom tidyr separate drop_na
#' @importFrom reshape2 melt
#' @importFrom stats pt quantile sd median var
prob_sup_bglr <- function(extr,
                          int = 0.2, 
                          increase = TRUE,
                          verbose = FALSE) {
  # Check input class
  if (!inherits(extr, "extr_bglr")) {
    stop("Input 'extr' must be of class 'extr_bglr'.")
  }
  check <- extr$check
  
  # Extract posterior chains and align as [samples x genotypes]
  post_g    <- t(extr$post$chain_int)
  post_ge   <- t(extr$post$chain_ge)
  post_gge  <- t(extr$post$chain_gge)
  
  # Names and dimensions
  gen_names   <- colnames(post_g)
  num_gen     <- length(gen_names)
  num_samples <- nrow(post_g)
  
  ge_cols   <- colnames(post_ge)
  ge_split  <- strsplit(ge_cols, ":")
  gen_ge    <- sapply(ge_split, `[`, 1)
  loc_ge    <- sapply(ge_split, `[`, 2)
  
  # Legacy gradient [1, l] and mean parameters
  X <- extr$gradient
  G <- extr$G_param
  
  # 1. Across-environment performance probabilities (intercept)
  perfo_across <- if (increase) {
    colMeans(post_g >= quantile(post_g, 1 - int))
  } else {
    colMeans(post_g <= quantile(post_g, 1 - int))
  }
  perfo_df <- data.frame(ID = gen_names, prob = perfo_across)
  perfo_df <- perfo_df[order(-perfo_df$prob), ]
  
  # 2. Pairwise performance probabilities
  pair_perfo <- matrix(NA, num_gen, num_gen, dimnames = list(gen_names, gen_names))
  for (i in seq_along(gen_names)) {
    for (j in seq_along(gen_names)) {
      pair_perfo[i, j] <- if (increase) {
        mean(post_g[, i] > post_g[, j])
      } else {
        mean(post_g[, i] < post_g[, j])
      }
    }
  }
  
  # 3. Across-environment stability: variance of centered GE per sample, per genotype
  var_ge <- sapply(gen_names, function(gen) {
    cols <- gen_ge == gen
    if (sum(cols) < 2) {
      warning("Genotype ", gen, " has fewer than 2 environments. Variance not computed.")
      return(rep(NA, num_samples))
    }
    apply(post_ge[, cols, drop = FALSE], 1, var)
  })
  stabi_across <- colMeans(var_ge <= quantile(var_ge, int, na.rm = TRUE))
  stabi_df <- data.frame(ID = gen_names, prob = stabi_across)
  stabi_df <- stabi_df[order(-stabi_df$prob), ]
  
  # 4. Pairwise stability probabilities
  pair_stabi <- matrix(NA, num_gen, num_gen, dimnames = list(gen_names, gen_names))
  for (i in seq_along(gen_names)) {
    for (j in seq_along(gen_names)) {
      pair_stabi[i, j] <- mean(var_ge[, i] < var_ge[, j], na.rm = TRUE)
    }
  }
  
  # 5. Joint probability (heuristic product)
  joint_prob <- perfo_across * stabi_across
  joint_df <- data.frame(ID = gen_names, prob = joint_prob)
  joint_df <- joint_df[order(-joint_df$prob), ]
  
  # 6. Within-environment performance probabilities (using uncentered G+GE)
  env <- unique(loc_ge)
  within_perfo_df <- data.frame()
  for (i in env) {
    gge_i <- post_gge[, grepl(i, colnames(post_gge))]
    within_perfo_i <- if (increase) {
      colMeans(gge_i >= quantile(gge_i, 1 - int))
    } else {
      colMeans(gge_i <= quantile(gge_i, 1 - int))
    }
    within_perfo_df_i <- data.frame(ID = gen_names, prob = within_perfo_i)
    within_perfo_df_i$env <- i
    within_perfo_df_i <- within_perfo_df_i[order(-within_perfo_df_i$prob), ]
    within_perfo_df <- rbind(within_perfo_df, within_perfo_df_i)
  }
  
  # Intercept HPD summaries
  hpd_g <- data.frame(
    ID = gen_names,
    Median  = apply(post_g, 2, median),
    HPD2.5  = apply(post_g, 2, quantile, probs = 0.025),
    HPD68   = apply(post_g, 2, quantile, probs = 0.68),
    HPD95   = apply(post_g, 2, quantile, probs = 0.95),
    HPD97.5 = apply(post_g, 2, quantile, probs = 0.975),
    SE      = apply(post_g, 2, sd),
    check   = gen_names %in% check
  )
  
  # 7. Probability of being better than the best check (intercepts)
  hpd_check <- hpd_g[hpd_g$ID %in% check, ]
  best <- hpd_check[order(hpd_check$Median, decreasing = TRUE), ]$ID[1]
  post_best <- post_g[, colnames(post_g) %in% best]
  prob_best <- colMeans(post_g > post_best)
  best_probs <- data.frame(ID = names(prob_best), prob = prob_best) |>
    dplyr::arrange(desc(prob))
  
  # Approximate t-style comparison vs best check
  N_obs <- extr$N_obs
  hpd_g <- dplyr::left_join(hpd_g, N_obs, by = "ID")
  blup_check <- hpd_g[hpd_g$ID == best, ]
  t_result <- hpd_g |> 
    dplyr::rowwise() |> 
    dplyr::mutate(
      diff    = Median - blup_check$Median,
      se_diff = sqrt((SE^2 / Obs) + (blup_check$SE^2 / blup_check$Obs)),
      t_stat  = diff / se_diff,
      df      = Obs + blup_check$Obs - 2,
      p_value = 2 * pt(-abs(t_stat), df = df)
    ) |> 
    dplyr::ungroup() |> 
    dplyr::select(ID, t_stat, p_value) |> 
    dplyr::arrange(p_value)
  
  hpd_g <- merge(hpd_g, best_probs, by = "ID") |> 
    dplyr::left_join(t_result, by = "ID") |> 
    dplyr::select(ID, Median, HPD2.5, HPD68, HPD95, HPD97.5, SE, prob, t_stat, p_value, check, Obs) |> 
    dplyr::mutate(p_value = round(p_value, 4)) |> 
    dplyr::arrange(desc(prob))
  
  # HPD summaries for G×E predictions (uncentered)
  hpd_gge <- data.frame(
    ID      = colnames(post_gge),
    Median  = apply(post_gge, 2, median),
    HPD2.5  = apply(post_gge, 2, quantile, probs = 0.025),
    HPD97.5 = apply(post_gge, 2, quantile, probs = 0.975)
  ) |> 
    tidyr::separate(col = ID, into = c("ENV", "ID"), sep = ":")
  
  # 8. Posterior summaries for all Legendre coefficients (degrees 1..D)
  hpd_coef <- data.frame()
  if (!is.null(extr$post$chain_coef) && length(extr$post$chain_coef) >= 1) {
    for (d in seq_along(extr$post$chain_coef)) {
      post_d <- t(extr$post$chain_coef[[d]])  # [samples x genotypes]
      df_d <- data.frame(
        ID      = colnames(post_d),
        degree  = paste0("slope", d),
        Median  = apply(post_d, 2, median),
        SE      = apply(post_d, 2, sd),
        HPD2.5  = apply(post_d, 2, quantile, probs = 0.025),
        HPD97.5 = apply(post_d, 2, quantile, probs = 0.975)
      )
      hpd_coef <- rbind(hpd_coef, df_d)
    }
  }
  
  # 9. Regression profiles over the environment gradient using Legendre basis
  # Build Beta means for deg0..degD (deg0 = intercept; degd = slope{d})
  D <- ncol(extr$gradient_basis) - 1
  slope_means <- NULL
  if (D >= 1) {
    slope_cols <- paste0("slope", 1:D)
    present <- intersect(slope_cols, colnames(G))
    slope_means <- as.matrix(G[, present, drop = FALSE])
    # Ensure ordering slope1..slopeD
    slope_means <- slope_means[, paste0("slope", 1:D), drop = FALSE]
  }
  Beta_means <- cbind(intercept = G$intercept, slope_means)
  colnames(Beta_means) <- paste0("deg", 0:D)
  # Predictions on the adaptive grid: [nGrid x nGen]
  Yhat_grid_gen <- extr$gradient_basis %*% t(Beta_means)
  reg <- reshape2::melt(Yhat_grid_gen)
  reg$gradient <- rep(extr$gradient[, 2], nrow(Beta_means)) # reuse legacy l-grid
  reg <- reg |>
    dplyr::select(-Var1) |>
    dplyr::rename(ID = Var2) |>
    dplyr::mutate(check = ID %in% check)
  
  # Consolidated Bayesian stats table (intercept, degree-1 slope, top-prob, best-check-prob)
  # Join degree-1 summaries if available
  if (nrow(hpd_coef) > 0 && "slope1" %in% hpd_coef$degree) {
    hpd_slope <- hpd_coef[hpd_coef$degree == "slope1", c("ID","Median","SE")]
  } else {
    hpd_slope <- data.frame(ID = gen_names, Median = NA_real_, SE = NA_real_)
  }
  
  bayes_pred <- hpd_g |> 
    dplyr::select(ID, check, Obs, Median, SE, prob) |> 
    dplyr::rename(prob_best = prob) |> 
    dplyr::left_join(perfo_df, by = "ID") |> 
    dplyr::rename(prob_top = prob) |> 
    dplyr::left_join(hpd_slope, by = "ID") |> 
    dplyr::rename(intercept   = Median.x,
                  sd_intercept = SE.x,
                  slope        = Median.y,
                  sd_slope     = SE.y) |> 
    dplyr::select(ID,
                  check,
                  Obs,
                  intercept,
                  sd_intercept,
                  slope,
                  sd_slope,
                  prob_top,
                  prob_best) |> 
    dplyr::arrange(-intercept)
  
  # Output structure (probabilities of slope > 1 removed)
  structure(
    list(
      across = list(
        g_hpd          = hpd_g,
        perfo          = perfo_df,
        pair_perfo     = pair_perfo,
        stabi          = stabi_df,
        pair_stabi     = pair_stabi,
        best_probs     = best_probs,
        joint_prob     = joint_df,
        reg            = reg,
        bayesian_stats = bayes_pred,
        hpd_coef       = hpd_coef
      ),
      within = list(
        perfo    = within_perfo_df,
        hpd_gge  = hpd_gge,
        # For backward compatibility: keep degree-1 HPD summary without slope>1 probs
        hpd_slope = if (nrow(hpd_coef) > 0) hpd_coef[hpd_coef$degree == "slope1", ] else NULL
      )
    ),
    class = "probsup_bglr",
    control = list(
      intensity = int,
      increase  = increase
    )
  )
}
