# Finlay-Wilkinson wrapper, now delegating both fits to the bbglr() engine.
# Simulate a MET with genotype-specific intercepts (a) and sensitivities (b)
# to a linear environment gradient, and check the two-step recovery + the
# downstream extraction / probability layer.

simulate_met <- function(n_gen = 30, n_env = 8, seed = 7) {
  set.seed(seed)
  M <- scale(matrix(rbinom(n_gen * 300, 2, 0.3), n_gen, 300))
  G <- tcrossprod(M) / ncol(M); diag(G) <- diag(G) + 1e-4
  dimnames(G) <- list(paste0("g", seq_len(n_gen)), paste0("g", seq_len(n_gen)))
  a <- crossprod(chol(G), rnorm(n_gen, 0, 2))[, 1]
  b <- crossprod(chol(G), rnorm(n_gen, 0, 0.6))[, 1]
  envs    <- paste0("E", seq_len(n_env))
  env_idx <- stats::setNames(scale(seq_len(n_env))[, 1] * 2, envs)
  d <- expand.grid(gen = rownames(G), env = envs, stringsAsFactors = FALSE)
  d$y <- 50 + env_idx[d$env] * 3 + a[d$gen] + b[d$gen] * env_idx[d$env] +
    rnorm(nrow(d), 0, 1)
  d$check <- d$gen %in% paste0("g", 1:3)
  list(data = d, G = G, a = a, b = b)
}

test_that("FW_bglr fits via bbglr and recovers genotype intercepts", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  sim <- simulate_met()
  fit <- FW_bglr(sim$data, gen = "gen", env = "env", trait = "y", check = "check",
                 kinship.matrix = sim$G, nIter = 3000, burnIn = 1000, order = 1,
                 verbose = FALSE)

  expect_s3_class(fit, "bglr_met")
  # engine keys were resolved for every reconstructed quantity
  expect_true(all(c("int", "slope", "envfix", "l") %in% names(fit$keys)))
  expect_equal(nrow(fit$INT), nrow(sim$G))
  expect_equal(ncol(fit$INT), ncol(fit$SLOPE))
  expect_equal(nrow(fit$GE), nrow(sim$G) * length(unique(sim$data$env)))

  # intercepts are strongly recovered; slope signal is weaker but positive
  G_int   <- rowMeans(fit$INT)
  G_slope <- rowMeans(fit$SLOPE)
  expect_gt(cor(G_int[names(sim$a)], sim$a), 0.9)
  expect_gt(cor(G_slope[names(sim$b)], sim$b), 0.2)
})

test_that("extr_outs_bglr and prob_sup_bglr consume the refactored object", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  sim <- simulate_met()
  fit <- FW_bglr(sim$data, gen = "gen", env = "env", trait = "y", check = "check",
                 kinship.matrix = sim$G, nIter = 2000, burnIn = 600, order = 1,
                 verbose = FALSE)

  ext <- extr_outs_bglr(fit, gen = "gen")
  expect_s3_class(ext, "extr_bglr")
  expect_true(all(c("intercept", "l", "slope", "error") %in% ext$variances$effect))
  expect_true(all(ext$variances$var > 0))
  expect_true(all(c("mu", "intercept", "l", "slope1") %in% names(ext$diagnostic_trace)))

  ps <- prob_sup_bglr(ext, int = 0.2, increase = TRUE, verbose = FALSE)
  expect_true(all(c("across", "within") %in% names(ps)))
})
