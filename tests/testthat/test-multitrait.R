# Multi-trait and factor-analytic fits on simulated correlated traits.
# Three traits are used deliberately: the covariance-trace vech ordering only
# differs from a naive column-major layout when t >= 3, so this guards it.

simulate_mt <- function(n_gen = 40, n_rep = 4, seed = 11) {
  set.seed(seed)
  ids   <- rep(paste0("g", seq_len(n_gen)), each = n_rep)
  env   <- factor(rep(rep(c("E1", "E2"), each = n_gen * n_rep / 2)))
  g_eff <- rnorm(n_gen, 0, 2)[as.integer(factor(ids))]
  data.frame(
    gen   = ids,
    env   = env,
    yield = 10 + g_eff       + rnorm(length(ids)),
    y2    =  5 + 0.9 * g_eff + rnorm(length(ids)),
    y3    = 20 + 0.5 * g_eff + rnorm(length(ids))
  )
}

test_that("factor-analytic multi-trait fit yields valid variances and h2", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  dat <- simulate_mt()
  td  <- tempfile(); dir.create(td)
  fit <- bbglr(cbind(yield, y2, y3) ~ env, random = ~ fa(gen, 1),
               data = dat, nIter = 2500, burnIn = 800, thin = 5,
               nChains = 1, verbose = FALSE, saveAt = td)
  expect_s3_class(fit, "breedRB_fit")
  expect_true(fit$response$multitrait)

  vc <- varcomp(fit)
  # every genetic and residual variance must be positive (guards vech ordering)
  vars <- vc[grepl("^genetic:|^residual:", vc$term), ]
  expect_true(all(vars$mean > 0))
  # genetic correlations must lie in [-1, 1]
  rg <- vc[grepl("^rg:", vc$term), ]
  expect_true(all(rg$mean >= -1 & rg$mean <= 1))
  # traits share a single latent factor -> strong positive rg
  expect_true(all(rg$mean > 0.5))

  h2 <- heritability(fit)
  expect_s3_class(h2, "breedRB_h2")
  expect_true(all(h2$summary$mean > 0 & h2$summary$mean < 1))
})

test_that("fa() requires a multi-trait (cbind) response", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  dat <- simulate_mt()
  expect_error(
    bbglr(yield ~ env, random = ~ fa(gen, 1), data = dat,
          nIter = 200, burnIn = 50, verbose = FALSE, saveAt = tempfile()),
    "multi-trait"
  )
})
