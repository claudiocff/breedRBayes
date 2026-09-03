# End-to-end fit on simulated data with a known genetic-variance ratio.
# Kept small; skipped on CRAN.

simulate_gblup <- function(n_gen = 60, n_rep = 6, h2 = 0.5, seed = 1) {
  set.seed(seed)
  M <- matrix(rbinom(n_gen * 400, 2, 0.3), n_gen, 400)
  M <- scale(M)
  G <- tcrossprod(M) / ncol(M)
  diag(G) <- diag(G) + 1e-4
  dimnames(G) <- list(paste0("g", 1:n_gen), paste0("g", 1:n_gen))
  vg <- 10; ve <- vg * (1 - h2) / h2
  u  <- crossprod(chol(G), rnorm(n_gen, 0, sqrt(vg)))
  ids <- rep(rownames(G), each = n_rep)
  y   <- u[ids, 1] + rnorm(length(ids), 0, sqrt(ve))
  list(data = data.frame(gen = ids, y = as.numeric(y)), G = G, h2 = h2)
}

test_that("bbglr fits GBLUP and recovers heritability", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  sim <- simulate_gblup(h2 = 0.5)
  td <- tempfile(); dir.create(td)
  fit <- bbglr(y ~ 1, random = ~ vm(gen, G), data = sim$data,
               relmat = list(G = sim$G), nIter = 4000, burnIn = 1500, thin = 5,
               nChains = 2, verbose = FALSE, saveAt = td)

  expect_s3_class(fit, "breedRB_fit")

  vc <- varcomp(fit)
  expect_true(all(c("vm_gen__G_", "varE") %in% vc$term))

  h2 <- heritability(fit)
  expect_s3_class(h2, "breedRB_h2")
  # posterior credible interval should bracket the true h2
  expect_gt(h2$summary$upper, 0.5)
  expect_lt(h2$summary$lower, 0.5)

  d <- mcmc_diag(fit)
  expect_true(all(c("n_eff", "Rhat") %in% names(d)))
  expect_true(all(d$Rhat < 1.3, na.rm = TRUE))

  s <- solution(fit, term = "vm(gen, G)", type = "random")
  expect_equal(nrow(s), nrow(sim$G))
  expect_true(all(c("effect", "solution", "sd", "lower", "upper") %in% names(s)))

  # add_mu shifts every solution onto the intercept scale by a constant mean
  smu <- solution(fit, term = "vm(gen, G)", type = "random", add_mu = TRUE)
  shift <- smu$solution - s$solution[match(smu$effect, s$effect)]
  expect_equal(shift, rep(shift[1], nrow(smu)), tolerance = 1e-6)
  expect_gt(abs(shift[1]), 0)

  # gebv() still works but is deprecated and delegates to solution()
  g <- suppressWarnings(gebv(fit))
  expect_warning(gebv(fit), "deprecated")
  expect_equal(nrow(g), nrow(sim$G))
  expect_true(all(c("ID", "gebv", "sd") %in% names(g)))
})
