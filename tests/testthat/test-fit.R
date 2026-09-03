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

  # pr(): probability of ranking in the top 20%
  p <- pr(fit, term = "vm(gen, G)", type = "random", threshold = 0.20)
  expect_equal(nrow(p), nrow(sim$G))
  expect_true(all(c("effect", "solution", "prob") %in% names(p)))
  expect_true(all(p$prob >= 0 & p$prob <= 1))
  # exactly k levels are flagged in every draw => probabilities sum to k
  expect_equal(sum(p$prob), attr(p, "k"), tolerance = 1e-8)
  # top-ranked genotype by posterior mean should have high top-20% probability
  expect_gt(p$prob[1], 0.5)

  # pr(pair = TRUE): pairwise P(A > B)
  pp <- pr(fit, term = "vm(gen, G)", type = "random", pair = TRUE)
  expect_equal(nrow(pp), choose(nrow(sim$G), 2))
  expect_true(all(c("A", "B", "prob") %in% names(pp)))
  expect_true(all(pp$prob >= 0 & pp$prob <= 1))

  # gebv() still works but is deprecated and delegates to solution()
  g <- suppressWarnings(gebv(fit))
  expect_warning(gebv(fit), "deprecated")
  expect_equal(nrow(g), nrow(sim$G))
  expect_true(all(c("ID", "gebv", "sd") %in% names(g)))
})

# ---------------------------------------------------------------------------
# mrk(): automatic GBLUP / RR-BLUP selection + solve_SNP() back-solve.
# ---------------------------------------------------------------------------

simulate_markers <- function(n_gen, n_mrk, n_rep = 3, seed = 1) {
  set.seed(seed)
  M <- matrix(rbinom(n_gen * n_mrk, 2, 0.3), n_gen, n_mrk)
  rownames(M) <- paste0("g", seq_len(n_gen))
  colnames(M) <- paste0("m", seq_len(n_mrk))
  Mc <- scale(M, center = TRUE, scale = FALSE)
  gv <- as.numeric(scale(as.numeric(Mc %*% rnorm(n_mrk))) * 3)
  names(gv) <- rownames(M)
  ids <- rep(rownames(M), each = n_rep)
  df  <- data.frame(gen = ids, y = gv[ids] + rnorm(length(ids), 0, 1))
  list(M = M, data = df, gv = gv)
}

test_that("mrk() auto-selects GBLUP when markers >= genotypes and solve_SNP back-solves", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  sim <- simulate_markers(n_gen = 40, n_mrk = 200)
  fit <- bbglr(y ~ 1, random = ~ mrk(gen, M), data = sim$data,
               relmat = list(M = sim$M), nIter = 3000, burnIn = 1000,
               nChains = 1, verbose = FALSE)
  key <- .vm_keys(fit)
  expect_length(key, 1L)
  expect_equal(fit$meta[[key]]$components[[1]]$method, "GBLUP")

  s <- solution(fit, term = key, type = "random")
  expect_equal(nrow(s), nrow(sim$M))
  expect_gt(cor(s$solution, sim$gv[s$effect]), 0.8)

  # GBLUP fit needs M to back-solve marker effects
  expect_error(solve_SNP(fit), "required")
  snp <- solve_SNP(fit, sim$M)
  expect_equal(nrow(snp), ncol(sim$M))
  expect_equal(attr(snp, "method"), "GBLUP")
  # marker effects reproduce the GEBVs: Mc %*% b == u
  Mc <- scale(sim$M[s$effect, ], center = TRUE, scale = FALSE)
  recon <- as.numeric(Mc %*% snp$effect[match(colnames(Mc), snp$marker)])
  expect_gt(cor(recon, s$solution), 0.999)
})

test_that("mrk() auto-selects RR-BLUP when genotypes > markers and solve_SNP reads fitted effects", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  sim <- simulate_markers(n_gen = 200, n_mrk = 60)
  fit <- bbglr(y ~ 1, random = ~ mrk(gen, M), data = sim$data,
               relmat = list(M = sim$M), nIter = 3000, burnIn = 1000,
               nChains = 1, verbose = FALSE)
  key <- .vm_keys(fit)
  expect_equal(fit$meta[[key]]$components[[1]]$method, "RRBLUP")

  s <- solution(fit, term = key, type = "random")
  expect_equal(nrow(s), nrow(sim$M))
  expect_gt(cor(s$solution, sim$gv[s$effect]), 0.8)

  # RR-BLUP estimates marker effects directly; M is not required
  snp <- solve_SNP(fit)
  expect_equal(nrow(snp), ncol(sim$M))
  expect_equal(attr(snp, "method"), "RRBLUP")
  Mc <- scale(sim$M[s$effect, ], center = TRUE, scale = FALSE)
  recon <- as.numeric(Mc %*% snp$effect[match(colnames(Mc), snp$marker)])
  expect_gt(cor(recon, s$solution), 0.999)
})
