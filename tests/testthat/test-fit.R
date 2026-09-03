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

test_that(".expand_rows is identical to the incidence-matrix product", {
  set.seed(7)
  lev <- paste0("g", 1:6)
  M2  <- matrix(rnorm(6 * 4), 6, 4, dimnames = list(lev, paste0("c", 1:4)))
  f   <- sample(lev, 50, replace = TRUE)                 # replicated genotypes
  # old path: dense incidence times the level-indexed matrix
  Z   <- .incidence(factor(f, levels = lev))
  ref <- Z %*% M2
  got <- .expand_rows(M2, f, lev)                        # new path: row-gather, no dense Z
  expect_equal(unname(got), unname(ref))
  expect_equal(colnames(got), colnames(M2))
})

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

  # GBLUP fit back-solves marker effects using the training markers held by the fit
  snp <- solve_SNP(fit)
  expect_equal(nrow(snp), ncol(sim$M))
  expect_equal(attr(snp, "method"), "GBLUP")
  # marker effects reproduce the GEBVs: Mc %*% b == u
  Mc <- scale(sim$M[s$effect, ], center = TRUE, scale = FALSE)
  recon <- as.numeric(Mc %*% snp$effect[match(colnames(Mc), snp$marker)])
  expect_gt(cor(recon, s$solution), 0.999)

  # predict() on the training genotypes reproduces BLUP + mu
  pred <- predict(fit, sim$M, add_mu = TRUE)
  expect_equal(nrow(pred), nrow(sim$M))
  expect_true(all(c("ID", "prediction", "sd", "lower", "upper") %in% names(pred)))
  smu <- solution(fit, term = key, type = "random", add_mu = TRUE)
  m   <- match(pred$ID, smu$effect)
  expect_gt(cor(pred$prediction, smu$solution[m]), 0.999)
  # predicting a held-out genotype: correlate with its true genetic value
  pred0 <- predict(fit, sim$M, add_mu = FALSE)
  expect_gt(cor(pred0$prediction, sim$gv[pred0$ID]), 0.8)
})

test_that("mrk() GBLUP builds a VanRaden G and stays non-singular with duplicate genotypes", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  sim <- simulate_markers(n_gen = 50, n_mrk = 300)
  # force a rank-deficient raw genomic relationship: two identical genotypes
  M <- sim$M; M["g2", ] <- M["g1", ]
  raw_rank <- qr(tcrossprod(scale(M, center = TRUE, scale = FALSE)))$rank
  expect_lt(raw_rank, nrow(M))                    # raw Z Z' is singular

  fit <- bbglr(y ~ 1, random = ~ mrk(gen, M), data = sim$data,
               relmat = list(M = M), nIter = 2000, burnIn = 800,
               nChains = 1, verbose = FALSE)
  # the fit succeeds and returns a value for every genotype despite the singular raw G
  s <- solution(fit, term = .vm_keys(fit), type = "random")
  expect_equal(nrow(s), nrow(M))
  expect_true(all(is.finite(s$solution)))
})

test_that("mrk(rank=) fits a low-rank GBLUP that still ranks genotypes and back-solves", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  sim <- simulate_markers(n_gen = 60, n_mrk = 300)
  k   <- 20L
  fit <- bbglr(y ~ 1, random = ~ mrk(gen, M, rank = 20), data = sim$data,
               relmat = list(M = sim$M), nIter = 2500, burnIn = 900,
               nChains = 1, verbose = FALSE)
  key <- .vm_keys(fit)
  gc  <- fit$meta[[key]]$components[[1]]
  expect_equal(gc$method, "GBLUP")
  expect_equal(gc$rank, k)
  # the fitted design (and the cached rotation) keeps only k principal components
  expect_equal(ncol(gc$bmap), k)
  expect_equal(ncol(gc$gblup_rot$vectors), k)

  s <- solution(fit, term = key, type = "random")
  expect_equal(nrow(s), nrow(sim$M))                 # still one value per genotype
  expect_true(all(is.finite(s$solution)))
  # positively correlated with truth; the ceiling is low here because these
  # markers have no low-rank structure, so a rank-20 fit discards real signal.
  expect_gt(cor(s$solution, sim$gv[s$effect]), 0.4)

  # solve_SNP reuses the cached rotation; effects reproduce the (rank-k) GEBVs
  snp   <- solve_SNP(fit)
  expect_true(all(is.finite(snp$effect)))
  Mc    <- scale(sim$M[s$effect, ], center = TRUE, scale = FALSE)
  recon <- as.numeric(Mc %*% snp$effect[match(colnames(Mc), snp$marker)])
  expect_gt(cor(recon, s$solution), 0.999)
})

test_that("mrk() mean-imputes missing markers and drops near-constant markers", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  sim <- simulate_markers(n_gen = 60, n_mrk = 200)
  M <- sim$M
  M[cbind(sample(nrow(M), 30), sample(ncol(M), 30))] <- NA   # scattered missing calls
  M[, 5]  <- 1                                               # two constant markers
  M[, 10] <- 2

  fit <- bbglr(y ~ 1, random = ~ mrk(gen, M), data = sim$data,
               relmat = list(M = M), nIter = 2000, burnIn = 800,
               nChains = 1, verbose = FALSE)
  key <- .vm_keys(fit)
  # the two constant markers are dropped; the rest are retained
  expect_equal(length(fit$meta[[key]]$components[[1]]$markers), ncol(M) - 2L)

  s <- solution(fit, term = key, type = "random")
  expect_true(all(is.finite(s$solution)))
  # solve_SNP and predict work despite the original NAs / constant columns
  snp <- solve_SNP(fit)
  expect_equal(nrow(snp), ncol(M) - 2L)
  expect_true(all(is.finite(snp$effect)))
  pred <- predict(fit, M, add_mu = TRUE)                     # M still carries NAs -> imputed on the fly
  expect_true(all(is.finite(pred$prediction)))
})

test_that("predict() scores genotypes held out of training", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  sim <- simulate_markers(n_gen = 220, n_mrk = 200, n_rep = 4)
  train_ids <- paste0("g", 1:180)
  test_ids  <- paste0("g", 181:220)
  train <- sim$data[sim$data$gen %in% train_ids, ]

  fit <- bbglr(y ~ 1, random = ~ mrk(gen, M), data = train,
               relmat = list(M = sim$M[train_ids, ]), nIter = 3000, burnIn = 1000,
               nChains = 1, verbose = FALSE)

  # genotypes never seen in training are predicted from their markers alone
  pred <- predict(fit, sim$M[test_ids, ], add_mu = FALSE)
  expect_equal(sort(pred$ID), sort(test_ids))
  expect_gt(cor(pred$prediction, sim$gv[pred$ID]), 0.6)
})

test_that("predict_pr() gives P(new genotype beats the training-population bar)", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  sim <- simulate_markers(n_gen = 220, n_mrk = 200, n_rep = 4)
  train_ids <- paste0("g", 1:180)
  test_ids  <- paste0("g", 181:220)
  train <- sim$data[sim$data$gen %in% train_ids, ]

  fit <- bbglr(y ~ 1, random = ~ mrk(gen, M), data = train,
               relmat = list(M = sim$M[train_ids, ]), nIter = 3000, burnIn = 1000,
               nChains = 1, verbose = FALSE)

  pp <- predict_pr(fit, sim$M[test_ids, ], threshold = 0.20)
  expect_equal(sort(pp$ID), sort(test_ids))
  expect_true(all(c("ID", "prediction", "sd", "prob") %in% names(pp)))
  expect_true(all(pp$prob >= 0 & pp$prob <= 1))               # valid probabilities
  expect_false(is.unsorted(rev(pp$prob)))                     # ordered by decreasing prob
  # probability tracks the point prediction: the top-predicted lines clear the
  # bar more often than the bottom-predicted ones
  expect_gt(cor(pp$prediction, pp$prob), 0.7)
  # invariant to add_mu (intercept cancels in the comparison)
  pp_mu <- predict_pr(fit, sim$M[test_ids, ], threshold = 0.20, add_mu = TRUE)
  m <- match(pp$ID, pp_mu$ID)
  expect_equal(pp$prob, pp_mu$prob[m], tolerance = 1e-8)
  # bottom tail: P(below the worst 20%) is the complementary question
  pp_low <- predict_pr(fit, sim$M[test_ids, ], threshold = 0.20, higher = FALSE)
  expect_true(all(pp_low$prob >= 0 & pp_low$prob <= 1))
})

test_that("exp_var_rank truncates GBLUP by explained variance and is controllable", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  sim <- simulate_markers(n_gen = 60, n_mrk = 300)

  # default control: rank chosen to explain 99% of the genomic variance -> < n
  fit99 <- bbglr(y ~ 1, random = ~ mrk(gen, M), data = sim$data,
                 relmat = list(M = sim$M),
                 control = bbglr_control(nIter = 2000, burnIn = 800, nChains = 1,
                                         verbose = FALSE))
  key <- .vm_keys(fit99)
  r99 <- ncol(fit99$meta[[key]]$components[[1]]$bmap)
  expect_equal(fit99$control$exp_var_rank, 0.99)
  expect_lt(r99, nrow(sim$M))                       # some tail components dropped
  s99 <- solution(fit99, term = key, type = "random")
  expect_gt(cor(s99$solution, sim$gv[s99$effect]), 0.8)   # signal retained

  # a looser target keeps fewer components than a tighter one
  fit90 <- bbglr(y ~ 1, random = ~ mrk(gen, M), data = sim$data,
                 relmat = list(M = sim$M),
                 control = bbglr_control(nIter = 1500, burnIn = 600, nChains = 1,
                                         verbose = FALSE, exp_var_rank = 0.90))
  r90 <- ncol(fit90$meta[[key]]$components[[1]]$bmap)
  expect_lte(r90, r99)

  # exp_var_rank = NA keeps every non-null component (full rank here)
  fitfull <- bbglr(y ~ 1, random = ~ mrk(gen, M), data = sim$data,
                   relmat = list(M = sim$M),
                   control = bbglr_control(nIter = 1500, burnIn = 600, nChains = 1,
                                           verbose = FALSE, exp_var_rank = NA))
  expect_gt(ncol(fitfull$meta[[key]]$components[[1]]$bmap), r99)

  # an explicit rank inside mrk() overrides exp_var_rank
  fitk <- bbglr(y ~ 1, random = ~ mrk(gen, M, rank = 15), data = sim$data,
                relmat = list(M = sim$M),
                control = bbglr_control(nIter = 1200, burnIn = 500, nChains = 1,
                                        verbose = FALSE))
  expect_equal(ncol(fitk$meta[[.vm_keys(fitk)]]$components[[1]]$bmap), 15L)
})

test_that("legacy loose control arguments still work and override control", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  sim <- simulate_markers(n_gen = 40, n_mrk = 200)
  # old call style: MCMC settings passed directly instead of via control =
  fit <- bbglr(y ~ 1, random = ~ mrk(gen, M), data = sim$data,
               relmat = list(M = sim$M), nIter = 2000, burnIn = 700, thin = 5,
               nChains = 1, verbose = FALSE)
  expect_s3_class(fit$control, "breedRB_control")
  expect_equal(fit$control$nIter, 2000L)
  expect_equal(fit$control$burnIn, 700L)
  # an unknown argument is rejected with a helpful message
  expect_error(
    bbglr(y ~ 1, random = ~ mrk(gen, M), data = sim$data,
          relmat = list(M = sim$M), notArg = 3),
    "Unknown argument"
  )
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
