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
  expect_true(all(c("effect", "solution", "sd", "lower", "upper",
                    "pev", "reliability") %in% names(s)))
  # PEV is the posterior variance of the effect (= sd^2) and reliability is in [0, 1]
  expect_equal(s$pev, s$sd^2, tolerance = 1e-8)
  expect_true(all(s$reliability >= 0 & s$reliability <= 1))

  # add_mu shifts every solution onto the intercept scale by a constant mean
  smu <- solution(fit, term = "vm(gen, G)", type = "random", add_mu = TRUE)
  shift <- smu$solution - s$solution[match(smu$effect, s$effect)]
  expect_equal(shift, rep(shift[1], nrow(smu)), tolerance = 1e-6)
  expect_gt(abs(shift[1]), 0)
  # PEV / reliability are a property of the effect, unchanged by the mu shift
  m <- match(s$effect, smu$effect)
  expect_equal(s$pev, smu$pev[m], tolerance = 1e-8)
  expect_equal(s$reliability, smu$reliability[m], tolerance = 1e-8)

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

# ---------------------------------------------------------------------------
# Genomic random regression (reaction norms): per-genotype and per-marker
# intercept + slope, for both GBLUP and RR-BLUP genomic bases.
# ---------------------------------------------------------------------------

simulate_rr <- function(n_gen, n_mrk, n_env = 12, seed = 1) {
  set.seed(seed)
  M <- matrix(rbinom(n_gen * n_mrk, 2, 0.3), n_gen, n_mrk)
  rownames(M) <- paste0("g", seq_len(n_gen)); colnames(M) <- paste0("m", seq_len(n_mrk))
  Ms <- scale(M)
  u_int <- as.numeric(Ms %*% rnorm(n_mrk)); u_slp <- as.numeric(Ms %*% rnorm(n_mrk))
  names(u_int) <- names(u_slp) <- rownames(M)
  gg <- rep(rownames(M), each = n_env)
  xx <- rep(seq(-1, 1, length.out = n_env), n_gen)
  y  <- u_int[gg] + u_slp[gg] * xx + rnorm(n_gen * n_env, 0, 0.5)
  list(M = M, u_int = u_int, u_slp = u_slp,
       data = data.frame(gen = gg, x = xx, y = as.numeric(y)))
}

test_that("solution()/solve_SNP() give per-genotype and per-marker intercept + slope for a genomic RR", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  int_term <- "mrk(gen, M)"; slp_term <- "mrk(gen, M):leg(x, 1)"
  for (sim in list(GBLUP  = simulate_rr(n_gen = 50,  n_mrk = 300),   # markers >= genotypes
                   RRBLUP = simulate_rr(n_gen = 120, n_mrk = 60))) { # genotypes > markers
    td <- tempfile(); dir.create(td)
    # environmental covariate enters BOTH the fixed (population) and random parts
    fit <- bbglr(y ~ 1 + leg(x, 1), random = ~ mrk(gen, M) + mrk(gen, M):leg(x, 1),
                 data = sim$data, relmat = list(M = sim$M),
                 nIter = 4000, burnIn = 1500, thin = 5, nChains = 2,
                 verbose = FALSE, saveAt = td)

    # intercept from the main effect, slope from the interaction — one row per genotype
    ii <- solution(fit, term = int_term, type = "random")
    ss <- solution(fit, term = slp_term, type = "random")
    expect_equal(nrow(ii), nrow(sim$M))
    expect_equal(nrow(ss), nrow(sim$M))                  # leg order 1 => one slope per genotype
    expect_gt(cor(ii$solution, sim$u_int[ii$effect]), 0.9)
    expect_gt(cor(ss$solution, sim$u_slp[ss$effect]), 0.9)
    expect_true(all(c("pev", "reliability") %in% names(ss)))  # RR term still gets PEV/reliability

    # marker effects for intercept and slope, and the GEBV = Mc b round-trip per coef
    Mc <- scale(sim$M[ii$effect, ], center = TRUE, scale = FALSE)
    attr(Mc, "scaled:center") <- NULL
    for (tm in list(list(int_term, ii), list(slp_term, ss))) {
      snp <- solve_SNP(fit, term = tm[[1]])
      expect_equal(nrow(snp), ncol(sim$M))
      b <- snp$effect; names(b) <- snp$marker
      recon <- as.numeric(Mc[, names(b)] %*% b)
      geno  <- tm[[2]]
      expect_gt(cor(recon, geno$solution[match(ii$effect, geno$effect)]), 0.999)
    }
  }
})

test_that("an order>=2 random regression splits into one variance component per Legendre degree", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  set.seed(21); n_gen <- 40; n_mrk <- 250; n_env <- 10
  M <- matrix(rbinom(n_gen * n_mrk, 2, 0.3), n_gen, n_mrk)
  rownames(M) <- paste0("g", seq_len(n_gen)); colnames(M) <- paste0("m", seq_len(n_mrk))
  Ms  <- scale(M)
  u0  <- as.numeric(Ms %*% rnorm(n_mrk))          # intercept
  u1  <- as.numeric(Ms %*% rnorm(n_mrk)) * 1.0    # linear (larger variance)
  u2  <- as.numeric(Ms %*% rnorm(n_mrk)) * 0.4    # quadratic (smaller variance)
  names(u0) <- names(u1) <- names(u2) <- rownames(M)
  gg <- rep(rownames(M), each = n_env)
  xx <- rep(seq(-1, 1, length.out = n_env), n_gen)
  Bo <- legendre_basis(.scale_unit(xx, range(xx)), order = 2, orthonormal = TRUE)
  y  <- u0[gg] + u1[gg] * Bo[, 2] + u2[gg] * Bo[, 3] + rnorm(n_gen * n_env, 0, 0.5)
  dat <- data.frame(gen = gg, x = xx, y = as.numeric(y))

  td <- tempfile(); dir.create(td)
  fit <- bbglr(y ~ 1 + leg(x, 2), random = ~ mrk(gen, M) + mrk(gen, M):leg(x, 2),
               data = dat, relmat = list(M = M),
               nIter = 4000, burnIn = 1500, thin = 2, nChains = 1,
               verbose = FALSE, saveAt = td)

  key <- .resolve_term(fit, "mrk(gen, M):leg(x, 2)")
  meta <- fit$meta[[key]]
  # the RR term is now two BGLR blocks, one per degree, each with its own varB
  expect_true(isTRUE(meta$rr_split))
  expect_equal(length(meta$eta_keys), 2L)
  expect_match(meta$eta_keys, "deg[12]$", all = TRUE)

  vc <- varcomp(fit)
  # two independent variance components exist and are estimated separately
  vB <- vc$mean[match(meta$eta_keys, vc$term)]
  expect_equal(length(vB), 2L)
  expect_true(all(is.finite(vB)))
  expect_true(vB[1] != vB[2])                     # genuinely distinct, not one shared component
  expect_gt(vB[1], vB[2])                          # linear variance > quadratic (as simulated)

  # solution() reassembles the split blocks into the per-genotype deg1/deg2 layout
  s  <- solution(fit, term = key, type = "random")
  D  <- attr(s, "draws")
  expect_equal(ncol(D), 2L * n_gen)
  expect_true(all(grepl(":deg[12]$", colnames(D))))
  expect_setequal(unique(sub(".*:", "", colnames(D))), c("deg1", "deg2"))

  # per-coefficient heritability returns intercept + one h2 per degree
  h2 <- heritability(fit, genetic = key)
  expect_equal(nrow(h2$summary), 3L)

  # an order-1 RR is NOT split (identical to the single-block build)
  fit1 <- bbglr(y ~ 1 + leg(x, 1), random = ~ mrk(gen, M) + mrk(gen, M):leg(x, 1),
                data = dat, relmat = list(M = M),
                nIter = 2000, burnIn = 800, thin = 2, nChains = 1,
                verbose = FALSE, saveAt = tempfile())
  m1 <- fit1$meta[[.resolve_term(fit1, "mrk(gen, M):leg(x, 1)")]]
  expect_false(isTRUE(m1$rr_split))
  expect_equal(length(m1$eta_keys), 1L)
})

test_that("a solo fixed factor uses cell-means coding (every level shown); interactions keep contrasts", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  set.seed(3)
  env <- factor(rep(c("E1", "E2", "E3", "E4"), each = 60))
  rep_ <- factor(rep(rep(c("R1", "R2", "R3"), each = 20), 4))
  gen  <- factor(rep(paste0("g", 1:20), 12))
  mu_env <- c(E1 = 2, E2 = 5, E3 = 9, E4 = 4)
  ug     <- stats::setNames(rnorm(20, 0, 1), paste0("g", 1:20))
  y <- mu_env[as.character(env)] + ug[as.character(gen)] + rnorm(240, 0, 0.7)
  dat <- data.frame(env = env, rep = rep_, gen = gen, y = y)

  fit <- bbglr(y ~ env + env:rep, random = ~ gen, data = dat,
               nIter = 3000, burnIn = 1000, nChains = 1, verbose = FALSE)

  s <- solution(fit, term = "env", type = "fixed")
  expect_equal(nrow(s), nlevels(env))                     # ALL levels, none dropped
  expect_true(all(levels(env) %in% s$effect))             # incl. the old reference E1
  # add_mu recovers the per-level means (cell means), which track the truth
  sm <- solution(fit, term = "env", type = "fixed", add_mu = TRUE)
  cm <- stats::setNames(sm$solution, sm$effect)[names(mu_env)]
  expect_gt(cor(cm, mu_env), 0.99)

  # a fixed factor inside an interaction keeps treatment contrasts (identifiable)
  si <- solution(fit, term = "env:rep", type = "fixed")
  expect_equal(nrow(si), (nlevels(env) - 1L) * (nlevels(rep_) - 1L))
})

test_that("reaction_norm() evaluates per-genotype curves across the gradient", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  sim <- simulate_rr(n_gen = 50, n_mrk = 300, n_env = 12)
  fit <- bbglr(y ~ 1 + leg(x, 1), random = ~ mrk(gen, M) + mrk(gen, M):leg(x, 1),
               data = sim$data, relmat = list(M = sim$M),
               nIter = 4000, burnIn = 1500, thin = 5, nChains = 2, verbose = FALSE)

  rn <- reaction_norm(fit, term = "mrk(gen, M):leg(x, 1)", type = "random",
                      plot = FALSE, n_grid = 60)
  expect_true(all(c("id", "gradient", "value") %in% names(rn)))
  expect_equal(nrow(rn), nrow(sim$M) * 60L)                 # one row per genotype x grid
  expect_equal(length(unique(rn$id)), nrow(sim$M))
  expect_equal(range(rn$gradient), c(-1, 1))                 # Legendre [-1, 1] domain

  # per-genotype fitted slope (linear fit of the curve) recovers the true slope
  sl <- vapply(split(rn, rn$id),
               function(d) stats::coef(stats::lm(value ~ gradient, d))[2], numeric(1))
  expect_gt(cor(sl[names(sim$u_slp)], sim$u_slp), 0.9)

  # add_fixed_reg = FALSE returns genotype deviations centred on ~0 across the grid
  rnd <- reaction_norm(fit, term = "mrk(gen, M):leg(x, 1)", plot = FALSE,
                       add_fixed_reg = FALSE, n_grid = 60)
  expect_lt(abs(mean(rnd$value)), 1e-6)
  # the fixed regression shifts the curves by a non-zero population mean
  expect_gt(abs(mean(rn$value) - mean(rnd$value)), 0)

  # leg_basis = FALSE puts the gradient back on the original covariate scale
  rn2 <- reaction_norm(fit, term = "mrk(gen, M):leg(x, 1)", plot = FALSE,
                       leg_basis = FALSE, n_grid = 60)
  expect_equal(range(rn2$gradient), range(sim$data$x), tolerance = 1e-8)

  # a non-RR term is rejected; the plot object is returned when plot = TRUE
  expect_error(reaction_norm(fit, term = "mrk(gen, M)", plot = FALSE),
               "random-regression")
  rp <- reaction_norm(fit, term = "mrk(gen, M):leg(x, 1)", plot = TRUE, n_grid = 20)
  expect_s3_class(attr(rp, "plot"), "ggplot")
})

test_that("reaction_norm() works for a plain-factor (non-genomic) random regression", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  set.seed(2)
  n_gen <- 30; n_env <- 14
  gg <- rep(paste0("g", 1:n_gen), each = n_env)
  xx <- rep(seq(-2, 2, length.out = n_env), n_gen)
  a  <- stats::setNames(rnorm(n_gen, 5, 1),   paste0("g", 1:n_gen))
  b  <- stats::setNames(rnorm(n_gen, 0, 0.8), paste0("g", 1:n_gen))
  y  <- a[gg] + b[gg] * xx + rnorm(n_gen * n_env, 0, 0.5)
  dat <- data.frame(gen = gg, x = xx, y = y)

  # gen is a plain random factor (no relmat): gen:leg(x) is a factor x leg RR
  fit <- bbglr(y ~ leg(x), random = ~ gen + gen:leg(x), residual = ~ units,
               data = dat, nIter = 4000, burnIn = 1500, thin = 5, nChains = 2,
               verbose = FALSE)

  rn <- reaction_norm(fit, term = "gen:leg(x)", type = "random", plot = FALSE,
                      n_grid = 50)
  expect_true(all(c("id", "gradient", "value") %in% names(rn)))
  expect_equal(length(unique(rn$id)), n_gen)
  expect_equal(nrow(rn), n_gen * 50L)
  expect_equal(range(rn$gradient), c(-1, 1))

  # each genotype's fitted slope recovers the simulated reaction-norm slope
  sl <- vapply(split(rn, rn$id),
               function(d) stats::coef(stats::lm(value ~ gradient, d))[2], numeric(1))
  expect_gt(cor(sl[names(b)], b), 0.9)

  # original-scale axis matches the covariate range
  rn2 <- reaction_norm(fit, term = "gen:leg(x)", plot = FALSE, leg_basis = FALSE,
                       n_grid = 50)
  expect_equal(range(rn2$gradient), range(dat$x), tolerance = 1e-8)

  # ---- goodness-of-fit helpers on the same RR model --------------------------
  fv <- fitted(fit)
  rs <- residuals(fit)
  expect_equal(length(fv), nrow(fit$data))
  expect_false(is.null(attr(fv, "sd")))
  expect_equal(as.numeric(fv + rs), fit$data$y, tolerance = 1e-8)  # obs = fitted + resid

  mf <- model_fit(fit)
  expect_s3_class(mf, "breedRB_modelfit")
  expect_true(all(c("n", "r2", "rmse", "dic", "pD", "varE", "reliability") %in% names(mf)))
  expect_gt(mf$r2, 0.5)
  expect_true(all(c("gen", "gen:leg(x)") %in% mf$reliability$term))

  # DIC prefers the reaction norm over a slope-free model on data with real GxE
  fit0 <- bbglr(y ~ leg(x), random = ~ gen, residual = ~ units, data = dat,
                nIter = 3000, burnIn = 1000, nChains = 1, verbose = FALSE)
  expect_lt(model_fit(fit)$dic, model_fit(fit0)$dic)

  # ---- rr_gradient(): across-gradient analytics -----------------------------
  g <- rr_gradient(fit, term = "gen:leg(x)", n_grid = 10L, threshold = 0.2,
                   plot = FALSE)
  expect_s3_class(g, "breedRB_rrgradient")
  # correlation surface: symmetric, unit diagonal, in [-1, 1]
  expect_equal(dim(g$gcor), c(10L, 10L))
  expect_equal(unname(diag(g$gcor)), rep(1, 10L), tolerance = 1e-8)
  expect_true(all(g$gcor >= -1 - 1e-8 & g$gcor <= 1 + 1e-8))
  expect_equal(unname(g$gcor), unname(t(g$gcor)), tolerance = 1e-8)
  # covariance surface = Phi K Phi' -> its diagonal is the genetic variance,
  # so h2 = varg / (varg + varE) is a valid fraction in (0, 1)
  expect_true(all(g$h2$mean > 0 & g$h2$mean < 1))
  # K is the 2x2 (intercept, slope) coefficient (co)variance: symmetric, PSD,
  # positive variances
  expect_equal(dim(g$K), c(2L, 2L))
  expect_equal(g$K[1, 2], g$K[2, 1], tolerance = 1e-12)
  expect_true(all(diag(g$K) > 0))
  expect_gte(min(eigen(g$K, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
  # reliability / probability are per-genotype x gradient, bounded to [0, 1]
  expect_equal(dim(g$reliability), c(n_gen, 10L))
  expect_equal(dim(g$prob_top),   c(n_gen, 10L))
  expect_true(all(g$reliability >= 0 & g$reliability <= 1))
  expect_true(all(g$prob_top    >= 0 & g$prob_top    <= 1))
  # each gradient column: exactly k = round(0.2*n) genotypes expected in the top
  # (probabilities sum to k across genotypes at every point)
  expect_equal(unname(colSums(g$prob_top)), rep(round(0.2 * n_gen), 10L),
               tolerance = 0.05)
  expect_named(g$plots, c("cor", "h2", "reliability", "prob"))
  expect_true(all(vapply(g$plots, ggplot2::is.ggplot, logical(1))))
})

test_that("model_fit() breaks out reliability per Legendre degree for a q>1 RR", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  set.seed(5)
  n_gen <- 40; n_env <- 12
  gg  <- rep(paste0("g", 1:n_gen), each = n_env)
  xx  <- rep(seq(-1, 1, length.out = n_env), n_gen)
  L   <- legendre_basis(seq(-1, 1, length.out = n_env), 2, orthonormal = TRUE)
  a   <- rnorm(n_gen); b <- rnorm(n_gen, 0, 0.7); cc <- rnorm(n_gen, 0, 0.3)
  y   <- 10 + a[factor(gg)] +
    b[factor(gg)]  * rep(L[, 2], n_gen) +
    cc[factor(gg)] * rep(L[, 3], n_gen) + rnorm(n_gen * n_env, 0, 0.5)
  dat <- data.frame(gen = gg, x = xx, y = y)

  fit <- bbglr(y ~ leg(x, 2), random = ~ gen + gen:leg(x, 2), residual = ~ units,
               data = dat, nIter = 3000, burnIn = 1000, thin = 2, nChains = 2,
               verbose = FALSE)
  mf <- model_fit(fit)
  # the order-2 interaction is reported as two per-degree rows (intercept stays single)
  expect_true("gen" %in% mf$reliability$term)
  expect_true(all(c("gen:leg(x, 2):deg1", "gen:leg(x, 2):deg2") %in%
                    mf$reliability$term))
  expect_false("gen:leg(x, 2)" %in% mf$reliability$term)   # no collapsed-across-degree row
  # each degree row counts exactly the genotypes (not genotypes x degrees)
  degrows <- mf$reliability[grepl(":deg[0-9]+$", mf$reliability$term), ]
  expect_true(all(degrows$n_effects == n_gen))

  # varcomp() exposes the per-coefficient variances (diagonal of K) as var() rows
  # -- not visible from the single shared BGLR interaction component -- plus the
  # covariance off-diagonal, and they equal the K reported by rr_gradient()
  vc <- varcomp(fit)
  expect_true(all(c("var(gen)", "var(gen:leg(x, 2):deg1)",
                    "var(gen:leg(x, 2):deg2)") %in% vc$term))
  expect_true("cov(gen:leg(x, 2):deg1, gen:leg(x, 2):deg2)" %in% vc$term)
  g <- rr_gradient(fit, term = "gen:leg(x, 2)", n_grid = 5L, plot = FALSE)
  vrows <- vc$mean[grepl("^var\\(", vc$term)]
  expect_equal(vrows, diag(g$K), tolerance = 1e-6)

  # term resolution is whitespace-insensitive: the deparsed label carries a
  # space ("gen:leg(x, 2)") but users routinely type it without ("gen:leg(x,2)")
  g2 <- rr_gradient(fit, term = "gen:leg(x,2)", n_grid = 5L, plot = FALSE)
  expect_equal(g2$K, g$K, tolerance = 1e-12)
  rn_ns <- reaction_norm(fit, term = "gen:leg(x,2)", plot = FALSE, n_grid = 10L)
  rn_sp <- reaction_norm(fit, term = "gen:leg(x, 2)", plot = FALSE, n_grid = 10L)
  expect_equal(rn_ns$value, rn_sp$value, tolerance = 1e-12)

  # gxe(): per-genotype adaptability / responsiveness / stability table
  gx <- gxe(fit, term = "gen:leg(x, 2)")
  expect_s3_class(gx, "breedRB_gxe")
  expect_true(all(c("id", "adaptability", "responsiveness", "cv_ge",
                    "adaptability_rank", "stability_rank") %in% names(gx)))
  expect_equal(nrow(gx), n_gen)
  expect_true(all(gx$cv_ge >= 0))
  # ordered by adaptability (best first, higher = better by default)
  expect_equal(gx$adaptability, sort(gx$adaptability, decreasing = TRUE))
  # stability rank 1 == smallest coefficient of variation
  expect_equal(gx$id[gx$stability_rank == 1L], gx$id[which.min(gx$cv_ge)])
  # adaptability equals the mean of each genotype's phenotype-scale curve
  rn_ph <- reaction_norm(fit, term = "gen:leg(x, 2)", plot = FALSE, n_grid = 100L)
  amean <- vapply(split(rn_ph$value, rn_ph$id), mean, numeric(1))
  expect_equal(unname(gx$adaptability[match(names(amean), gx$id)]),
               unname(amean), tolerance = 1e-8)

  # heritability() of a single RR interaction -> per-coefficient h2 rows
  h2 <- heritability(fit, genetic = "gen:leg(x, 2)")
  expect_s3_class(h2, "breedRB_h2")
  expect_equal(nrow(h2$summary), 3L)                       # intercept + 2 degrees
  expect_true(all(c("h2(gen)", "h2(gen:leg(x, 2):deg1)",
                    "h2(gen:leg(x, 2):deg2)") %in% h2$summary$quantity))
  expect_true(all(h2$summary$mean > 0 & h2$summary$mean < 1))
  # whitespace-insensitive too
  expect_equal(heritability(fit, genetic = "gen:leg(x,2)")$summary$mean,
               h2$summary$mean, tolerance = 1e-12)
  # per-coefficient h2 == K_jj / (K_jj + varE) from rr_gradient()/varcomp()
  vE <- vc$mean[vc$term == "varE"]
  expect_equal(unname(colMeans(h2$draws)),
               unname(diag(g$K) / (diag(g$K) + vE)), tolerance = 0.02)
})
