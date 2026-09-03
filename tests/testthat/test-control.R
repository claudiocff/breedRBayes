# Unit tests for bbglr_control() and the explained-variance rank helper.

test_that("bbglr_control() validates and defaults sensibly", {
  ctrl <- bbglr_control()
  expect_s3_class(ctrl, "breedRB_control")
  expect_equal(ctrl$nIter, 5000L)
  expect_equal(ctrl$exp_var_rank, 0.99)

  # NULL exp_var_rank means "keep all components" (stored as NA)
  expect_true(is.na(bbglr_control(exp_var_rank = NULL)$exp_var_rank))
  expect_true(is.na(bbglr_control(exp_var_rank = NA)$exp_var_rank))

  expect_error(bbglr_control(exp_var_rank = 0),   "exp_var_rank")
  expect_error(bbglr_control(exp_var_rank = 1.5), "exp_var_rank")
  expect_error(bbglr_control(burnIn = 6000, nIter = 5000), "burnIn")
  expect_error(bbglr_control(nIter = 2.5), "nIter")
  expect_error(bbglr_control(thin = 0), "thin")
})

test_that(".svds_expvar reaches the target without the full spectrum", {
  skip_if_not_installed("RSpectra")
  set.seed(11)
  # strong low-rank structure: 5 latent factors + tiny noise
  n <- 200L; p <- 400L; r0 <- 5L
  Mc <- scale(matrix(rnorm(n * r0), n, r0) %*% matrix(rnorm(r0 * p), r0, p) +
                0.01 * matrix(rnorm(n * p), n, p),
              center = TRUE, scale = FALSE)

  # ground truth from a full eigendecomposition
  lam_full <- pmax(eigen(tcrossprod(Mc), symmetric = TRUE)$values, 0)
  k_ref    <- .rank_for_expvar(lam_full, 0.99)

  hit <- .svds_expvar(Mc, 0.99, total = sum(rowSums(Mc^2)))
  expect_false(is.null(hit))                      # structure -> succeeds incrementally
  expect_equal(ncol(hit$vectors), k_ref)          # same rank as the full-spectrum route
  expect_lt(ncol(hit$vectors), n)                 # far fewer than n components
  # the returned eigenvalues match the leading full-spectrum eigenvalues
  expect_equal(hit$evals, lam_full[seq_len(k_ref)], tolerance = 1e-4)
})

test_that(".svds_expvar gives up (NULL) when the target needs most of the spectrum", {
  skip_if_not_installed("RSpectra")
  set.seed(12)
  # i.i.d. markers: no low-rank structure, 99% needs almost all components
  Mc <- scale(matrix(rnorm(60 * 400), 60, 400), center = TRUE, scale = FALSE)
  expect_null(.svds_expvar(Mc, 0.99, total = sum(rowSums(Mc^2))))
})

test_that(".rank_for_expvar picks the smallest leading set", {
  ev <- c(90, 6, 3, 1)                    # descending eigenvalues, total 100
  expect_equal(.rank_for_expvar(ev, 0.90), 1L)   # 90% -> first alone
  expect_equal(.rank_for_expvar(ev, 0.95), 2L)   # 96% at k=2
  expect_equal(.rank_for_expvar(ev, 0.99), 3L)   # 99% at k=3
  expect_equal(.rank_for_expvar(ev, 1.0),  4L)   # everything
  expect_equal(.rank_for_expvar(c(0, 0, 0), 0.9), 1L)  # degenerate -> at least 1
})
