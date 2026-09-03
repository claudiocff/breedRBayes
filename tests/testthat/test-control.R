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

test_that(".rank_for_expvar picks the smallest leading set", {
  ev <- c(90, 6, 3, 1)                    # descending eigenvalues, total 100
  expect_equal(.rank_for_expvar(ev, 0.90), 1L)   # 90% -> first alone
  expect_equal(.rank_for_expvar(ev, 0.95), 2L)   # 96% at k=2
  expect_equal(.rank_for_expvar(ev, 0.99), 3L)   # 99% at k=3
  expect_equal(.rank_for_expvar(ev, 1.0),  4L)   # everything
  expect_equal(.rank_for_expvar(c(0, 0, 0), 0.9), 1L)  # degenerate -> at least 1
})
