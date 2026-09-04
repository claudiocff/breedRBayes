test_that("MAD rule flags a clear outlier in a numeric vector", {
  v <- c(1, 2, 3, 4, 5, 100)

  fl <- outlier_rm(v, action = "flag", verbose = FALSE)
  expect_type(fl, "logical")
  expect_equal(which(fl), 6L)

  na <- outlier_rm(v, action = "na", verbose = FALSE)
  expect_true(is.na(na[6]))
  expect_equal(na[1:5], v[1:5])

  rm <- outlier_rm(v, action = "remove", verbose = FALSE)
  expect_equal(as.numeric(rm), v[1:5])

  info <- attr(fl, "outliers")
  expect_equal(info$index, 6L)
  expect_equal(info$k, 3)
})

test_that("grouping screens strata independently (formula RHS and by=)", {
  set.seed(1)
  df <- data.frame(env = rep(c("A", "B"), each = 20),
                   yield = c(rnorm(20, 10), rnorm(20, 50)))
  df$yield[5]  <- 1000                    # extreme within env A
  df$yield[35] <- -900                    # extreme within env B

  cl <- outlier_rm(yield ~ env, data = df, k = 3, verbose = FALSE)
  fl <- attr(cl, "outliers")$flag
  expect_true(all(c(5L, 35L) %in% which(fl)))
  expect_equal(sum(is.na(cl$yield)), sum(fl))
  expect_true("yield_outlier" %in% names(cl))
  expect_equal(cl$yield_outlier, fl)

  # a global (ungrouped) screen would miss env-relative outliers because the
  # two environments differ so much in scale; grouping is what makes it work
  cl_by <- outlier_rm(yield ~ 1, data = df, by = "env", k = 3, verbose = FALSE)
  expect_equal(attr(cl_by, "outliers")$flag, fl)
})

test_that("remove drops rows; the outliers attribute still indexes original rows", {
  df <- data.frame(g = "x", y = c(1, 2, 3, 4, 5, 100))
  out <- outlier_rm(y ~ g, data = df, k = 3, action = "remove", verbose = FALSE)
  expect_equal(nrow(out), 5L)
  expect_equal(attr(out, "outliers")$index, 6L)   # attribute describes the original data
})

test_that("tiny strata and constant vectors are left untouched", {
  # fewer than min_n observations -> not scored, nothing flagged
  df <- data.frame(g = c("a", "a", "b", "b"), y = c(1, 2, 3, 999))
  out <- outlier_rm(y ~ g, data = df, k = 3, action = "flag", min_n = 5, verbose = FALSE)
  expect_false(any(attr(out, "outliers")$flag))

  # constant vector: MAD and MeanAD both zero -> no division by zero, no flags
  expect_false(any(outlier_rm(rep(7, 10), action = "flag", verbose = FALSE)))
})

test_that("NA responses never get flagged and input validation fires", {
  v <- c(1, 2, 3, NA, 5, 100)
  fl <- outlier_rm(v, action = "flag", verbose = FALSE)
  expect_false(fl[4])
  expect_true(fl[6])

  expect_error(outlier_rm(1:5, k = -1), "positive")
  expect_error(outlier_rm(letters), "numeric")
  expect_error(outlier_rm(y ~ g), "supply `data`")
})

test_that("residual-based screen on a fitted model returns a cleaned data frame", {
  skip_on_cran()
  skip_if_not_installed("BGLR")

  set.seed(42)
  n_gen <- 40; n_rep <- 5
  M <- matrix(rbinom(n_gen * 200, 2, 0.3), n_gen, 200)
  rownames(M) <- paste0("g", seq_len(n_gen))
  u   <- as.numeric(scale(M) %*% rnorm(200)) ; u <- 3 * u / sd(u)
  ids <- rep(rownames(M), each = n_rep)
  dat <- data.frame(gen = ids, y = u[match(ids, rownames(M))] + rnorm(length(ids)))
  dat$y[7]  <- dat$y[7]  + 40      # gross positive residual outlier
  dat$y[100] <- dat$y[100] - 40    # gross negative residual outlier

  td <- tempfile(); dir.create(td)
  fit <- bbglr(y ~ 1, random = ~ mrk(gen, M), data = dat,
               relmat = list(M = M), nIter = 2000, burnIn = 800, thin = 2,
               nChains = 1, verbose = FALSE, saveAt = td)

  clean <- outlier_rm(fit, k = 4, verbose = FALSE)     # default action = "remove"
  expect_s3_class(clean, "data.frame")
  expect_equal(names(clean), names(dat))               # no marker column on remove
  expect_lt(nrow(clean), nrow(dat))                    # some rows dropped
  info <- attr(clean, "outliers")
  expect_true(all(c(7L, 100L) %in% info$index))        # the injected outliers caught
  expect_equal(nrow(clean), nrow(dat) - length(info$index))

  # action = "na" keeps every row, blanks the response, adds a marker column
  na <- outlier_rm(fit, k = 4, action = "na", verbose = FALSE)
  expect_equal(nrow(na), nrow(dat))
  expect_true("y_outlier" %in% names(na))
  expect_true(is.na(na$y[7]) && is.na(na$y[100]))
  expect_true(all(na$y_outlier[c(7, 100)]))

  # action = "flag" changes nothing but marks the rows
  fl <- outlier_rm(fit, k = 4, action = "flag", verbose = FALSE)
  expect_equal(fl$y, dat$y)
  expect_true(all(attr(fl, "outliers")$flag[c(7, 100)]))
})
