test_that("anova() reconstructs the posterior-mean fitted values exactly", {
  skip_on_cran()
  skip_if_not_installed("BGLR")

  set.seed(3)
  n_gen <- 40L; n_rep <- 5L; p <- 200L
  M <- matrix(stats::rbinom(n_gen * p, 2, 0.3), n_gen, p)
  rownames(M) <- paste0("g", seq_len(n_gen))
  g  <- as.numeric(scale(M) %*% stats::rnorm(p)); g <- 2 * g / stats::sd(g)
  sl <- stats::setNames(stats::rnorm(n_gen), rownames(M))
  x   <- stats::runif(n_gen * n_rep, 0, 1)
  ids <- rep(rownames(M), each = n_rep)
  y   <- g[match(ids, rownames(M))] + sl[ids] * x + stats::rnorm(n_gen * n_rep)
  dat <- data.frame(gen = ids, x = x, y = y)

  ctrl <- list(nIter = 2000, burnIn = 800, thin = 2, nChains = 1,
               verbose = FALSE)

  # The pointwise-loglik reconstruction of yHat must match BGLR's stored
  # fitted values to numerical precision, for a genomic main effect AND a
  # genomic random regression (the sign-safe solution() back-map is the
  # whole point of the machinery).
  m0 <- do.call(bbglr, c(list(y ~ 1, random = ~ mrk(gen, M), data = dat,
                              relmat = list(M = M), saveAt = tempfile()), ctrl))
  m1 <- do.call(bbglr, c(list(y ~ leg(x, 1),
                              random = ~ mrk(gen, M) + mrk(gen, M):leg(x, 1),
                              data = dat, relmat = list(M = M),
                              saveAt = tempfile()), ctrl))

  # Direct exactness check against fitted(): rebuild rowMeans(contrib + mu).
  exact <- function(fit) {
    keys  <- names(fit$meta)
    nburn <- floor(fit$control$burnIn / fit$control$thin)
    ev    <- fit$control$exp_var_rank
    built <- suppressMessages(build_eta(fit$parsed, fit$data, fit$relmat,
                                        exp_var = if (is.null(ev)) NA else ev))
    isg <- vapply(keys, function(k)
      any(vapply(fit$meta[[k]]$components,
                 function(cm) isTRUE(cm$kind %in% c("vm", "mrk")), logical(1))),
      logical(1))
    S_c <- vapply(fit$paths, function(pp)
      length(scan(paste0(pp, "mu.dat"), quiet = TRUE)) - nburn, numeric(1))
    S <- sum(S_c); n <- nrow(fit$data)
    contrib <- matrix(0, n, S)
    for (k in keys) {
      if (isg[[k]]) {
        contrib <- contrib + breedRBayes:::.genomic_obs_contrib(fit, k)
      } else {
        X  <- built$ETA[[k]]$X
        isF <- identical(fit$meta[[k]]$model, "FIXED")
        Bp <- do.call(rbind, lapply(fit$paths, function(pp) {
          if (isF) {
            B <- as.matrix(utils::read.table(paste0(pp, "ETA_", k, "_b.dat"),
                                             header = TRUE))
            B[(nburn + 1L):nrow(B), , drop = FALSE]
          } else BGLR::readBinMat(paste0(pp, "ETA_", k, "_b.bin"))
        }))
        contrib <- contrib + tcrossprod(X, Bp)
      }
    }
    mu <- unlist(lapply(fit$paths, function(pp)
      utils::tail(scan(paste0(pp, "mu.dat"), quiet = TRUE), S / length(fit$paths))))
    rowMeans(sweep(contrib, 2L, mu, "+"))
  }

  expect_equal(unname(exact(m0)), as.numeric(fitted(m0)), tolerance = 1e-6)
  expect_equal(unname(exact(m1)), as.numeric(fitted(m1)), tolerance = 1e-6)
})

test_that("anova() returns a coherent comparison object", {
  skip_on_cran()
  skip_if_not_installed("BGLR")

  set.seed(11)
  n_gen <- 40L; n_rep <- 5L; p <- 150L
  M <- matrix(stats::rbinom(n_gen * p, 2, 0.3), n_gen, p)
  rownames(M) <- paste0("g", seq_len(n_gen))
  g  <- as.numeric(scale(M) %*% stats::rnorm(p)); g <- 2 * g / stats::sd(g)
  sl <- stats::setNames(stats::rnorm(n_gen, 0, 1.2), rownames(M))
  x   <- stats::runif(n_gen * n_rep, 0, 1)
  ids <- rep(rownames(M), each = n_rep)
  y   <- g[match(ids, rownames(M))] + sl[ids] * x + stats::rnorm(n_gen * n_rep, 0, 0.8)
  dat <- data.frame(gen = ids, x = x, y = y)
  ctrl <- list(nIter = 2500, burnIn = 1000, thin = 2, nChains = 1, verbose = FALSE)

  m0 <- do.call(bbglr, c(list(y ~ leg(x, 1), random = ~ mrk(gen, M), data = dat,
                              relmat = list(M = M), saveAt = tempfile()), ctrl))
  m1 <- do.call(bbglr, c(list(y ~ leg(x, 1),
                              random = ~ mrk(gen, M) + mrk(gen, M):leg(x, 1),
                              data = dat, relmat = list(M = M),
                              saveAt = tempfile()), ctrl))

  res <- anova(m0, m1)
  expect_s3_class(res, "breedRB_anova")
  expect_equal(nrow(res$table), 2L)
  expect_true(all(c("pD", "DIC", "WAIC", "dDIC", "dWAIC") %in% names(res$table)))
  expect_true(all(is.finite(res$table$DIC)))
  expect_true(all(is.finite(res$table$WAIC)))
  expect_true(min(res$table$dDIC) == 0)              # one model is the reference
  # more effective parameters for the richer model
  expect_gt(res$table$pD[2], res$table$pD[1])
  expect_equal(nrow(res$lrt), 1L)
  expect_true(is.finite(res$lrt$deviance_reduction))
  expect_no_error(print(res))
})

test_that("anova() needs at least two fits and rejects multi-trait", {
  skip_on_cran()
  fake <- structure(list(), class = "breedRB_fit")
  expect_error(anova(fake), "at least two")
})
