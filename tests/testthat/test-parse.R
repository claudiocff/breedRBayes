test_that("parse_model splits fixed/random/residual and specials", {
  p <- parse_model(yield ~ 1 + env,
                   ~ vm(gen, G) + leg(idx, 2):vm(gen, G),
                   residual = ~ dsum(~units | env))
  expect_false(p$response$multitrait)
  expect_identical(p$response$traits, "yield")
  expect_true(p$intercept)
  expect_length(p$random, 2)
  expect_true(p$residual$hetero)
  expect_identical(p$residual$group, "env")

  # a vm component is recognised in one of the random terms
  kinds <- unlist(lapply(p$random, function(t) vapply(t$components, `[[`, character(1), "type")))
  expect_true("vm" %in% kinds)
  expect_true("leg" %in% kinds)
})

test_that("cbind LHS triggers multi-trait", {
  p <- parse_model(cbind(y1, y2) ~ env, ~ vm(gen, G))
  expect_true(p$response$multitrait)
  expect_identical(p$response$traits, c("y1", "y2"))
})

test_that("unknown special errors", {
  expect_error(parse_model(y ~ 1, ~ wibble(gen)), "Unknown special")
})
