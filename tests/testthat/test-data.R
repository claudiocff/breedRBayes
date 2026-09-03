test_that("bundled wheat example data is well formed", {
  data(wheat, package = "breedRBayes")
  data(wheat_M, package = "breedRBayes")

  expect_s3_class(wheat, "data.frame")
  expect_equal(names(wheat), c("gen", "env", "yield"))
  expect_equal(nrow(wheat), 599L * 4L)
  expect_equal(nlevels(wheat$gen), 599L)
  expect_equal(nlevels(wheat$env), 4L)
  expect_true(is.numeric(wheat$yield) && all(is.finite(wheat$yield)))

  expect_true(is.matrix(wheat_M))
  expect_equal(dim(wheat_M), c(599L, 1279L))
  expect_true(all(wheat_M %in% c(0, 1)))
  # markers keyed to the genotype levels so mrk(gen, wheat_M) lines up
  expect_equal(rownames(wheat_M), levels(wheat$gen))
})
