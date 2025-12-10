test_that("Delta method works for ASREML", {
  skip_if_not_installed("asreml")
  skip_on_cran()

  asreml_model_random <- readRDS(test_path("fixtures/asreml_model_random.rds"))
  asreml_model_fixed <- readRDS(test_path("fixtures/asreml_model_fixed.rds"))
  asreml_model_both <- readRDS(test_path("fixtures/asreml_model_both.rds"))

  target <- "gen"

  expect_equal(H2_Delta(asreml_model_random, target), 0.8090841)
  expect_error(H2_Delta(asreml_model_fixed, target, aggregate = "harmonic"))
  expect_lt(H2_Delta(asreml_model_random, target, type = "BLUP", aggregate = "harmonic"), H2_Delta(asreml_model_random, target, type = "BLUP"))
  expect_error(h2_Delta(asreml_model_both, target, aggregate = "harmonic"))

  res_ls <- H2_Delta_by_genotype(asreml_model_random, target = "gen", type = "BLUP")
  expect_named(H2_Delta_by_genotype(asreml_model_random, target, type = "BLUP"), levels(asreml_model_random$mf$gen))
  expect_true(length(res_ls) == length(levels(asreml_model_random$mf[["gen"]])))

  H2_mat <- H2_Delta_pairwise(asreml_model_random, target = "gen", type = "BLUP")
  expect_true(is.matrix(as.matrix(H2_mat)))
  expected <- rowMeans(as.matrix(H2_mat), na.rm = TRUE)

  # compare numeric values
  expect_equal(as.numeric(unlist(res_ls)), as.numeric(expected), tolerance = 1e-7)
})
