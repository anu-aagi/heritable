test_that("Inner checks are triggered", {
  skip_if_not_installed("asreml")
  skip_on_cran()

  asreml_model_random <- readRDS(test_path("fixtures/asreml_model_random.rds"))
  asreml_model_fixed <- readRDS(test_path("fixtures/asreml_model_fixed.rds"))
  target <- "gen"

  # Won't converge
  asreml_model_failed_converge <- asreml_model_random
  asreml_model_failed_converge$converge <- FALSE

  expect_error(H2(model = c(asreml_model_random, asreml_model_fixed), target = target))
  expect_warning(H2(asreml_model_failed_converge, target = target, method = "Oakey"))
  expect_error(H2(model = asreml_model_random, target = "foo"))

  # Target level
  expect_false(check_target_random(asreml_model_fixed, target))
  expect_true(check_target_random(asreml_model_random, target))

  # Method level
  expect_error(H2(asreml_model_fixed, target, "Oakey"))
  expect_error(H2(asreml_model_fixed, target, "Cullis"))
})
