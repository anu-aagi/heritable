test_that("Inner checks are triggered", {
  skip_if_not_installed("asreml")
  skip_on_cran()

  lme4_lettuce <- readRDS(test_path("fixtures/lettuce_lme4.rds"))
  lmer_model_random <- readRDS(test_path("fixtures/lmer_model_random.rds"))
  asreml_model_random <- readRDS(test_path("fixtures/asreml_model_random.rds"))
  asreml_model_fixed <- readRDS(test_path("fixtures/asreml_model_fixed.rds"))
  target <- "gen"

  # Won't converge
  lmer_model_failed_converge <- lmer_model_random
  lmer_model_failed_converge@optinfo$conv$opt <- 10
  asreml_model_failed_converge <- asreml_model_random
  asreml_model_failed_converge$converge <- FALSE

  expect_error(H2(model = c(asreml_model_random, asreml_model_fixed), target = target))
  expect_warning(H2(lmer_model_failed_converge, target = target, method = "Oakey"))
  expect_warning(H2(asreml_model_failed_converge, target = target, method = "Oakey"))
  expect_error(H2(model = asreml_model_random, target = "foo"))
  expect_error(H2(model = asreml_model_fixed, target = "gen"))

  # Target level
  expect_false(check_target_random(asreml_model_fixed, target))
  expect_true(check_target_random(asreml_model_random, target))
  expect_error(H2(model = asreml_model_random, target = c("foo", "rep")))

  # Method level
  expect_error(H2(asreml_model_fixed, target, "Oakey"))
  expect_error(H2(asreml_model_fixed, target, "Cullis"))

  # Model level
  expect_true(check_single_random_effect(pull_terms(lme4_lettuce)))
  expect_false(check_single_random_effect(pull_terms(lmer_model_random)))
})
