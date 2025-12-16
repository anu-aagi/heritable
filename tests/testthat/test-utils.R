test_that("Helper functions working with asreml", {
  skip_if_not_installed("asreml")
  skip_on_ci()
  skip_on_cran()

  asreml_model_random <- readRDS(test_path("fixtures/asreml_model_random.rds"))
  asreml_model_fixed <- fit_counterpart_model(asreml_model_random, target = "gen")

  # Can we pull out the correct terms?
  # Can we fit the counter model?
  expect_named(pull_terms(asreml_model_random), c("fixed", "random"))
  expect_true("gen" %in% pull_terms(asreml_model_random)$random)
  expect_true("gen" %in% pull_terms(asreml_model_fixed)$fixed)
})


test_that("Helper functions working with lme4", {
  skip_on_cran()
  skip_if_not_installed("lme4")

  lme4_lettuce <- readRDS(here::here("vignettes/fixtures/lettuce_lme4.rds"))
  lmer_model_random <- readRDS(test_path("fixtures/lmer_model_random.rds"))
  lmer_model_fixed <- fit_counterpart_model(lmer_model_random, target = "gen")

  # Can we pull out the correct terms?
  # Can we fit the counter model
  expect_named(pull_terms(lmer_model_random), c("fixed", "random"))
  expect_true("gen" %in% pull_terms(lmer_model_random)$random)
  expect_true("gen" %in% pull_terms(lmer_model_fixed)$fixed)

  lm_model <- fit_counterpart_model(lme4_lettuce)
})



