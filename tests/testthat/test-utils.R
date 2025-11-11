test_that("Helper functions working with lme4", {
  skip_on_cran()
  skip_if_not_installed("lme4")

  dat <- agridat::john.alpha
  target <- "gen"

  model_ran_lmer <- lme4::lmer(data = dat, formula = yield ~ rep + (1 | gen) + (1 | rep:block))
  model_fixed_lmer <- lme4::lmer(data = dat, formula = yield ~ rep + gen + (1 | rep:block))

  # Pull out the terms correctly?
  model_ran_lmer_terms <- pull_terms(model_ran_lmer)

  expect_named(model_ran_lmer_terms, c("fixed", "random"))

  # Fit counterpart models?
  model_lmer_counter <- fit_counterpart_model.lmerMod(model_ran_lmer, target = "gen")
  model_lmer_fixed_counter <- fit_counterpart_model.lmerMod(model_fixed_lmer, target = "gen")

  expect_true(any(grepl(target, pull_terms(model_lmer_counter)$fixed)))
  expect_true(any(grepl(target, pull_terms(model_lmer_fixed_counter)$random)))
})


test_that("Helper functions working with asreml", {
  skip_on_cran()
  skip_if_not_installed("asreml")

  dat <- agridat::john.alpha
  target <- "gen"

  # random genotype effect
  model_ran_asreml <- asreml::asreml(
    data = dat,
    fixed = yield ~ rep,
    random = ~ gen + rep:block,
    trace = FALSE
  )

  model_fixed_asreml <- asreml::asreml(
    data = dat,
    fixed = yield ~ rep + gen,
    random = ~ rep:block,
    trace = FALSE
  )

  # can it pull out the correct terms?
  model_ran_asreml_terms <- pull_terms(model_ran_asreml)

  expect_named(model_ran_asreml_terms, c("fixed", "random"))

  # can it fit counterpart models?
  model_asreml_counter <- fit_counterpart_model.asreml(model_ran_asreml, target = "gen")
  model_asreml_fixed_counter <- fit_counterpart_model.asreml(model_fixed_asreml, target = "gen")

  expect_true(any(grepl(target, pull_terms(model_asreml_counter)$fixed)))
  expect_true(any(grepl(target, pull_terms(model_asreml_fixed_counter)$random)))
})
