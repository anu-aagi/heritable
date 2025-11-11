test_that("Helper functions working with lme4", {
  skip_on_cran()
  skip_if_not_installed("lme4")

  dat <- agridat::john.alpha

  model_ran_lmer <- lme4::lmer(data = dat, formula = yield ~ rep + (1 | gen) + (1 | rep:block))
  model_fixed_lmer <- lme4::lmer(data = dat, formula = yield ~ rep + gen + (1 | rep:block))

  # Pull out the terms correctly?
  model_ran_lmer <- pull_terms.lmerMod(model_ran_lmer)

  expect_equal(model_ran_asreml_terms$fixed, model_ran_lmer$fixed)
  expect_equal(model_ran_asreml_terms$random, model_ran_lmer$random)

  # Fit counterpart models?
  model_lmer_counter <- fit_counterpart_model.lmerMod(model_ran_lmer, target = "gen")
  model_lmer_fixed_counter <- fit_counterpart_model.lmerMod(model_fixed_lmer, target = "gen")

  expect_true(any(grepl(target, pull_terms.lmerMod(model_lmer_counter)$fixed)))
  expect_true(any(grepl(target, pull_terms.lmerMod(model_lmer_fixed_counter)$random)))
})


test_that("Helper functions working with asreml", {
  skip_on_cran()
  skip_if_not_installed("asreml")

  dat <- agridat::john.alpha

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
  model_ran_asreml_terms <- pull_terms.asreml(model_ran_asreml)

  expect_equal(model_ran_asreml_terms$fixed, model_ran_lmer$fixed)
  expect_equal(model_ran_asreml_terms$random, model_ran_lmer$random)

  # can it fit counterpart models?
  model_asreml_counter <- fit_counterpart_model.asreml(model_ran_asreml, target = "gen")
  model_asreml_fixed_counter <- fit_counterpart_model.asreml(model_fixed_asreml, target = "gen")

  expect_true(any(grepl(target, pull_terms.asreml(model_asreml_counter)$fixed)))
  expect_true(any(grepl(target, pull_terms.asreml(model_asreml_fixed_counter)$random)))
})
