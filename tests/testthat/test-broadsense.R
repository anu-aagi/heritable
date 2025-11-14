test_that("heritability for asreml works", {
  skip_if_not_installed("asreml")
  skip_on_cran()


# random gen --------------------------------------------------------------
  truth <-
    c(
      "Cullis" = 0.8090841,
      "Oakey" = 0.8090728,
      "Piepho" = 0.8029759,
      "Delta" = 0.8397128, # CHECK VALUE
      "Naive" = 0.6364751 # CHECK VALUE
    )

  model <- asreml::asreml(yield ~ rep,
                          random = ~ gen + rep:block,
                          data = agridat::john.alpha,
                          trace = FALSE)
  expect_equal(H2(model, target = "gen"), truth, tolerance = 1e-7)

  model <- lme4::lmer(data = dat, formula = yield ~ rep + (1 | gen) + (1 | rep:block))
  expect_equal(H2(model, target = "gen"), truth, tolerance = 1e-7)

  # fixed gen ---------------------------------------------------------------
  truth <- c("Cullis" = NA,
             "Oakey" =  NA,
             "Piepho" = 0.8029759,
             "Delta" = 0.8030227, # CHECK VALUE
             "Naive" = NA)

  model <- asreml::asreml(yield ~ rep + gen,
                          random = ~ rep:block,
                          data = agridat::john.alpha,
                          trace = FALSE
  )
  expect_equal(H2(model, target = "gen"), truth, tolerance = 1e-7)

  model <- lme4::lmer(data = dat, formula = yield ~ rep + gen + (1 | rep:block))
  expect_equal(H2(model, target = "gen"), truth, tolerance = 1e-7)

})

test_that("counterpart model can be fitted", {
  skip_if_not_installed("asreml")
  skip_on_cran()

  fit <- asreml::asreml(yield ~ rep,
    random = ~ gen + rep:block,
    data = agridat::john.alpha,
    trace = FALSE
  )

  fit_counter <- fit_counterpart_model.asreml(fit, target = "gen")

  expect_true("gen" %in% pull_terms(fit_counter)$fixed)
  expect_false("gen" %in% pull_terms(fit_counter)$random)

   fit_fixed <- asreml::asreml(yield ~ rep + gen,
    random = ~ rep:block,
    data = agridat::john.alpha,
    trace = FALSE
  )

  fit_counter <- fit_counterpart_model.asreml(fit_fixed, target = "gen")

  expect_false("gen" %in% pull_terms(fit_counter)$fixed)
  expect_true("gen" %in% pull_terms(fit_counter)$random)

})
