test_that("heritability for asreml works", {
  skip_if_not_installed("asreml")
  skip_on_cran()


# random gen --------------------------------------------------------------
  model <- asreml::asreml(yield ~ rep,
    random = ~ gen + rep:block,
    data = agridat::john.alpha,
    trace = FALSE
  )
  H2(model, target = "gen")
  expect_equal(H2(model, target = "gen"),
               c("Cullis" = 0.8090841,
                 "Oakey" =  0.8090728,
                 "Delta" = 0.8397128, # CHECK VALUE
                 "Piepho" = 0.8029759,
                 "Naive" = 0.6364751), # CHECK VALUE
               tolerance = 1e-7)


# fixed gen ---------------------------------------------------------------
  model <- asreml::asreml(yield ~ rep + gen,
                          random = ~ rep:block,
                          data = agridat::john.alpha,
                          trace = FALSE
  )
  expect_equal(H2(model, target = "gen"),
               c("Cullis" = NA,
                 "Oakey" =  NA,
                 "Delta" = 0.8030227, # CHECK VALUE
                 "Piepho" = 0.8029759,
                 "Naive" = NA),
               tolerance = 1e-7)

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
