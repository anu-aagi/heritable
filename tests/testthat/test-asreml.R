test_that("heritability for asreml works", {
  skip_if_not_installed("asreml")
  skip_on_cran()

  fit <- asreml::asreml(yield ~ rep,
    random = ~ gen + rep:block,
    data = agridat::john.alpha,
    trace = FALSE
  )

  expect_equal(H2(fit, target = "gen"), 0.8090841)
  expect_equal(H2(fit, target = "gen", method = "Piepho")), 0.8090841)
}

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
})