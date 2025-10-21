test_that("heritability for asreml works", {
  skip_if_not_installed("asreml")
  skip_on_cran()

  fit <- asreml::asreml(yield ~ rep,
    random = ~ gen + rep:block,
    data = agridat::john.alpha,
    trace = FALSE
  )

  expect_equal(H2(fit, target = "gen"), 0.8090841)
})
