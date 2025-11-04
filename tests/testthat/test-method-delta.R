test_that("delta method works", {
  skip_if_not_installed("asreml")
  skip_on_cran()

  # Genotype as random effect
  model <- asreml::asreml(
    fixed = yield ~ rep,
    random = ~ gen + rep:block,
    data = agridat::john.alpha,
    trace = FALSE
  )

  target <- "gen"
})
