test_that("delta method works", {
  skip_if_not_installed("asreml")
  skip_on_cran()

  # Genotype as random effect
  model_random <- asreml::asreml(
    fixed = yield ~ rep,
    random = ~ gen + rep:block,
    data = agridat::john.alpha,
    trace = FALSE
  )

  target <- "gen"

  H2_Delta.asreml(model_random, target)
})
