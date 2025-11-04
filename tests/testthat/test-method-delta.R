test_that("delta method works", {
  # Genotype as random effect
  model <- asreml::asreml(
    fixed = yield ~ rep,
    random = ~ gen + rep:block,
    data = agridat::john.alpha,
    trace = FALSE
  )

  target <- "gen"
})
