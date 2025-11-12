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

  # Genotype as fixed effect
  model_fixed <- asreml::asreml(
    fixed = yield ~ rep + gen,
    random = ~ rep:block,
    data = agridat::john.alpha,
    trace = FALSE
  )

  # Genotype as both fixed and random effect
  model_both <- asreml::asreml(
    fixed = yield ~ gen,
    random = ~ gen:block,
    data = agridat::john.alpha,
    trace = FALSE 
  )

  target <- "gen"
  
  expect_equal(H2_Delta.asreml(model_random, target), 0.8090841)
  expect_equal(H2_Delta.asreml(model_fixed, target, mean = "harmonic"),  0.8029759) 
  expect_lt(H2_Delta.asreml(model_random, target, mean = "harmonic"), H2_Delta.asreml(model_random, target)) 
  expect_error(H2_Delta.asreml(model_both, target, mean = "harmonic"))
})
