test_that("model level checks", {
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

  # Won't converge
  model_failed_converge <- asreml::asreml(
    fixed = yield ~ rep + gen,
    random = ~ rep:block, # Add complex interaction
    maxiter = 2, # Limit iterations
    data = agridat::john.alpha,
    trace = FALSE
  )

  target <- "gen"

  expect_error(H2(model = c(model_random, model_fixed), target = target))
  H2(model_failed_converge, target = target)
})
