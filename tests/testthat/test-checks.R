test_that("Inner checks are triggered", {
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

  # Won't converge
  model_failed_converge <- model_random
  model_failed_converge$converge <- FALSE

  target <- "gen"

  expect_error(H2(model = c(model_random, model_fixed), target = target))
  expect_warning(H2(model_failed_converge, target = target, method = "Oakey"))
  expect_error(H2(model = model_random, target = "tamago"))

  # Target level
  expect_false(check_target_random(model_fixed, target))
  expect_true(check_target_random(model_random, target))

  # Method level
  expect_message(H2(model_fixed, target, "Oakey"))
  expect_message(H2(model_fixed, target, "Cullis"))
})
