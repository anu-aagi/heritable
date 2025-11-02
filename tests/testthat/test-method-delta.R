test_that("delta method works", {

  #############
  # Fit model #
  #############
  # Genotype as random effect
  model <- asreml::asreml(fixed = yield ~ rep,
                          random= ~ gen + rep:block,
                          data= agridat::john.alpha,
                          trace = FALSE)

  target <-  "gen"

  H2(model, target, method = "Delta")
})
