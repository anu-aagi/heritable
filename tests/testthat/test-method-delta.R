test_that("delta method works", {

  #############
  # Fit model #
  #############
  # Genotype as random effect
  model <- asreml::asreml(fixed = yield ~ rep,
                          random=       ~ gen + rep:block,
                          data=agridat::john.alpha)

  target <-  "gen"

  H2_Delta(model, target)

})
