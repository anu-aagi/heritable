test_that("narrowsense heritability works", {

  model <- asreml::asreml(y ~ rep,
                          random=~ vm(gen, lettuce_GRM),
                          data = lettuce_phenotypes |>
                            subset(loc == "L2"))

  # should be 0.586
  h2_Oakey(model, target = "gen")

})
