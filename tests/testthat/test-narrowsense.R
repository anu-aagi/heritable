test_that("narrowsense heritability works", {

  model <- asreml::asreml(y ~ rep,
                          random=~ vm(gen, lettuce_GRM),
                          data = lettuce_phenotypes |>
                            subset(loc == "L2"))

  # should be 0.586 / getting 0.5788
  h2_Oakey(model, target = "gen")

  # should be 0.871 -- getting right
  h2_Delta(model, target = "gen", type = "BLUP")
  # should be 0.818 -- getting 0.8459519
  h2_Delta(model, target = "gen", type = "BLUE")

  h2_Delta_pairwise(model, target = "gen", type = "BLUP")
  h2_Delta_pairwise(model, target = "gen", type = "BLUE")

  h2_Delta_by_genotype(model, target = "gen", type = "BLUP")
})
