test_that("h2 heritability works", {
  skip_if_not_installed("asreml")
  skip_on_ci()
  skip_on_cran()
  
  asreml_model_random <- readRDS(file = test_path("fixtures/asreml_model_random.rds"))
  asreml_model_grm <- readRDS(file = test_path("fixtures/asreml_model_grm.rds"))

  # From Schmidt et al 2019 Fig 1. 
  # Delta 
  expect_equal(h2_Delta(asreml_model_grm, target = "gen", type = "BLUP"), 0.871, tolerance = 1e-3)
  #TODO: h2 BLUE values not matching
  expect_equal(h2_Delta(asreml_model_grm, target = "gen", type = "BLUE"), 0.818, tolerance = 1e-3)
  
  # Oakey
  #TODO: h2 Oakey values not matching
  expect_equal(h2_Oakey(asreml_model_grm, target = "gen"), 0.586, tolerance = 1e-3)
  
  # Structural checks
  expect_named(h2(asreml_model_grm, target = "gen"), c("Oakey", "Delta"))
  expect_s4_class(h2_Delta_pairwise(asreml_model_grm, target = "gen", type = "BLUP"), "dspMatrix")
  expect_type(h2_Delta_by_genotype(asreml_model_grm, target = "gen", type = "BLUP"), "list")
  expect_named(h2(asreml_model_grm, target = "gen"), c("Oakey", "Delta"))
  expect_length(h2(asreml_model_grm, target = "gen"), 2)
  expect_named(h2_Delta_by_genotype(asreml_model_grm, target = "gen", type = "BLUP"), levels(asreml_model_grm$mf$gen), ignore.order = TRUE)
  expect_equal(nrow(h2_Delta_pairwise(asreml_model_grm, target = "gen", type = "BLUP")), length(levels(asreml_model_grm$mf$gen)))
})
