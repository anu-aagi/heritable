test_that("h2 heritability works", {
  # skip_if_not_installed("asreml")
  # skip_on_ci()
  # skip_on_cran()
  #
  # asreml_model_random <- readRDS(file = test_path("fixtures/asreml_model_random.rds"))
  # asreml_model_grm <- readRDS(file = test_path("fixtures/asreml_model_grm.rds"))
  #
  # # Should be 0.871 -- getting right
  # expect_equal(h2_Delta(asreml_model_grm, target = "gen", type = "BLUP"), 0.871, tolerance = 1e-3)
  #
  # # TODO: Why are values not matching? How do I know 0.818 is correct?
  # # should be 0.818 -- getting 0.8459519
  # # expect_lt(h2_Delta(asreml_model_grm, target = "gen"), H2_Delta(asreml_model_random, target = "gen"))
  # expect_lt(h2_Delta(asreml_model_grm, target = "gen", type = "BLUE"), h2_Delta(asreml_model_grm, target = "gen", type = "BLUP"))
  # # should be 0.586 / getting 0.5788
  # expect_lt(h2_Oakey(asreml_model_grm, target = "gen"), H2_Oakey(asreml_model_random, target = "gen"))
  #
  # expect_named(h2(asreml_model_grm, target = "gen"), c("Oakey", "Delta"))
  #
  # expect_s4_class(h2_Delta_pairwise(asreml_model_grm, target = "gen", type = "BLUP"), "dspMatrix")
  # expect_type(h2_Delta_by_genotype(asreml_model_grm, target = "gen", type = "BLUP"), "list")
  #
  # expect_named(h2(asreml_model_grm, target = "gen"), c("Oakey", "Delta"))
  # expect_length(h2(asreml_model_grm, target = "gen"), 2)
  #
  # expect_named(h2_Delta_by_genotype(asreml_model_grm, target = "gen", type = "BLUP"), levels(asreml_model_grm$mf$gen), ignore.order = TRUE)
  # expect_equal(nrow(h2_Delta_pairwise(asreml_model_grm, target = "gen", type = "BLUP")), length(levels(asreml_model_grm$mf$gen)))
})
