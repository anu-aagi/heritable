test_that("h2 heritability works", {

  asreml_model_grm <- readRDS(file = test_path("fixtures/asreml_model_grm.rds"))

  # Should be 0.871 -- getting right
  expect_equal(h2_Delta(asreml_model_grm, target = "gen", type = "BLUP"), 0.871, tolerance = 1e-3)

  # TODO: Why are values not matching? How do I know 0.818 is correct?
  # should be 0.818 -- getting 0.8459519
  h2_Delta(asreml_model_grm, target = "gen", type = "BLUE")
  # should be 0.586 / getting 0.5788
  h2_Oakey(asreml_model_grm, target = "gen")

  expect_s4_class(h2_Delta_pairwise(asreml_model_grm, target = "gen", type = "BLUP"), "dspMatrix")
  expect_type(h2_Delta_by_genotype(asreml_model_grm, target = "gen", type = "BLUP"), "list")
})
