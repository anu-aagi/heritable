test_that("H2 for asreml works", {
  skip_if_not_installed("asreml")
  skip_on_cran()

  # ASREML Random "gen" --------------------------------------------------------
  asreml_model_random <- readRDS(test_path("fixtures/asreml_model_random.rds"))

  truth_random_asreml <- c(
    "Cullis" = 0.8090841,
    "Oakey" = 0.8090728,
    "Piepho" = 0.8029759,
    "Delta" = 0.8090841,
    "Standard" = 0.8400648
  )

  asreml_H2_all_methods <- H2(asreml_model_random, target = "gen")

  expect_equal(setNames(as.numeric(asreml_H2_all_methods), names(asreml_H2_all_methods)), truth_random_asreml, tolerance = 1e-5)
  expect_equal(H2_Delta(asreml_model_random, target = "gen", type = "BLUE"), 0.8030227, tolerance = 1e-5)


  # Different residual structure -----------------------------------------------
  asreml_model_R <- readRDS(file = test_path("fixtures/asreml_model_R.rds"))

  expect_no_error(H2(asreml_model_R, "gen"))

  # GxE models -----------------------------------------------------------------
  #TODO: Currently not yet implemented so should error
  asreml_model_g_by_e <- readRDS(file = test_path("fixtures/asreml_model_g_by_e.rds"))

  expect_error(H2(asreml_model_g_by_e, target = "gen"))
})


test_that("H2 works for lme4",{
  # lme4 Random "gen" ----------------------------------------------------------
  lmer_model_random <- readRDS(test_path("fixtures/lmer_model_random.rds"))

  # NOTE: in theory these values should be the same
  # Differences seem to arise because of different estimates of VCOMP
  truth_random_asreml <- c(
    "Cullis" = 0.8090841,
    "Oakey" = 0.8090728,
    "Piepho" = 0.8029759,
    "Delta" = 0.8090841,
    "Standard" = 0.8400648
  )

  truth_random_lme4 <- c(
    "Cullis" = 0.8091339,
    "Oakey" = 0.8091339,
    "Piepho" = 0.7966376,
    "Delta" = truth_random_asreml[["Delta"]],
    "Standard" = truth_random_asreml[["Standard"]]
  )

  lmer_H2_all_methods <- H2(lmer_model_random, target = "gen")

  expect_equal(setNames(as.numeric(lmer_H2_all_methods), names(lmer_H2_all_methods)), truth_random_lme4, tolerance = 1e-4)
  expect_equal(H2_Delta(lmer_model_random, target = "gen", type = "BLUE"), 0.7967, tolerance = 1e-5)

  # GxE models -----------------------------------------------------------------
  #TODO: Currently not yet implemented so should error
  lmer_model_g_by_e <- readRDS(test_path("fixtures/lmer_model_g_by_e.rds"))

  expect_error(H2(lmer_model_g_by_e, target = "gen"))
})
