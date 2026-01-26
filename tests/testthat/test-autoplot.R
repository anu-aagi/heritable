test_that("Autoplot", {
  skip_if_not_installed("asreml")
  skip_on_ci()
  skip_on_cran()

  asreml_model_random <- readRDS(test_path("fixtures/asreml_model_random.rds"))

  # Extract the random effects
  expect_named(extract_var_comps(asreml_model_random))

  # Visual testing
  vdiffr::expect_doppelganger("Bar plot", autoplot(asreml_model_random))
  vdiffr::expect_doppelganger("Stacked bar plot", autoplot(asreml_model_random, type = "stacked"))
  vdiffr::expect_doppelganger("Pie plot", autoplot(asreml_model_random, type = "pie"))
  vdiffr::expect_doppelganger("Donut plot", autoplot(asreml_model_random, type = "donut"))
})
