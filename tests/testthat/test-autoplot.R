test_that("Autoplot", {
  skip_if_not_installed("asreml")
  skip_on_ci()
  skip_on_cran()

  asreml_model_random <- readRDS(test_path("fixtures/asreml_model_random.rds"))

  # Extract the random effects
  H2 <- H2(asreml_model_random, "gen", "Standard")
  expect_named(extract_var_comps(H2))

  # Visual testing
  vdiffr::expect_doppelganger("Bar plot", autoplot(H2))
  vdiffr::expect_doppelganger("Stacked bar plot", autoplot(H2, type = "stacked"))
  vdiffr::expect_doppelganger("Pie plot", autoplot(H2, type = "pie"))
  vdiffr::expect_doppelganger("Donut plot", autoplot(H2, type = "donut"))
})
