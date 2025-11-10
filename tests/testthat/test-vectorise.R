test_that("H2 can handle multiple methods", {
  skip_if_not_installed("asreml")
  skip_on_cran()

  # Setup test model
  model <- asreml::asreml(yield ~ rep,
    random = ~ gen + rep:block,
    data = agridat::john.alpha,
    trace = FALSE
  )

  # Test multiple methods
  res_multi <- H2(model, target = "gen", method = c("Cullis", "Oakey"))
  expect_named(res_multi, c("Cullis", "Oakey"))
  expect_equal(res_multi[["Cullis"]], 0.8090841, tolerance = 1e-7)
  expect_equal(res_multi[["Oakey"]], 0.8090728, tolerance = 1e-7)
  
  # Test all methods
  res_all <- H2(model, target = "gen", 
                method = c("Cullis", "Oakey", "Delta", "Piepho", "Naive"))
  expect_named(res_all, c("Cullis", "Oakey", "Delta", "Piepho", "Naive"))
  expect_length(res_all, 5)
  expect_type(res_all, "double")
})