test_that("OAKEY heritability estimation works", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  skip_on_ci()

  asreml_model_random <- readRDS(test_path("fixtures/asreml_model_random.rds"))
  lettuce_asreml <- readRDS(test_path("fixtures/lettuce_asreml.rds"))

  # Method implemented by YD as per Schmid et al. 2009
  H2_Oakey(asreml_model_random, target = "gen")

  # Method implemented by ET as per Oakey et al. 2006
  H2_Oakey_old <- function(model, target) {
  n_g <- model$noeff[[target]]
  vc_g <- model$vparameters[[target]] * model$sigma2
  vcov_g <- predict(model,
    classify = target,
    only = target,
    vcov = TRUE,
    trace = FALSE
  )$vcov

  Gg_inv <- diag(1 / vc_g, nrow = n_g, ncol = n_g)
  M <- diag(n_g) - (Gg_inv %*% vcov_g)
  eM <- eigen(M)
  thres <- 1e-5

  mean(eM$values[eM$values > thres])
  }

  expect_equal(
    H2_Oakey(asreml_model_random, target = "gen"),
    H2_Oakey_old(asreml_model_random, target = "gen"),
    tolerance = 1e-3
  )
})
