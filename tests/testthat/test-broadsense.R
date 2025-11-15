test_that("heritability for asreml works", {
  skip_if_not_installed("asreml")
  skip_on_cran()


# random gen --------------------------------------------------------------
  truth_random_asreml <- c(
    "Cullis" = 0.8090841,
    "Oakey" = 0.8090728,
    "Piepho" = 0.8029759,
    "Delta" = 0.8090841,
    "Naive" = 0.8400648
  )

  model <- asreml::asreml(yield ~ rep,
                          random = ~ gen + rep:block,
                          data = agridat::john.alpha,
                          trace = FALSE)
  expect_equal(H2(model, target = "gen"), truth_random_asreml, tolerance = 1e-5)

  # NOTE: in theory these values should be the same
  # differences seem to arise because of different estimates of VCOMP
  truth_random_lme4 <- c(
    "Cullis" = 0.8091339,
    "Oakey" = 0.8091339,
    "Piepho" = 0.7966376,
    "Delta" = truth_asreml[["Delta"]],
    "Naive" = truth_asreml[["Naive"]]
  )

  model <- lme4::lmer(data = dat, formula = yield ~ rep + (1 | gen) + (1 | rep:block))
  expect_equal(H2(model, target = "gen"), truth_random_lme4, tolerance = 1e-4)

  # fixed gen ---------------------------------------------------------------
  truth_fixed_asreml <- c(
    "Cullis" = NA,
    "Oakey" = NA,
    "Piepho" = truth_random_asreml[["Piepho"]],
    "Delta" = 0.8030227, # arithmetic mean
    "Naive" = NA
  )

  model <- asreml::asreml(yield ~ rep + gen,
                          random = ~ rep:block,
                          data = agridat::john.alpha,
                          trace = FALSE
  )
  expect_equal(H2(model, target = "gen"), truth_fixed_asreml, tolerance = 1e-4)

  truth_fixed_lme4 <- c(
    "Cullis" = NA,
    "Oakey" = NA,
    "Piepho" = truth_random_lme4[["Piepho"]],
    "Delta" = 0.7967, # arithmetic mean # CHECK!
    "Naive" = NA
  )


  model <- lme4::lmer(data = dat, formula = yield ~ rep + gen + (1 | rep:block))
  expect_equal(H2(model, target = "gen"), truth_fixed_lme4, tolerance = 1e-5)


# GxE models --------------------------------------------------------------
  model <- asreml::asreml(yield ~ year:loc:trial,
                          random = ~ gen + gen:loc,
                          data = agridat::adugna.sorghum |>
                            transform(year = as.factor(year)),
                          trace = FALSE
  )
  H2(model, target = "gen")




})
