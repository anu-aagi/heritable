test_that("heritability for asreml works", {
  skip_if_not_installed("asreml")
  skip_on_cran()


# random gen --------------------------------------------------------------
  truth_random_asreml <- c(
    "Cullis" = 0.8090841,
    "Oakey" = 0.8090728,
    "Piepho" = 0.8029759,
    "Delta" = 0.8090841,
    "Standard" = 0.8400648
  )

  model <- asreml::asreml(yield ~ rep,
                          random = ~ gen + rep:block,
                          data = agridat::john.alpha,
                          trace = FALSE)
  expect_equal(H2(model, target = "gen"), truth_random_asreml, tolerance = 1e-5)
  expect_equal(H2_Delta(model, target = "gen", type = "BLUE"), 0.8030227, tolerance = 1e-5)

  # NOTE: in theory these values should be the same
  # differences seem to arise because of different estimates of VCOMP
  truth_random_lme4 <- c(
    "Cullis" = 0.8091339,
    "Oakey" = 0.8091339,
    "Piepho" = 0.7966376,
    "Delta" = truth_random_asreml[["Delta"]],
    "Standard" = truth_random_asreml[["Standard"]]
  )

  model <- lme4::lmer(data = dat, formula = yield ~ rep + (1 | gen) + (1 | rep:block))
  expect_equal(H2(model, target = "gen"), truth_random_lme4, tolerance = 1e-4)
  expect_equal(H2_Delta(model, target = "gen", type = "BLUE"), 0.7967, tolerance = 1e-5)


# Different R

  model <- asreml::asreml(yield ~ 1,
                        random=~ gen + col + row,
                        residual=~ar1(col):ar1(row),
                        data = agridat::gilmour.serpentine |>
                          transform(col = as.factor(col),
                                    row = as.factor(row)))
  H2(model, "gen")



# GxE models --------------------------------------------------------------
  model <- asreml::asreml(yield ~ year:loc:trial,
                          random = ~ gen + gen:loc,
                          data = agridat::adugna.sorghum,
                          trace = FALSE
  )
  expect_error(H2(model, target = "gen"))

  model <- lme4::lmer(yield ~ year:loc:trial + (1|gen) + (1|gen:loc),
                      data = agridat::adugna.sorghum)

  expect_error(H2(model, target = "gen"))



})
