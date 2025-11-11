test_that("Reproduce lme4 H2", {
  dat <- agridat::john.alpha

  # random genotype effect
  model_ran_lmer <- lme4::lmer(data = dat, formula = yield ~ rep + (1 | gen) + (1 | rep:block))

  expect_equal(unname(H2(model_ran_lmer, target = "gen", method = "Cullis")), 0.8091338, tolerance = 1e-7)
  expect_equal(unname(H2(model_ran_lmer, target = "gen", method = "Oakey")), 0.8091338, tolerance = 1e-7)
})


test_that("Reproduce lme4 Piepho", {
  requireNamespace("lme4", quietly = TRUE)

  dat <- agridat::john.alpha

  # random genotype effect
  model_ran <- lme4::lmer(data = dat, formula = yield ~ rep + (1 | gen) + (1 | rep:block))

  # fixed genotype effect
  g_fix <- lme4::lmer(
    data = dat,
    formula = yield ~ rep + gen + (1 | rep:block)
  )

  # Extract vc_g
  vc <- model_ran |>
    lme4::VarCorr() |>
    as.data.frame() # extract estimated variance components (vc)
  vc_g <- subset(vc, grp == "gen")$vcov # genotypic vc

  # Calculate mean variance of a difference between genotypes
  dBLUE <- emmeans::emmeans(g_fix, pairwise ~ gen)$contrasts |> as.data.frame()
  dBLUE$var <- dBLUE$SE^2

  vd_BLUE_avg <- mean(dBLUE$var) # 0.07295899
  # Fonti able to reproduce no rounding error

  # H2 Piepho ---------------------------------------------------------------
  H2_Piepho <- vc_g / (vc_g + vd_BLUE_avg / 2)
  H2_Piepho # 0.7966375
  # Fonti able to reproduce no rounding error
})

test_that("H2.lmerMod gives expected output", {
  requireNamespace("lme4", quietly = TRUE)

  dat <- agridat::john.alpha

  # random genotype effect
  model <- lme4::lmer(data = dat, formula = yield ~ rep + (1 | gen) + (1 | rep:block))

  output <- H2(model, method = "Naive", target = "gen")

  expect_length(output, 1) # if you expect one element
  expect_type(output, "double") # numeric type
  expect_equal(unname(output), 0.6364804, tolerance = 1e-7)
})
