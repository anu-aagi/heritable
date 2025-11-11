test_that("Reproduce lme4 H2", {
  dat <- agridat::john.alpha

  # random genotype effect
  model_ran_lmer <- lme4::lmer(data = dat, formula = yield ~ rep + (1 | gen) + (1 | rep:block))

  expect_equal(unname(H2(model_ran_lmer, target = "gen", method = "Cullis")), 0.8091338, tolerance = 1e-7)
  expect_equal(unname(H2(model_ran_lmer, target = "gen", method = "Oakey")), 0.8091338, tolerance = 1e-7)
  expect_equal(unname(H2(model_ran_lmer, target = "gen", method = "Naive")), 0.6364804, tolerance = 1e-7)
  expect_equal(unname(H2(model_ran_lmer, target = "gen", method = "Piepho")), 0.7966375, tolerance = 1e-7)
})


test_that("Reproduce lme4 Piepho", {

  dat <- agridat::john.alpha
  target = "gen"
  model_fix <- lme4::lmer(data = dat, formula = yield ~ rep + gen + (1 | rep:block))
  model_ran <- fit_counterpart_model(model_fix, target)

 # Extract vc_g
  vc <- model_ran |>
    lme4::VarCorr() |>
    as.data.frame() # extract estimated variance components (vc)
  
  vc_g <- subset(vc, grp == "gen")$vcov # genotypic vc

  # Calculate mean variance of a difference between genotypes
  dBLUE <- emmeans::emmeans(model_fix, specs = as.formula(paste("pairwise ~", target)))$contrasts |> as.data.frame()
  dBLUE$var <- dBLUE$SE^2 # Get variance

  vd_BLUE_avg <- mean(dBLUE$var) # 0.07295899
  # Fonti able to reproduce no rounding error

  # H2 Piepho ---------------------------------------------------------------
  H2_Piepho <- vc_g / (vc_g + vd_BLUE_avg / 2)
  H2_Piepho <- H2_Piepho_parameters(vc_g, vd_BLUE_avg)
  H2_Piepho # 0.7966375
  # Fonti able to reproduce no rounding error
})