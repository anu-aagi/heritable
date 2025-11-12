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

test_that("Implementing H2 BLUES Delta", {
  dat <- agridat::john.alpha
  target = "gen"
  model_fix <- lme4::lmer(data = dat, formula = yield ~ rep + gen + (1 | rep:block))
  model_ran <- fit_counterpart_model(model_fix, target)

 # Extract vc_g and vc_e
  vc <- lme4::VarCorr(model_ran)
  vc_g <- vc[[target]][1]

  # Calculate mean variance of a difference between genotypes
  dBLUE <- emmeans::emmeans(model_fix, specs = as.formula(paste("pairwise ~", target)))$contrasts |> as.data.frame()
  # Take contrasts and turn into variance-covariance matrix

    # Convert pairwise differences to variance-covariance matrix
  # Extract genotype levels
  lev_g <- levels(model_fix@frame[[target]])
  n_g <- length(lev_g)
  
  # Create variance-covariance matrix for genotypes (H2: covariance = 0)
  # TODO: For narrow sense, we will need to replace this from the kinship matrix
  cov_g <- matrix(0, nrow = n_g, ncol = n_g)
  diag(cov_g) <- vc_g  # Set diagonal to genotype variance
  dimnames(cov_g) <- list(lev_g, lev_g)

  # Start with empty variance matrix for differences 
  Vd_g <- matrix(0, nrow = n_g, ncol = n_g)
  dimnames(Vd_g) <- list(lev_g, lev_g)
  
    # Fill in the pairwise variances from dBLUE
  for (i in 1:nrow(dBLUE)) {
    # Extract genotype names from contrast column
    pair <- strsplit(as.character(dBLUE$contrast[i]), " - ")[[1]]
    g1 <- pair[1]
    g2 <- pair[2]
    
    # Variance of difference: Var(g1 - g2) = Var(g1) + Var(g2) - 2*Cov(g1, g2)
    # Get covariance between g1 and g2 (0 by default, but can be specified)
    Vd_g[g1, g2] <- dBLUE$var[i] - 2 * cov_g[g1, g2]
    Vd_g[g2, g1] <- dBLUE$var[i] - 2 * cov_g[g1, g2] # symmetric
  }
  }

  
  # H2 Delta ---------------------------------------------------------------
    BLUE_H2D_ij <- H2_Delta_BLUE_parameters(vc_g, cov = 0, Vd_g)
    BLUP_H2D_ij <- H2_Delta_BLUP_parameters(vc_g, cov = 0, Vd_g)
})