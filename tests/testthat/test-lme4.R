test_that("Reproduce lme4 H2", {
  dat <- agridat::john.alpha

  # random genotype effect
  model_ran_lmer <- lme4::lmer(data = dat, formula = yield ~ rep + (1 | gen) + (1 | rep:block))
  model_fix_lmer <- fit_counterpart_model(model_ran_lmer, target)

  expect_equal(unname(H2(model_ran_lmer, target = "gen", method = "Cullis")), 0.8091338, tolerance = 1e-7)
  expect_equal(unname(H2(model_ran_lmer, target = "gen", method = "Oakey")), 0.8091338, tolerance = 1e-7)
  expect_equal(unname(H2(model_ran_lmer, target = "gen", method = "Naive")), 0.6364804, tolerance = 1e-7)
  expect_equal(unname(H2(model_ran_lmer, target = "gen", method = "Piepho")), 0.7966375, tolerance = 1e-7)
  expect_named(H2(model_fix_lmer, target = "gen", method = "Delta"), "Delta")
  expect_named(H2(model_ran_lmer, target = "gen", method = "Delta"), "Delta")
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
  deltas <- emmeans::emmeans(model_fix, specs = as.formula(paste("pairwise ~", target)))$contrasts |> as.data.frame()
  deltas$var <- deltas$SE^2 # Get variance
  # Take contrasts and turn into variance-covariance matrix
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
  for (i in 1:nrow(deltas)) {
    # Extract genotype names from contrast column
    pair <- strsplit(as.character(deltas$contrast[i]), " - ")[[1]]
    g1 <- pair[1]
    g2 <- pair[2]
    
    # Variance of difference: Var(g1 - g2) = Var(g1) + Var(g2) - 2*Cov(g1, g2)
    # Get covariance between g1 and g2 (0 by default, but can be specified)
    Vd_g[g1, g2] <- deltas$var[i] - 2 * cov_g[g1, g2]
    Vd_g[g2, g1] <- deltas$var[i] - 2 * cov_g[g2, g1] # symmetric
  }
  Vd_g

  # H2 Delta ---------------------------------------------------------------
    BLUE_H2D_ij <- H2_Delta_BLUE_parameters(vc_g, cov = 0, Vd_g)
})

test_that("Implementing H2 BLUPS Delta",{
  dat <- agridat::john.alpha

  # random genotype effect
  model <- lme4::lmer(data = dat, formula = yield ~ rep + (1 | gen) + (1 | rep:block))
  
  vc <- lme4::VarCorr(model)
  ngrps <- lme4::ngrps(model)
  # Note the index and kronecker order needs to be followed careful downstream
  Glist <- lapply(names(vc), function(agrp) {
    Matrix::kronecker(vc[[agrp]], diag(ngrps[[agrp]]))
  })
  G <- do.call(Matrix::bdiag, Glist)

  n <- nrow(model@frame)
  R <- diag(n) * sigma(model)^2

  X <- as.matrix(lme4::getME(model, "X"))
  Z <- as.matrix(lme4::getME(model, "Z"))

  C11 <- t(X) %*% solve(R) %*% X
  C12 <- t(X) %*% solve(R) %*% Z
  C21 <- t(Z) %*% solve(R) %*% X
  C22 <- t(Z) %*% solve(R) %*% Z + solve(G)

  C <- rbind(
    cbind(C11, C12),
    cbind(C21, C22)
  )
  C_inv <- solve(C)
  gnames <- levels(model@flist[[target]])
  C22_g <- C_inv[gnames, gnames]
  
  n_g <- ngrps[[target]]
  vc_g <- vc[[target]][1]

# Compute variance of difference from PEV
var_diff_matrix <- outer(
  1:nrow(C22_g), 1:ncol(C22_g),
  Vectorize(function(i, j) C22_g[i, i] + C22_g[j, j] - 2 * C22_g[i, j])
)

diag(var_diff_matrix) <- NA
rownames(var_diff_matrix) <- colnames(var_diff_matrix) <- rownames(C22_g)

#Compute variance of difference between two genotypes
var_diff <- function(i, j, C22) {
  C22[i, i] + C22[j, j] - 2 * C22[i, j]
}

var_diffs <- apply(pairs, 2, function(p) var_diff(p[1], p[2], C22_g))

# Put in a tidy data frame
var_diffs_df <- data.frame(
  g1 = pairs[1, ],
  g2 = pairs[2, ],
  var_diff = var_diffs
)

expect_equal(var_diff_matrix[1,2], var_diffs_df$var_diff[1])

H2_Delta_BLUP_parameters(vc_g = vc_g, cov = 0, vd_BLUP_matrix = var_diff_matrix)


})