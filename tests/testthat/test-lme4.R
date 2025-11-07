test_that("Reproduce lme4 Cullis", {
  requireNamespace("lme4", quietly = TRUE)
  
  dat <- agridat::john.alpha
  
  # random genotype effect
  g.ran <- lme4::lmer(data = dat, formula = yield ~ rep + (1|gen) + (1|rep:block))

  # handle model estimates --------------------------------------------------
  # to my knowledge, lme4 does not offer a function to
  # extract variance-covariance-matrices for BLUPs (a.k.a. prediction error variance [PEV] matrix).
  # therefore, I here manually reconstruct mixed model equation for this specific example.
  # notice that this solution therefore only works for this specific model!

  vc <- g.ran |> lme4::VarCorr() |> as.data.frame() # extract estimated variance components (vc)

  # R = varcov-matrix for error term
  n    <- g.ran@frame |> nrow() # numer of observations
  vc_e <- subset(vc, grp=="Residual")$vcov     # error 
  R    <- diag(n)*vc_e                                      # R matrix = I_n * vc_e

  # G = varcov-matrx for all random effects
  # subset of G regarding genotypic effects
  n_g  <- summary(g.ran)$ngrps["gen"] # number of genotypes
  vc_g <- subset(vc, grp=="gen")$vcov # genotypic vc
  G_g  <- diag(n_g)*vc_g   # gen part of G matrix = I * vc.g

  # subset of G regarding incomplete block effects
  n_b  <- summary(g.ran)$ngrps["rep:block"] # number of incomplete blocks
  vc_b <- subset(vc, grp=="rep:block")$vcov # incomplete block vc
  G_b  <- diag(n_b)*vc_b  # incomplete block part of G matrix = I * vc.b

  G <- bdiag(G_g, G_b) # G is blockdiagonal with G.g and G.b in this example

  # Design Matrices
  X <- g.ran |> lme4::getME("X") |> as.matrix() # Design matrix fixed effects
  Z <- g.ran |> lme4::getME("Z") |> as.matrix() # Design matrix random effects

  # Mixed Model Equation (HENDERSON 1986; SEARLE et al. 2006)
  C11 <- t(X) %*% solve(R) %*% X
  C12 <- t(X) %*% solve(R) %*% Z
  C21 <- t(Z) %*% solve(R) %*% X
  C22 <- t(Z) %*% solve(R) %*% Z + solve(G) 

  C <- rbind(cbind(C11, C12),  
            cbind(C21, C22)) |> as.matrix() # Combine components into one matrix C

  # Mixed Model Equation Solutions 
  C_inv <- C |> solve()                             # Inverse of C
  C22_g <- C_inv[levels(dat$gen), levels(dat$gen)] # subset of C.inv that refers to genotypic BLUPs

  # Mean variance of BLUP-difference from C22 matrix of genotypic BLUPs
  one        <- matrix(1, nrow=n_g, ncol=1)      # vector of 1s
  P_mu       <- diag(n_g, n_g) - one %*% t(one)  # P_mu = matrix that centers for overall-mean
  vdBLUP_sum <- psych::tr(P_mu %*% C22_g)        # sum of all variance of differences = trace of P_mu*C22_g
  vdBLUP_avg <- vdBLUP_sum * (2/(n_g*(n_g-1)))   # mean variance of BLUP-difference = divide sum by number of genotype pairs


  # H2 Cullis ---------------------------------------------------------------
  H2Cullis <- 1 - (vdBLUP_avg / 2 / vc_g)
  H2Cullis #0.8091336
  # Fonti got 0.8091338, might be a rounding difference?
})


test_that("Reproduce lme4 Oakey", {
  requireNamespace("lme4", quietly = TRUE)
  
  dat <- agridat::john.alpha
  
  # random genotype effect
  g.ran <- lme4::lmer(data = dat, formula = yield ~ rep + (1|gen) + (1|rep:block))

  # Everything up to C22_g is the same as in Cullis test
    requireNamespace("lme4", quietly = TRUE)
  
  dat <- agridat::john.alpha
  
  # random genotype effect
  g.ran <- lme4::lmer(data = dat, formula = yield ~ rep + (1|gen) + (1|rep:block))

  # handle model estimates --------------------------------------------------
  # to my knowledge, lme4 does not offer a function to
  # extract variance-covariance-matrices for BLUPs (a.k.a. prediction error variance [PEV] matrix).
  # therefore, I here manually reconstruct mixed model equation for this specific example.
  # notice that this solution therefore only works for this specific model!

  vc <- g.ran |> lme4::VarCorr() |> as.data.frame() # extract estimated variance components (vc)

  # R = varcov-matrix for error term
  n    <- g.ran@frame |> nrow() # numer of observations
  vc_e <- subset(vc, grp=="Residual")$vcov     # error 
  R    <- diag(n)*vc_e                                      # R matrix = I_n * vc_e

  # G = varcov-matrx for all random effects
  # subset of G regarding genotypic effects
  n_g  <- summary(g.ran)$ngrps["gen"] # number of genotypes
  vc_g <- subset(vc, grp=="gen")$vcov # genotypic vc
  G_g  <- diag(n_g)*vc_g   # gen part of G matrix = I * vc.g

  # subset of G regarding incomplete block effects
  n_b  <- summary(g.ran)$ngrps["rep:block"] # number of incomplete blocks
  vc_b <- subset(vc, grp=="rep:block")$vcov # incomplete block vc
  G_b  <- diag(n_b)*vc_b  # incomplete block part of G matrix = I * vc.b

  G <- bdiag(G_g, G_b) # G is blockdiagonal with G.g and G.b in this example

  # Design Matrices
  X <- g.ran |> lme4::getME("X") |> as.matrix() # Design matrix fixed effects
  Z <- g.ran |> lme4::getME("Z") |> as.matrix() # Design matrix random effects

  # Mixed Model Equation (HENDERSON 1986; SEARLE et al. 2006)
  C11 <- t(X) %*% solve(R) %*% X
  C12 <- t(X) %*% solve(R) %*% Z
  C21 <- t(Z) %*% solve(R) %*% X
  C22 <- t(Z) %*% solve(R) %*% Z + solve(G) 

  C <- rbind(cbind(C11, C12),  
            cbind(C21, C22)) |> as.matrix() # Combine components into one matrix C

  # Mixed Model Equation Solutions 
  C_inv <- C |> solve()                             # Inverse of C
  C22_g <- C_inv[levels(dat$gen), levels(dat$gen)] # subset of C.inv that refers to genotypic BLUPs

  ED    <- diag(n_g) - (solve(G_g) %*% C22_g)       # [see p. 813 bottom left in Oakey (2006)]
  eM    <- ED |> eigen()                             # obtain eigenvalues

  # H2 Oakey ----------------------------------------------------------------
  # main method [see eq. (7) in Oakey (2006)]
  H2Oakey <- sum(eM$values)/(n_g-1) 
  H2Oakey # 0.8091336
  # Fonti got 0.8091338, might be a rounding difference

  # approximate method [see p. 813 top right in Oakey (2006)]
  H2Oakey.approx <- 1 - psych::tr( as.matrix(solve(G_g) %*% C22_g / n_g ) ) # or sum(diag(m), na.rm = TRUE)
  H2Oakey.approx # 0.7754197
  # Fonti able to reproduce 0.7754199 with some rounding error

})
