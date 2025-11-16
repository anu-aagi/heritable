

#' @export
H2_Standard.lmerMod <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  # Get genotype variance
  vc <- lme4::VarCorr(model)
  vc_g <- vc[[target]][1]

  # Get residual variance
  vc_e <- stats::sigma(model)^2

  n_r <- table(model@flist[[target]])

  H2_Standard <- H2_Standard_parameters(vc_g, vc_e, n_r)

  return(H2_Standard)
}

#' @export
H2_Cullis.lmerMod <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  C_inv <- PEV_from_lme4(model)
  g <- geno_components_from_lme4(model, target, C_inv)

  one <- matrix(1, nrow = g$n_g, ncol = 1)
  P_mu <- diag(g$n_g, g$n_g) - one %*% t(one)
  vdBLUP_sum <- sum(diag(P_mu %*% g$C22_g))
  vdBLUP_avg <- vdBLUP_sum * (2 / (g$n_g * (g$n_g - 1)))

  H2_Cullis_parameters(vdBLUP_avg, g$vc_g)
}

PEV_from_lme4 <- function(model) {
  vc <- lme4::VarCorr(model)
  ngrps <- lme4::ngrps(model)
  # Note the index and kronecker order needs to be followed careful downstream
  Glist <- lapply(names(vc), function(agrp) {
    Matrix::kronecker(vc[[agrp]], diag(ngrps[[agrp]]))
  })
  G <- do.call(Matrix::bdiag, Glist)

  n <- nrow(model@frame)
  R <- diag(n) * stats::sigma(model)^2

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
  solve(C)
}

geno_components_from_lme4 <- function(model, target, C_inv) {
  vc <- lme4::VarCorr(model)
  ngrps <- lme4::ngrps(model)
  gnames <- levels(model@flist[[target]])
  C22_g <- C_inv[gnames, gnames]
  n_g <- ngrps[[target]]
  vc_g <- vc[[target]][1]
  list(n_g = n_g, vc_g = vc_g, C22_g = C22_g, gnames = gnames)
}

#' @export
H2_Oakey.lmerMod <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  C_inv <- PEV_from_lme4(model)
  g <- geno_components_from_lme4(model, target, C_inv)
  Gg_inv <- diag(1 / g$vc_g, nrow = g$n_g, ncol = g$n_g)

  return(H2_Oakey_parameters(Gg_inv, g$C22_g))
}

#' @export
H2_Piepho.lmerMod <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is fixed or random
  # If random, refit fixed model
  if(check_target_random(model, target)){
    model_ran <- model
    model_fix <- fit_counterpart_model(model, target)
  }
  # If fixed, refit random model
  else if(!check_target_random(model, target)){
    model_fix <- model
    model_ran <- fit_counterpart_model(model, target)
  }

  # Get genotype variance
  vc <- lme4::VarCorr(model_ran)
  vc_g <- vc[[target]][1]

 # Get mean variance of a difference between genotypes
  d_BLUE <- emmeans::emmeans(model_fix, specs = as.formula(paste("pairwise ~", target)))$contrasts |> as.data.frame()
  d_BLUE$var <- d_BLUE$SE^2 # Get variance

  vd_BLUE_avg <- mean(d_BLUE$var)

  # vc_g / (vc_g + vdBLUE_avg / 2)
  H2_Piepho <- H2_Piepho_parameters(vc_g, vd_BLUE_avg)

  return(H2_Piepho)
}

#' @export
H2_Delta_BLUE_pairwise.lmerMod <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is random or fixed
  # TODO: How to handle if both fixed and random?
  if(check_target_both(model, target)) {
    cli::cli_abort("The target {.var {target}} is fitted as both fixed and random effect")
  }

  # If fixed, compute H2 delta with BLUES
  if(!check_target_random(model, target)) {

    model_fix <- model
    model_ran <- fit_counterpart_model(model, target)

    # Extract vc_g and vc_e
    vc <- lme4::VarCorr(model_ran)
    vc_g <- vc[[target]][1]

    # Calculate mean variance of a difference between genotypes
    deltas <- emmeans::emmeans(model_fix, specs = as.formula(paste("pairwise ~", target)))$contrasts |> as.data.frame()
    deltas$var <- deltas$SE^2 # Get variance

    # Take pairwise differences and turn into variance-covariance matrix
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

    # Fill in the pairwise variances from deltas
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
    } else if(check_target_random(model, target)) {
      # Abort and tell user to compute Delta with BLUPs
      cli::cli_abort("The target {.var {target}} is fitted as a random effect. See H2_Delta_BLUP")
    }

  # H2 Delta BLUE
  H2_Delta_BLUE <- H2_Delta_BLUE_parameters(vc_g, vc_g, cov = 0, Vd_g)

  return(H2_Delta_BLUE)
}

#' @export
H2_Delta_BLUP_pairwise.lmerMod <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is random or fixed
  # TODO: How to handle if both fixed and random?
  if(check_target_both(model, target)) {
    cli::cli_abort("The target {.var {target}} is fitted as both fixed and random effect")
  }

  if(check_target_random(model, target)){
    C_inv <- PEV_from_lme4(model)
    g <- geno_components_from_lme4(model, target, C_inv)

    # Compute variance of difference from PEV
    Vd_g <- outer(
      1:nrow(g$C22_g), 1:ncol(g$C22_g),
      Vectorize(function(i, j) g$C22_g[i, i] + g$C22_g[j, j] - 2 * g$C22_g[i, j])
    )

    diag(Vd_g) <- NA
    dimnames(Vd_g) <- list(g$gnames, g$gnames)

  } else if(!check_target_random(model, target)) {
    # Abort and tell user to compute Delta with BLUES
    cli::cli_abort("The target {.var {target}} is fitted as a fixed effect. See H2_Delta_BLUE.")
  }

  # H2 Delta BLUP
  H2_Delta_BLUP <- H2_Delta_BLUP_parameters(g$vc_g, g$vc_g, cov = 0, Vd_g)

  row.names(H2_Delta_BLUP) <- rownames(Vd_g)
  colnames(H2_Delta_BLUP) <- colnames(Vd_g)

  H2_Delta_BLUP
}


#' @export
H2_Delta_pairwise.lmerMod <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is random or fixed
  if(!check_target_random(model, target)) {
    H2_Delta <- H2_Delta_BLUE_pairwise.lmerMod(model, target, options)
  } else if(check_target_random(model, target)) {
    H2_Delta <- H2_Delta_BLUP_pairwise.lmerMod(model, target, options)
  }

  return(H2_Delta)
}

