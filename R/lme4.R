#' @noRd
#' @keywords internal
geno_components_from_lme4 <- function(model, target, calc_C22 = TRUE) {
  X <- as.matrix(lme4::getME(model, "X"))
  Z <- as.matrix(lme4::getME(model, "Z"))

  sigma2 <- sigma(model)^2
  Lambda <- lme4::getME(model, "Lambda")
  G <- tcrossprod(Lambda) * sigma2
  dimnames(G) <- list(colnames(Z), colnames(Z))

  gnames <- levels(model@flist[[target]])
  vc_g <- G[gnames, gnames, drop=FALSE]
  n_g <- length(gnames)

  if(calc_C22){
    R <- diag(nrow(X)) * sigma2
    V <- R + Z %*% G %*% t(Z)
    Vinv <- solve(V)
    P <- Vinv - Vinv %*% X %*% solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv
    C22 <- G - G %*% t(Z) %*% P %*% Z %*% G
    dimnames(C22) <- list(colnames(Z), colnames(Z))
    C22_g <- C22[gnames, gnames, drop=FALSE]
  } else {
    C22_g <- NULL
  }

  list(n_g = n_g, vc_g = vc_g, C22_g = C22_g, gnames = gnames)
}


#' @noRd
#' @export
H2_Standard.lmerMod <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  # Get genotype variance
  vc_g <- geno_components_from_lme4(model, target, calc_C22 = FALSE)$vc_g
  vc_g <- mean(Matrix::diag(vc_g))

  # Get residual variance
  vc_e <- stats::sigma(model)^2

  n_r <- table(model@flist[[target]])

  H2_Standard <- H2_Standard_parameters(vc_g, vc_e, n_r)

  return(H2_Standard)
}

#' @noRd
#' @export
H2_Cullis.lmerMod <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  g <- geno_components_from_lme4(model, target)
  vc_g <- mean(Matrix::diag(g$vc_g))

  P_mu <- Matrix::Diagonal(n = g$n_g, x = g$n_g) - 1
  vdBLUP_sum <- sum(Matrix::diag(P_mu %*% g$C22_g))
  vdBLUP_avg <- vdBLUP_sum * (2 / (g$n_g * (g$n_g - 1)))

  return(H2_Cullis_parameters(vdBLUP_avg, vc_g))
}

#' @noRd
#' @export
H2_Oakey.lmerMod <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }
  g <- geno_components_from_lme4(model, target)
  Gg_inv <- Matrix::chol2inv(chol(g$vc_g))

  return(H2_Oakey_parameters(Gg_inv, g$C22_g))
}

#' @noRd
#' @export
H2_Piepho.lmerMod <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  model_ran <- model
  model_fix <- fit_counterpart_model(model, target)

  # Get genotype variance
  vc_g <- geno_components_from_lme4(model, target, calc_C22 = FALSE)$vc_g
  vc_g <- mean(Matrix::diag(vc_g))

 # Get mean variance of a difference between genotypes
  d_BLUE <- emmeans::emmeans(model_fix, specs = as.formula(paste("pairwise ~", target)))$contrasts |> as.data.frame()
  d_BLUE$var <- d_BLUE$SE^2 # Get variance

  vd_BLUE_avg <- mean(d_BLUE$var)

  # vc_g / (vc_g + vdBLUE_avg / 2)
  H2_Piepho <- H2_Piepho_parameters(vc_g, vd_BLUE_avg)

  return(H2_Piepho)
}

#' @noRd
#' @export
H2_Delta_pairwise.lmerMod <- function(model, target = NULL, type = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is random or fixed
  if(type == "BLUE") {
    H2_Delta <- H2_Delta_BLUE_pairwise.lmerMod(model, target, options)
  } else if(type == "BLUP") {
    H2_Delta <- H2_Delta_BLUP_pairwise.lmerMod(model, target, options)
  }

  return(H2_Delta)
}

#' @keywords internal
H2_Delta_BLUE_pairwise.lmerMod <- function(
  model,
  target = NULL,
  options = NULL
) {
  initial_checks(model, target, options)

  model_fix <- fit_counterpart_model(model, target)
  model_ran <- model

  # Extract vc_g and vc_e
  vc_g <- geno_components_from_lme4(model, target, calc_C22 = FALSE)$vc_g
  vc_g <- mean(Matrix::diag(vc_g))

  # Calculate mean variance of a difference between genotypes
  deltas <- emmeans::emmeans(
    model_fix,
    specs = as.formula(paste("pairwise ~", target))
  )$contrasts |>
    as.data.frame()
  deltas$var <- deltas$SE^2 # Get variance

  # Take pairwise differences and turn into variance-covariance matrix
  lev_g <- levels(model_fix@frame[[target]])
  n_g <- length(lev_g)

  # Create variance-covariance matrix for genotypes (H2: covariance = 0)
  # TODO: For narrow sense, we will need to replace this from the kinship matrix
  cov_g <- matrix(0, nrow = n_g, ncol = n_g)
  diag(cov_g) <- vc_g # Set diagonal to genotype variance
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

  # H2 Delta BLUE
  H2_Delta_BLUE <- H2_Delta_BLUE_parameters(vc_g, Vd_g)

  return(H2_Delta_BLUE)
}

#' @keywords internal
H2_Delta_BLUP_pairwise.lmerMod <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  g <- geno_components_from_lme4(model, target)
  vc_g <- mean(Matrix::diag(g$vc_g))
  C22_g <- g$C22_g

  # Compute variance of difference from PEV
  Vd_g <- outer(
    1:nrow(C22_g), 1:ncol(C22_g),
    Vectorize(function(i, j) C22_g[i, i] + C22_g[j, j] - 2 * C22_g[i, j])
  )

  diag(Vd_g) <- NA
  dimnames(Vd_g) <- list(g$gnames, g$gnames)


  # H2 Delta BLUP
  H2_Delta_BLUP <- H2_Delta_BLUP_parameters(vc_g, Vd_g)

  row.names(H2_Delta_BLUP) <- rownames(Vd_g)
  colnames(H2_Delta_BLUP) <- colnames(Vd_g)

  H2_Delta_BLUP
}




