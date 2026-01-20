#' @noRd
#' @export
H2_Standard.lmerMod <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  if(options$check %||% TRUE){
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense")
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  # Get genotype variance
  G_g <- var_comp(model, target, calc_C22 = FALSE)$G_g
  s2_g <- mean(diag(G_g))

  # Get residual variance
  s2_eps <- stats::sigma(model)^2

  n_r <- table(model@flist[[target]])

  H2_Standard <- H2_Standard_parameters(s2_g, s2_eps, n_r)

  return(H2_Standard)
}

#' @noRd
#' @export
H2_Cullis.lmerMod <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  if(options$check %||% TRUE){
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense")
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  g <- var_comp(model, target)
  s2_g <- mean(diag(g$G_g))
  n <- g$n_g
  C22_g <- g$C22_g

  # This is equivalent to delta <- var_diff(C22_g); delta_avg = mean(delta[lower.tri(delta)])
  delta_avg <-  (2 / (n * (n - 1))) * (n * sum(diag(C22_g)) - sum(C22_g))

  return(H2_Cullis_parameters(delta_avg, s2_g))
}

#' @noRd
#' @export
H2_Oakey.lmerMod <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  if(options$check %||% TRUE){
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense")
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }
  g <- var_comp(model, target)
  G_g_inv <- Matrix::chol2inv(chol(g$G_g))

  return(H2_Oakey_parameters(G_g_inv, g$C22_g))
}

#' @noRd
#' @export
H2_Piepho.lmerMod <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  if(options$check %||% TRUE){
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense")
  }

  conterpart <- fit_counterpart_model(model, target)

  # Get genotype variance
  G_g <- var_comp(model, target, calc_C22 = FALSE)$G_g
  s2_g <- mean(diag(G_g))

  # Get mean variance of a difference between genotypes
  frm <-  as.formula(paste("pairwise ~", target))
  EMM_fit <- emmeans::emmeans(conterpart, specs = frm)$contrasts
  delta_avg <- mean(data.frame(EMM_fit)$SE^2)  # Get variance

  # s2_g / (s2_g + delta_avg / 2)
  H2_Piepho <- H2_Piepho_parameters(s2_g, delta_avg)

  return(H2_Piepho)
}

#' @noRd
#' @export
H2_Delta_pairwise.lmerMod <- function(model, target = NULL, type = NULL, options = NULL) {

  initial_checks(model, target, options)

  if(options$check %||% TRUE){
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense")
  }

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

  conterpart <- fit_counterpart_model(model, target)

  # Extract vc_g and vc_e
  g <- var_comp(model, target, calc_C22 = FALSE)
  s2_g <- mean(diag(g$G_g))

  # Calculate mean variance of a difference between genotypes
  frm <-  as.formula(paste("pairwise ~", target))
  EMM_fit <- emmeans::emmeans(conterpart, specs = frm)$contrasts
  EMM_fit <- data.frame(EMM_fit)  # Get variance
  EMM_fit$var <- EMM_fit$SE^2

  # Take pairwise differences and turn into variance-covariance matrix
  gnames <- g$gnames
  n_g <- g$n_g

  # Start with empty variance matrix for differences
  delta <- matrix(0, nrow = n_g, ncol = n_g)
  dimnames(delta) <- list(gnames, gnames)

  # Fill in the pairwise variance of differences
  for (i in 1:nrow(EMM_fit)) {
    # Extract genotype names from contrast column
    pair <- strsplit(as.character(EMM_fit$contrast[i]), " - ")[[1]]
    g1 <- pair[1]
    g2 <- pair[2]

    delta[g1, g2] <- EMM_fit$var[i]
    delta[g2, g1] <- EMM_fit$var[i] # symmetric
  }

  # H2 Delta BLUE
  H2_Delta_BLUE <- H2_Delta_BLUE_parameters(s2_g, delta)

  return(H2_Delta_BLUE)
}

#' @keywords internal
H2_Delta_BLUP_pairwise.lmerMod <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  g <- var_comp(model, target)
  s2_g <- mean(diag(g$G_g))
  C22_g <- g$C22_g

  # Compute variance of difference from PEV
  delta <- var_diff(C22_g)
  diag(delta) <- NA
  dimnames(delta) <- list(g$gnames, g$gnames)

  # H2 Delta BLUP
  H2_Delta_BLUP <- H2_Delta_BLUP_parameters(s2_g, delta)

  dimnames(H2_Delta_BLUP) <- dimnames(delta)

  H2_Delta_BLUP
}




