

#' @export
H2_Oakey.asreml <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  n_g <- model$noeff[[target]]
  vc_g <- summary(model)$varcomp[target, "component"]
  vcov_g <- asreml::predict.asreml(model, classify = target, only = target, vcov = TRUE)$vcov

  H2_Oakey <- H2_Oakey_parameters(n_g, vc_g, vcov_g)

  return(H2_Oakey)
}

#' @export
H2_Cullis.asreml <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  # Get genotype variance
  vc_g <- asreml::summary.asreml(model)$varcomp[target, "component"]

  vdBLUP_mat <- asreml::predict.asreml(model,
    classify = target,
    only = target,
    sed = TRUE
  )$sed^2

  vd_BLUP_avg <- mean(vdBLUP_mat[upper.tri(vdBLUP_mat, diag = FALSE)])

  H2_Cullis <- H2_Cullis_parameters(vd_BLUP_avg, vc_g)

  return(H2_Cullis)
}

#' @export
H2_Piepho.asreml <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is random or fixed
  # TODO: 3 way check here, if both do nothing, if random fit fixed, if fixed fit random
  if (check_target_random(model, target)) {
    model_fix <- fit_counterpart_model.asreml(model, target)
    model_ran <- model
  } else {
    model_fix <- model
    model_ran <- fit_counterpart_model.asreml(model, target)
  }

  # Calculate the mean variance of a difference of two genotypic BLUEs
  # Get genotype variance
  vc_g <- asreml::summary.asreml(model_ran)$varcomp[target, "component"]

  vdBLUE_mat <- asreml::predict.asreml(model_fix,
    classify = target,
    only = target,
    sed = TRUE
  )$sed^2

  vdBLUE_avg <- mean(vdBLUE_mat[upper.tri(vdBLUE_mat, diag = FALSE)])

  # Calculate Piepho's H2
  H2_Piepho <- H2_Piepho_parameters(vc_g, vdBLUE_avg)

  return(H2_Piepho)
}

#' @export
H2_Delta_pairwise.asreml <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is random or fixed
  # TODO: How to handle if both fixed and random?
  if(check_target_both(model, target)) {
    cli::cli_abort("The target {.var {target}} is fitted as both fixed and random effect")
  }
  gpred <- asreml::predict.asreml(model, classify = target, only = target, sed = TRUE)
  Vd_g <- gpred$sed^2  # Variance of difference

  genotype_names <- gpred$pvals[[target]] # list of genotype names
  ngeno <- length(genotype_names) # number of genotypes
  dimnames(Vd_g) <- list(genotype_names, genotype_names) # name the covariance matrix

  # If fixed, compute H2 delta with BLUEs
  if(!check_target_random(model, target)) {
    # Fit counterpart model with target as random for vc_g
    model <- fit_counterpart_model.asreml(model, target)
    vc_g <- model$vparameters[[target]] * model$sigma2 # varcomp of geno
    H2_Delta_BLUE_parameters(vc_g, vc_g, cov = 0, Vd_g)
  } else {
    vc_g <- model$vparameters[[target]] * model$sigma2
    H2_Delta_BLUP_parameters(vc_g, vc_g, cov = 0, Vd_g)
  }
}

#' @export
H2_Naive.asreml <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }
  vc_g <- asreml::summary.asreml(model)$varcomp[target, "component"]
  vc_e <- asreml::summary.asreml(model)$varcomp["units!R", "component"]
  # TODO: may need to remove observations where phenotype is NA
  n_r <- table(model$mf[[target]])

  H2_Naive <- H2_Naive_parameters(vc_g, vc_e, n_r)

  return(H2_Naive)
}
