# h2_Cullis.asreml <- function(model, target = NULL, options = NULL) {
#
#   initial_checks(model, target, options)
#
#   # Check if target is random or fixed
#   if (!check_target_random(model, target)) {
#     return(NA)
#   }
#
#   vm <- target_vm_term_asreml(model, target)
#
#   n_g <-  model$noeff[[vm$target_vm]]
#   vc_g <- model$vparameters[[vm$target_vm]] * model$sigma2 * semivariance(vm$GRM)
#
#
#   vdBLUP_mat <- asreml::predict.asreml(model,
#                                        classify = target,
#                                        only = target,
#                                        sed = TRUE,
#                                        trace = FALSE,
#   )$sed^2
#
#   vd_BLUP_avg <- mean(vdBLUP_mat[upper.tri(vdBLUP_mat, diag = FALSE)])
#
#   H2_Cullis_parameters(vd_BLUP_avg, vc_g)
# }

#' @export
h2_Oakey.asreml <- function(model, target = NULL, options = NULL) {
  initial_checks(model, target, options)

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  vm <- target_vm_term_asreml(model, target)
  n_g <- model$noeff[[vm$target_vm]]
  Gg_inv <- 1 / (model$vparameters[[vm$target_vm]] * model$sigma2) * vm$GRMinv
  vcov_g <- asreml::predict.asreml(model,
    classify = target,
    only = target,
    vcov = TRUE,
    trace = FALSE
  )$vcov


  H2_Oakey_parameters(Gg_inv, vcov_g)
}


#' @export
h2_Delta_pairwise.asreml <- function(model, target = NULL, type = NULL, options = NULL) {
  initial_checks(model, target, options)

  vm <- target_vm_term_asreml(model, target)
  n_g <- model$noeff[[vm$target_vm]]
  Gg <- model$vparameters[[vm$target_vm]] * model$sigma2 * solve(vm$GRMinv)

  if (type == "BLUP") {
    gpred <- asreml::predict.asreml(model, classify = target, sed = TRUE, trace = FALSE)
    Vd_g <- gpred$sed^2 # Variance of difference
    genotype_names <- gpred$pvals[[target]] # list of genotype names
    dimnames(Vd_g) <- list(genotype_names, genotype_names) # name the covariance matrix
    h2_Delta_BLUP_parameters(Gg, Vd_g)
  } else if (type == "BLUE") {
    model_fix <- fit_counterpart_model.asreml(model, target)
    gpred <- asreml::predict.asreml(model_fix, classify = target, sed = TRUE, trace = FALSE)
    Vd_g <- gpred$sed^2 # Variance of difference
    genotype_names <- gpred$pvals[[target]] # list of genotype names
    dimnames(Vd_g) <- list(genotype_names, genotype_names) # name the covariance matrix
    h2_Delta_BLUE_parameters(Gg, Vd_g)
  }
}

get_vc_g_asreml <- function(model, target) {
  model$vparameters[[target]] * model$sigma2
}


#' @export
H2_Cullis.asreml <- function(model, target = NULL, options = NULL) {
  initial_checks(model, target, options)

  # Get genotype variance
  vc_g <- get_vc_g_asreml(model, target)

  vdBLUP_mat <- asreml::predict.asreml(model,
    classify = target,
    only = target,
    sed = TRUE,
    trace = FALSE
  )$sed^2

  vd_BLUP_avg <- mean(vdBLUP_mat[upper.tri(vdBLUP_mat, diag = FALSE)])

  H2_Cullis <- H2_Cullis_parameters(vd_BLUP_avg, vc_g)

  return(H2_Cullis)
}

#' @export
H2_Oakey.asreml <- function(model, target = NULL, options = NULL) {
  initial_checks(model, target, options)

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  n_g <- model$noeff[[target]]
  vc_g <- get_vc_g_asreml(model, target)
  vcov_g <- asreml::predict.asreml(model,
    classify = target,
    only = target,
    vcov = TRUE,
    trace = FALSE
  )$vcov

  Gg_inv <- diag(1 / vc_g, nrow = n_g, ncol = n_g)

  H2_Oakey <- H2_Oakey_parameters(Gg_inv, vcov_g)

  return(H2_Oakey)
}


#' @export
H2_Piepho.asreml <- function(model, target = NULL, options = NULL) {
  initial_checks(model, target, options)

  model_fix <- fit_counterpart_model.asreml(model, target)
  model_ran <- model

  # Calculate the mean variance of a difference of two genotypic BLUEs
  # Get genotype variance
  vc_g <- get_vc_g_asreml(model_ran, target)

  vdBLUE_mat <- asreml::predict.asreml(model_fix,
    classify = target,
    sed = TRUE,
    trace = FALSE
  )$sed^2

  vdBLUE_avg <- mean(vdBLUE_mat[upper.tri(vdBLUE_mat, diag = FALSE)])

  # Calculate Piepho's H2
  H2_Piepho <- H2_Piepho_parameters(vc_g, vdBLUE_avg)

  return(H2_Piepho)
}

#' @export
H2_Delta_pairwise.asreml <- function(model, target = NULL, type = NULL, options = NULL) {
  initial_checks(model, target, options)
  vc_g <- get_vc_g_asreml(model, target)
  if (type == "BLUP") {
    gpred <- asreml::predict.asreml(model, classify = target, sed = TRUE, trace = FALSE)
    Vd_g <- gpred$sed^2 # Variance of difference
    genotype_names <- gpred$pvals[[target]] # list of genotype names
    dimnames(Vd_g) <- list(genotype_names, genotype_names) # name the covariance matrix
    H2_Delta_BLUP_parameters(vc_g, Vd_g)
  } else if (type == "BLUE") {
    model_fix <- fit_counterpart_model.asreml(model, target)
    gpred <- asreml::predict.asreml(model_fix, classify = target, sed = TRUE, trace = FALSE)
    Vd_g <- gpred$sed^2 # Variance of difference
    genotype_names <- gpred$pvals[[target]] # list of genotype names
    dimnames(Vd_g) <- list(genotype_names, genotype_names) # name the covariance matrix
    H2_Delta_BLUE_parameters(vc_g, Vd_g)
  }
}

#' @export
H2_Standard.asreml <- function(model, target = NULL, options = NULL) {
  initial_checks(model, target, options)

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }
  vc_g <- get_vc_g_asreml(model, target)
  vc_e <- model$sigma2
  # TODO: may need to remove observations where phenotype is NA
  n_r <- table(model$mf[[target]])

  H2_Standard <- H2_Standard_parameters(vc_g, vc_e, n_r)

  return(H2_Standard)
}
