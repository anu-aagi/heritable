H2_Delta.asreml <- function(model, target = NULL, mean = c("arithmetic", "harmonic")) {
  mean <- match.arg(mean)
  # If model has not converged, warn
  check_model_convergence(model)

  # Check if target is in model, if not throw error
  check_target_exists(model, target)

  # Check if target is random or fixed
  # TODO: How to handle if both fixed and random?
  if(check_target_both(model, target)) {
    cli::cli_abort("The target {.var {target}} is fitted as both fixed and random effect")
  }
  # If fixed, compute H2 delta with BLUEs
  if(!check_target_random(model, target)) { 
    # Obtain BLUES
    g_lsm <- asreml::predict.asreml(model, classify = target, sed = TRUE) 
    Vd_g <- g_lsm$sed^2 # cov

    # Fit counterpart model with target as random for vc_g
    model <- fit_counterpart_model.asreml(model, target)
  } else if(check_target_random(model, target)) {  # If random, compute H2 delta with BLUPs
    # Obtain BLUPS
    g_pred <- asreml::predict.asreml(model,
    classify = target,
    only = target,
    sed = TRUE,
    vcov = TRUE
    )

    Vd_g <- g_pred$sed^2 # cov
  }

  genotype_names <- levels(model$mf[[target]]) # list of genotype names
  ngeno <- length(genotype_names) # number of genotypes
  dimnames(Vd_g) <- list(genotype_names, genotype_names) # name the covariance matrix

  vc_g <- asreml::summary.asreml(model)$varcomp[target, "component"] # varcomp of geno
  # TODO: Change below to use elements of Kinship/relationship as necessary for narrowsense
  var1 <- vc_g
  var2 <- vc_g
  cov <- 0
  v <- var1 + var2 - 2 * cov
  
  # For Delta BLUES
  # Check if g_lsm exists
  if(exists("g_lsm")){
    H2D_ij <- 1 / (1 + (Vd_g / v))
  } else {
    # For Delta BLUPS
    H2D_ij <- 1 - Vd_g / v
  }
  
  if(mean == "arithmetic") {
    H2D_ij <- mean(H2D_ij[upper.tri(H2D_ij)], na.rm = TRUE)
  } else if (mean == "harmonic") {
    H2D_ij <- length(H2D_ij[upper.tri(H2D_ij)]) / sum(1 / H2D_ij[upper.tri(H2D_ij)], na.rm = TRUE)
  }

  return(setNames(H2D_ij, "H2_Delta"))
}

H2_Delta_by_genotype.asreml <- function(model, target = NULL) {

}

H2_Delta_pairwise.asreml <- function(model, target = NULL) {

}





H2_delta_blup <- function(model, target, by = "all") {
  # TODO Check if target exists in any terms
  if (!target %in% c(pull_terms(model)$fixed, pull_terms(model)$random)) {
    cli::cli_abort("{.var {target}} is not fitted in the model. Cannot calculate H2 Delta.")
  }
  # TODO If in fixed, fit counter model
  if (target %in% pull_terms(model)$fixed) {
    model_ran <- fit_counterpart_model.asreml(model, target)
  } else
  # TODO If in random, proceed
  if (target %in% pull_terms(model)$random) {
    model_ran <- model
  }

  # BLUPs for genotype main effect
  g_pred <- asreml::predict.asreml(model_ran,
    classify = target,
    only = target,
    sed = TRUE,
    vcov = TRUE
  )

  genotype_names <- levels(model_ran$mf[[target]]) # list of genotype names
  ngeno <- length(genotype_names) # number of genotypes

  vc_g <- asreml::summary.asreml(model_ran)$varcomp[target, "component"] # varcomp of geno
  # TODO: Change below to use elements of Kinship/relationship as necessary for narrowsense
  var1 <- vc_g
  var2 <- vc_g
  cov <- 0
  Vd_g <- g_pred$sed^2 # cov
  dimnames(Vd_g) <- list(genotype_names, genotype_names)

  v <- var1 + var2 - 2 * cov
  H2D_ij <- 1 - Vd_g / v

  if (by == "all") {
    H2D_ij[upper.tri(H2D_ij)] |> mean()
  } else if (by == "genotype") {
    # Take the mean of get row of matrix
    H2D_i <- as.matrix(H2D_ij) |>
      rowMeans(na.rm = TRUE) |>
      data.frame()

    setNames(H2D_i, "H2D_i")
  } else if (by == "pairwise") {
    H2D_ij
  }
}



H2_delta_blue <- function(model, target, by = "all") {
  # TODO Check if target exists in any terms
  if (!target %in% c(pull_terms(model)$fixed, pull_terms(model)$random)) {
    cli::cli_abort("{.var {target}} is not fitted in the model. Cannot calculate H2 Delta.")
  }
  # TODO If in random, fit counter model
  if (target %in% pull_terms(model)$random) {
    model_ran <- model
    model_fix <- fit_counterpart_model.asreml(model, target)
  } else
  # TODO If in fixed, fit counter model
  if (target %in% pull_terms(model)$fixed) {
    model_fix <- model
    model_ran <- fit_counterpart_model.asreml(model, target)
  }
  # TODO If in both fixed and random, proceed?

  genotype_names <- levels(model$mf[[target]]) # list of genotype names
  ngeno <- length(genotype_names) # number of genotypes

  # Fit genotype as fixed effect model to get BLUEs
  g_lsm <- asreml::predict.asreml(model_fix, classify = target, sed = TRUE)

  vc_g <- asreml::summary.asreml(model_ran)$varcomp[target, "component"] # varcomp of geno
  # TODO: Change below to use elements of Kinship/relationship as necessary for narrowsense
  var1 <- vc_g
  var2 <- vc_g
  cov <- 0
  Vd_g <- g_lsm$sed^2 # cov
  dimnames(Vd) <- list(genotype_names, genotype_names)

  # H2D_ij = Numerator / Denominator
  # var1 + var2 - 2*cov / var1 + var2 - 2*cov + Vd
  v <- var1 + var2 - 2 * cov
  H2D_ij <- 1 - Vd_g / v

  if (by == "all") {
    H2D_ij[upper.tri(H2D_ij)] |> mean()
  } else if (by == "genotype") {
    # TODO: Take the mean of get row of matrix
    upper.tri(H2D_ij)
  } else if (by == "pairwise") {
    H2D_ij
  }
}
