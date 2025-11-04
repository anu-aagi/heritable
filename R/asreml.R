#' @export
H2.asreml <- function(model, target = NULL, method = c("Cullis", "Oakey", "Delta", "BLUE", "BLUP", "Piepho", "Reg", "SumDiv", "Naive")) {
  method <- match.arg(method)

  # If model has not converged, warn
  check_model_convergence(model)

  # TODO: Check if model is of class asreml
  # In the event someone calls the asreml function directly on a different class
  check_model_class_asreml(model)

  # Check if target is in model, if not throw error
  check_target_exists(model, target)

  H2 <- switch(method,
    Cullis = H2_Cullis.asreml(model, target),
    Oakey = H2_Oakey.asreml(model, target),
    Piepho = H2_Piepho.asreml(model, target),
    Delta_blup = H2_Delta.asreml(model, target),
    Naive = H2_Naive.asreml(model, target),
    H2.default(model)
  )

  structure(H2, class = c("heritable", class(H2)))
}

#' @export
H2_Oakey.asreml <- function(model, target = NULL) {
  n_g <- model$noeff[[target]]
  vc_g <- summary(model)$varcomp[target, "component"]
  Gg_inv <- diag(1 / vc_g, nrow = n_g, ncol = n_g)
  C22_g <- asreml::predict.asreml(model, classify = target, only = target, vcov = TRUE)$vcov
  M <- diag(n_g) - (Gg_inv %*% C22_g)
  eM <- eigen(M)

  sum(eM$values) / (n_g - 1)
}

#' @export
H2_Cullis.asreml <- function(model, target = NULL) {
  vc_g <- asreml::summary.asreml(model)$varcomp[target, "component"]

  vdBLUP_mat <- asreml::predict.asreml(model,
    classify = target,
    only = target,
    sed = TRUE
  )$sed^2

  vdBLUP_avg <- mean(vdBLUP_mat[upper.tri(vdBLUP_mat, diag = FALSE)])

  1 - (vdBLUP_avg / 2 / vc_g)
}

#' @export
H2_Piepho.asreml <- function(model, target = NULL) {
  # Obtain the requested random effect
  vc_g <- asreml::summary.asreml(model)$varcomp[target, "component"]

  # Calculate the mean variance of a difference of two genotypic BLUEs
  model_fix <- fit_counterpart_model.asreml(model, target)

  vdBLUE_mat <- asreml::predict.asreml(model_fix,
    classify = target,
    sed = TRUE
  )$sed^2

  vdBLUE_avg <- mean(vdBLUE_mat[upper.tri(vdBLUE_mat, diag = FALSE)])

  # Calculate Piepho's H2
  vc_g / (vc_g + (vdBLUE_avg / 2))
}

#' @export
H2_Reg.asreml <- function(model, target = NULL) {
  # Obtain BLUES
}

#' @export
H2_Delta.asreml <- function(model, target = NULL, type = c("blup", "blue"), by = "all") {
  # TODO Check if target is in fixed or random

  # TODO If neither, throw error
  if (!target %in% pull_terms(model)$random) {
    cli::cli_abort("{.var {target}} is not fitted as a random effect in the model. Cannot calculate H2 based on delta BLUPs.")
  }

  # TODO If in random, proceed

  # TODO If in fixed, fit counter model

  # BLUPs for genotype main effect
  g_pred <- asreml::predict.asreml(model,
    classify = target,
    only = target,
    sed = TRUE,
    vcov = TRUE
  )

  genotype_names <- levels(model$mf[[target]]) # list of genotype names
  ngeno <- length(genotype_names) # number of genotypes


  vc_g <- asreml::summary.asreml(model)$varcomp[target, "component"] # varcomp of geno
  # TODO: Change below to use elements of Kinship/relationship as necessary for narrowsense
  var1 <- vc_g
  var2 <- vc_g
  cov <- 0
  Vd_g <- g_pred$sed^2 # cov
  dimnames(Vd_g) <- list(genotype_names, genotype_names)

  v <- var1 + var2 - 2 * cov
  H2D_ij <- 1 - Vd_g / v
  H2D_ij
}

#' @export
H2_Naive.asreml <- function(model, target = NULL) {
  vc_g <- asreml::summary.asreml(model)$varcomp[target, "component"]
  vc_e <- asreml::summary.asreml(model)$varcomp["units!R", "component"]

  vc_g / (vc_g + vc_e)
}
