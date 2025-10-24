#' @export
H2.asreml <- function(model, target = NULL, method = c("Cullis", "Oakey", "BLUE", "BLUP", "Piepho", "Reg", "SumDiv")) {
  method <- match.arg(method)
  # If model has not converged, warn
  if (!model$converge) cli::cli_warn("The input model has not converged")

  # TODO: Check if target is in model, if not throw error

  # TODO Not all output is suppressed, even when I added surppressMessages()

  H2 <- switch(method,
    Cullis = H2_Cullis.asreml(model, target),
    Oakey = H2_Oakey.asreml(model, target),
    Piepho = H2_Piepho.asreml(model, target),
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

# TODO: How to best handle multiple model objects in this workflow?
# What happens if user has target in both fixed and random
# If G has an interaction with another parameter, needs to be isolated
#' @export

H2_Piepho.asreml <- function(model, target = NULL) {
  # Obtain the requested random effect
  vc_g <- asreml::summary.asreml(model)$varcomp[target, "component"]

  # Calculate the mean variance of a difference of two genotypic BLUEs
  model_fix <- fit_counterpart_model.asreml(model, target)

  vdBLUE.mat <- asreml::predict.asreml(model_fix,
    classify = target,
    sed = TRUE
  )$sed^2

  vdBLUE.avg <- mean(vdBLUE.mat[upper.tri(vdBLUE.mat, diag = FALSE)])

  # Calculate Piepho's H2
  vc_g / (vc_g + (vdBLUE.avg / 2))
}

#' @export
H2_Reg.asreml <- function(model, target = NULL) {
  # Obtain BLUES

}

#' @export
H2_delta_BLUPS.asreml <- function(model, target = NULL) {
  # Check if target is in random effects
  if(! target %in% pull_terms(model)$random){
    cli::cli_abort("{.var {target}} is not fitted as a random effect in the model. Cannot calculate H2 based on delta BLUPs.")
  }

  # Get BLUPS

  
}