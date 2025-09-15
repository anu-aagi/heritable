

#' @export
H2.asreml <- function(model, target = NULL, method = c("Cullis", "Oakey", "BLUE", "BLUP", "Piepho", "Reg", "SumDiv")) {
  method <- match.arg(method)
  # if model has not converged, warn
  if(!model$converge) cli::cli_warn("The input model has not converged")
  H2 <- switch(method,
               Cullis = H2_Cullis.asreml(model, target),
               Oakey = H2_Oakey(model, target),
               H2.default())
  structure(H2, class = c("heritable", class(H2)))
}

#' @export
H2_Oakey.asreml <- function(model, target = NULL) {
  n_g   <- model$noeff[[target]]
  vc_g  <- summary(model)$varcomp[target, 'component']
  Gg_inv   <- diag(1/vc_g, nrow = n_g, ncol = n_g)
  C22_g <- asreml::predict.asreml(model, classify = target, only = target, vcov=TRUE)$vcov
  M     <- diag(n_g) - (Gg_inv %*% C22_g)
  eM    <- eigen(M)

  sum(eM$values) / (n_g - 1)
}

#' @export
H2_Cullis.asreml <- function(model, target = NULL) {
  vc_g <- asreml::summary.asreml(model)$varcomp[target, 'component']
  # need to supress output from below
  vdBLUP_mat <- asreml::predict.asreml(model,
                                       classify = target,
                                       only = target,
                                       sed = TRUE)$sed^2
  vdBLUP_avg <- mean(vdBLUP_mat[upper.tri(vdBLUP_mat, diag = FALSE)])

  1 - (vdBLUP_avg / 2 / vc_g)
}
