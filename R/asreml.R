

#' @export
H2.asreml <- function(model, target = NULL, method = c("cullis", "oakey")) {
  method <- match.arg(method)
  # if model has not converged, warn
  if(!model$converge) cli::cli_warn("The input model has not converged")
  H2 <- switch(method,
               cullis = H2_cullis.asreml(model, target),
               H2.default())
  structure(H2, class = c("heritable", class(H2)))
}


#' @export
H2_cullis.asreml <- function(model, target = NULL) {
  vc_g <- asreml::summary.asreml(model)$varcomp[target, 'component']
  # need to supress output from below
  vdBLUP_mat <- asreml::predict.asreml(model,
                                       classify = target,
                                       only = target,
                                       sed = TRUE)$sed^2
  vdBLUP_avg <- mean(vdBLUP_mat[upper.tri(vdBLUP_mat, diag = FALSE)])

  1 - (vdBLUP_avg / 2 / vc_g)
}
