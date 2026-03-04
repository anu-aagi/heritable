check_glmmtmb_family <- function(model) {
  if (model$modelInfo$family$family != "gaussian")
    cli::cli_abort("Only the Gaussian response distribution is currently supported!")
}
