#' Fixed-effects-only fitted values from an asreml model
#'
#' @description
#' Returns the fitted values based on fixed effects only
#' (\eqn{\hat{y} = X\hat{\beta}}), excluding all random effects.
#'
#' @param model An object of class \code{"asreml"}, fitted with
#' \code{model.frame = TRUE}.
#' @param source The known genomic relationship matrix (GRM) used in `model` fitted using `asreml::vm()`, provided as a named list.
#' When not provided (an empty list by default), the GRM variable used for `vm` calling will be searched in the global environment.
#'
#' @return
#' A numeric vector of length \eqn{N}, giving the fixed-effects-only
#' fitted value for each observation.
#'
#' @details
#' This function reconstructs the fixed-effect design matrix from the
#' stored model frame and multiplies it by the estimated fixed
#' coefficients. Sparse fixed terms (if any) are included.
#'
#' Random effects (BLUPs) are not included.
#'
#' @export
get_fixed_fit_asreml <- function(model, source = list()) {
  if (!inherits(model, "asreml")) {
    stop("`model` must be an `asreml` object.")
  }

  # Get model frame
  model <- check_design_exists(model, build_mf = FALSE, source = source)
  design <- model$design

  N <- nrow(design)

  # Fixed and sparse term labels
  fixed_terms <- attr(model$formulae$fixed, "term.labels")
  sparse_terms <- attr(model$formulae$sparse, "term.labels")
  has_int <- isTRUE(attr(model$formulae$fixed, "intercept") == 1)

  # Get estimates
  term_names <- c(
    rownames(model$coefficients$fixed),
    rownames(model$coefficients$sparse)
  )
  beta <- c(
    model$coefficients$fixed,
    model$coefficients$sparse
  )

  as.numeric(design[, term_names, drop = FALSE] %*% beta)
}
