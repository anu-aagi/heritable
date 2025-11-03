#' Calculate narrow-sense heritability
#'
#' @param model A fitted model object. Currently only supports models with class `asreml`
#' @param target The name of the random effect for which heritability is to be calculated.
#' @param method The method to use for calculating heritability. Options are "Cullis", "Oakey", "BLUE", "BLUP", "Piepho", "Reg", and "SumDiv". Default is "Cullis".
#' @param ... Additional arguments passed to specific methods. See Details
#' @details The following heritability methods are currently implemented:
#'
#' - Cullis: [#TODO Insert equation here]
#' - Oakey: [#TODO Insert equation here]
#' @references
#' - Cullis et al. (2006) #TODO
#' - Oakey et al. (2006) #TODO
#' @export
h2 <- function(model, ...) {
  UseMethod("h2")
}

#' @export
h2.default <- function(model, ...) {
  cli::cli_abort("{.fn h2} is not implemented for class{?es} {.code {class(model)}}")
}

#' Calculate broad-sense heritability
#' @inheritParams h2
#'
#' @export
H2 <- function(model, ...) {
  UseMethod("H2")
}

#' Calculate broad-sense heritability using Cullis method
#'
#' @inheritParams h2
#'
#' @export
H2_Cullis <- function(model, ...) {
  UseMethod("H2_Cullis")
}

#' @export
H2_Oakey <- function(model, ...) {
  UseMethod("H2_Oakey")
}

#' @export
H2_BLUE <- function(model, ...) {
  UseMethod("H2_BLUE")
}

#' @export
H2_Piepho <- function(model, ...) {
  UseMethod("H2_Piepho")
}

#' @export
H2_Reg <- function(model, ...) {
  UseMethod("H2_Reg")
}

#' @export
H2_SumDiv <- function(model, ...) {
  UseMethod("H2_SumDiv")
}

#' @export
H2_Delta <- function(model, ...) {
  UseMethod("H2_Delta")
}

#' @export
H2_Naive <- function(model, ...) {
  UseMethod("H2_Naive")
}

#' @export
H2.default <- function(model, ...) {
  cli::cli_abort("{.fn H2} is not implemented for class{?es} {.code {class(model)}}")
}
