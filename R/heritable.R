
#' Calculate narrow-sense heritability
#'
#' @param model The fitted model.
#'
#' @export
h2 <- function(model, ...) {
  UseMethod("h2")
}


#' @export
h2.default <- function(model, ...) {
  cli::cli_abort("{.fn h2} is not implemented for class {class(model)}")
}

#' Calculate broad-sense heritability
#'
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
H2.default <- function(model, ...) {
  cli::cli_abort("{.fn H2} is not implemented for class {class(model)}")
}



