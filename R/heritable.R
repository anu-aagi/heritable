
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
H2_cullis <- function(model, ...) {
  UseMethod("H2_cullis")
}

#' @export
H2.default <- function(model, ...) {
  cli::cli_abort("{.fn H2} is not implemented for class {class(model)}")
}



