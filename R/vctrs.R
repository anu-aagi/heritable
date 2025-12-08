#' @importFrom vctrs vec_ptype2
#' @export
vec_ptype2.heritable.heritable <- function(x, y, ...) double()

#' @export
vec_ptype2.heritable.double <- function(x, y, ...) double()

#' @export
vec_ptype2.double.heritable <- function(x, y, ...) double()

#' @importFrom vctrs vec_cast
#' @export
vec_cast.heritable.heritable <- function(x, to, ...) x

#' @export
vec_cast.double.heritable <- function(x, to, ...) x

#' @export
vec_cast.heritable.double <- function(x, to, ...) x
