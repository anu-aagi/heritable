#' Class coercion with heritable objects
#' @importFrom vctrs vec_ptype2
#' @noRd
#' @export
#' @examples
#' # Create a simple heritable object for testing
#' h1 <- structure(c(0.8, 0.75, 0.82),
#'                 class = c("heritable", "numeric"),
#'                 names = c("Cullis", "Oakey", "Delta"))
#'
#' h2 <- structure(c(0.85, 0.78, 0.80),
#'                 class = c("heritable", "numeric"),
#'                 names = c("Cullis", "Oakey", "Delta"))
#'
#' vec_ptype2(h1, h2)
#' vec_cast(h1, h2)
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
