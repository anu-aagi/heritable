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
h2.default <- function(
  model,
  target = NULL,
  method = c("Cullis", "Oakey"),
  ...
) {

  method <- match.arg(method, several.ok = TRUE)

  initial_checks(model, target, options = NULL)

  # Calculate H2 for each method
  h2_values <- sapply(method, function(m) {
    switch(
      m,
      Cullis = h2_Cullis(model, target, options = list(check = FALSE)),
      Oakey = h2_Oakey(model, target, options = list(check = FALSE)),
      #Piepho = h2_Piepho(model, target, options = list(check = FALSE)),
      #Delta = h2_Delta(model, target, options = list(check = FALSE)),
      #Standard = h2_Standard(model, target, options = list(check = FALSE)),
      cli::cli_abort(
        "{.fn h2} is not implemented for method {.value m} of class{?es} {.code {class(model)}}"
      )
    )
  })

  # Set names and class
  h2_values <- stats::setNames(h2_values, method)
  structure(h2_values, class = c("heritable", class(h2_values)))
}

#' @export
h2_Cullis <- function(model, ...) {
  UseMethod("h2_Cullis")
}

#' @export
h2_Oakey <- function(model, ...) {
  UseMethod("h2_Oakey")
}


#' Calculate broad-sense heritability
#' @inheritParams h2
#' @param method Character vector of methods to calculate heritability.
#'        Options are "Cullis", "Oakey", "Delta", "Piepho", and "Standard".
#' @param target The name of the random effect for which heritability is to be calculated.
#' @export
H2 <- function(model, target = NULL,
               method = c("Cullis", "Oakey", "Delta", "Piepho", "Standard"),
               ...) {
  UseMethod("H2")
}

#' @importFrom stats setNames
#' @export
H2.default <- function(model, target = NULL, method = c("Cullis", "Oakey", "Piepho", "Delta", "Standard")) {
  # TODO: This will change if we want to vectorise over multiple methods
  method <- match.arg(method, several.ok = TRUE)

  initial_checks(model, target, options = NULL)

  # Calculate H2 for each method
  H2_values <- sapply(method, function(m) {
    switch(m,
           Cullis = H2_Cullis(model, target, options = list(check = FALSE)),
           Oakey = H2_Oakey(model, target, options = list(check = FALSE)),
           Piepho = H2_Piepho(model, target, options = list(check = FALSE)),
           Delta = H2_Delta(model, target, options = list(check = FALSE)),
           Standard = H2_Standard(model, target, options = list(check = FALSE)),
           cli::cli_abort("{.fn H2} is not implemented for method {.value m} of class{?es} {.code {class(model)}}")
    )
  })

  # Set names and class
  H2_values <- stats::setNames(H2_values, method)
  structure(H2_values, class = c("heritable", class(H2_values)))
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
H2_Piepho <- function(model, ...) {
  UseMethod("H2_Piepho")
}

#' @export
H2_Delta <- function(model, ...) {
  UseMethod("H2_Delta")
}

#' @export
H2_Delta.default <- function(model, target = NULL, aggregate = c("arithmetic", "harmonic"), options = NULL) {
  aggregate <- match.arg(aggregate)

  H2D_ij <- H2_Delta_pairwise(model, target)
  delta_values <- H2D_ij[upper.tri(H2D_ij)]

  switch(aggregate,
         "arithmetic" = mean(delta_values),
         "harmonic" = length(delta_values) / sum(1 / delta_values))
}

#' @export
H2_Delta_by_genotype <- function(model, ...) {
  UseMethod("H2_Delta_by_genotype")
}

#' @export
H2_Delta_by_genotype.default <- function(model, target = NULL, options = NULL) {
  H2D_ij <- H2_Delta_pairwise(model, target)

  H2D_i <- as.matrix(H2D_ij) |>
    rowMeans(na.rm = TRUE) |>
    data.frame()

  H2D_i <- setNames(H2D_i, "H2D_i")

  H2D_i_list <- split(H2D_i, rownames(H2D_i))

  return(H2D_i_list)
}


#' @export
H2_Delta_pairwise <- function(model, ...) {
  UseMethod("H2_Delta_pairwise")
}

#' @export
H2_Standard <- function(model, ...) {
  UseMethod("H2_Standard")
}

