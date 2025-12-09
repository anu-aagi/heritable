#' Calculate narrow-sense heritability
#'
#' @param model A fitted model object. Currently only supports models with class `asreml`
#' @param target The name of the random effect for which heritability is to be calculated.
#' @param method The method to use for calculating heritability. Options are "Cullis", "Oakey", "BLUE", "BLUP", "Piepho", "Reg", and "SumDiv". Default is "Cullis".
#' @param ... Additional arguments passed to specific methods. See Details
#' @export
h2 <- function(model, ...) {
  UseMethod("h2")
}

#' @export
h2.default <- function(
    model,
    target = NULL,
    method = c("Oakey", "Delta"),
    ...) {
  method <- match.arg(method, several.ok = TRUE)

  initial_checks(model, target, options = NULL)

  h2_values <- sapply(method, function(m) {
    switch(m,
      # Cullis = h2_Cullis(model, target, options = list(check = FALSE)),
      Oakey = h2_Oakey(model, target, options = list(check = FALSE)),
      # Piepho = h2_Piepho(model, target, options = list(check = FALSE)),
      # Delta = h2_Delta(model, target, options = list(check = FALSE)),
      # Standard = h2_Standard(model, target, options = list(check = FALSE)),
      cli::cli_abort(
        "{.fn h2} is not implemented for method {.value m} of class{?es} {.code {class(model)}}"
      )
    )
  })

  # Set names and class
  h2_values <- stats::setNames(h2_values, method)
  structure(h2_values,
    class = c("heritable", class(h2_values)),
    model = model, target = target
  )
}


#' @export
h2_Oakey <- function(model, target, options) {
  UseMethod("h2_Oakey")
}

#' @export
h2_Delta <- function(model, target = NULL, type = c("BLUE", "BLUP"), aggregate, options) {
  UseMethod("h2_Delta")
}

#' @export
h2_Delta.default <- function(model,
                             target = NULL,
                             type = c("BLUP", "BLUE"),
                             aggregate = c("arithmetic", "harmonic"),
                             options = NULL) {
  aggregate <- match.arg(aggregate)
  type <- match.arg(type)

  H2D_ij <- h2_Delta_pairwise(model, target, type = type)
  delta_values <- H2D_ij[upper.tri(H2D_ij)]

  switch(aggregate,
    "arithmetic" = mean(delta_values),
    "harmonic" = length(delta_values) / sum(1 / delta_values)
  )
}


#' @export
h2_Delta_pairwise <- function(model, target, type = c("BLUE", "BLUP"), options) {
  UseMethod("h2_Delta_pairwise")
}



#' Calculate broad-sense heritability
#' @param method Character vector of methods to calculate heritability.
#'        Options are "Cullis", "Oakey", "Delta", "Piepho", and "Standard".
#' @param target The name of the random effect for which heritability is to be calculated.
#' @details The following heritability methods are currently implemented:
#'
#' - Cullis: \deqn{H^2_{Cullis} = 1 - \frac{PEV^{BLUP}_{\overline\Delta ij}}{2\sigma^2_g}}
#' - Oakey: \deqn{H^2_{Oakey} = \frac{\sum_{i = n_z+1}^{n_g} \lambda_i}{\sum_{n_g}^{\lambda_i\neq 0}}}
#' @references
#' - Cullis, B. R., Smith, A. B., & Coombes, N. E. (2006). On the design of early generation variety trials with correlated data. Journal of Agricultural, Biological, and Environmental Statistics, 11(4), 381–393. https://doi.org/10.1198/108571106X154443
#' Oakey, H., Verbyla, A., Pitchford, W., Cullis, B., & Kuchel, H. (2006). Joint modeling of additive and non-additive genetic line effects in single field trials. Theoretical and Applied Genetics, 113(5), 809–819. https://doi.org/10.1007/s00122-006-0333-z
#' @export
H2 <- function(model, target = NULL,
               method = c("Cullis", "Oakey", "Delta", "Piepho", "Standard")
               ) {
  UseMethod("H2")
}

#' @importFrom stats setNames
#' @export
H2.default <- function(model, target = NULL, method = c("Cullis", "Oakey", "Piepho", "Delta", "Standard")) {
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
      cli::cli_abort("{.fn H2} is not implemented for method {.val {m}} of class{?es} {.code {class(model)}}")
    )
  })

  # Set names and class
  H2_values <- stats::setNames(H2_values, method)
  structure(H2_values,
    class = c("heritable", class(H2_values)),
    model = model, target = target
  )
}

#' Calculate broad-sense heritability using Cullis method
#'
#'
#' @export
H2_Cullis <- function(model, target = NULL, options) {
  UseMethod("H2_Cullis")
}

#' @export
H2_Oakey <- function(model, target = NULL, options) {
  UseMethod("H2_Oakey")
}

#' @export
H2_Piepho <- function(model, target = NULL, options) {
  UseMethod("H2_Piepho")
}

#' @export
H2_Delta <- function(
    model,
    target = NULL,
    type = c("BLUP", "BLUE"),
    aggregate = c("arithmetic", "harmonic"),
    options) {
  UseMethod("H2_Delta")
}

#' @export
H2_Delta.default <- function(model,
                             target = NULL,
                             type = c("BLUP", "BLUE"),
                             aggregate = c("arithmetic", "harmonic"),
                             options = NULL) {
  aggregate <- match.arg(aggregate)
  type <- match.arg(type)

  H2D_ij <- H2_Delta_pairwise(model, target, type = type)
  delta_values <- H2D_ij[upper.tri(H2D_ij)]

  switch(aggregate,
    "arithmetic" = mean(delta_values),
    "harmonic" = length(delta_values) / sum(1 / delta_values)
  )
}

#' @export
H2_Delta_by_genotype <- function(model, target, type, options) {
  UseMethod("H2_Delta_by_genotype")
}

#' @export
H2_Delta_by_genotype.default <- function(model, target = NULL, type = NULL, options = NULL) {
  H2D_ij <- H2_Delta_pairwise(model, target, type, options)

  H2D_i <- as.matrix(H2D_ij) |>
    rowMeans(na.rm = TRUE) |>
    data.frame()

  H2D_i <- setNames(H2D_i, "H2D_i")

  H2D_i_list <- split(H2D_i, rownames(H2D_i))

  return(H2D_i_list)
}


#' @export
H2_Delta_pairwise <- function(model, target, type, options) {
  UseMethod("H2_Delta_pairwise")
}

#' @export
H2_Standard <- function(model, target, options) {
  UseMethod("H2_Standard")
}
