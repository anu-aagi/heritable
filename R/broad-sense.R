#' Calculate heritability from mixed model object
#' @description
#' A case-specific wrapper for calculating heritability.
#'
# - The lowercase prefix `h2_` refers to the wrapper or subfunctions e.g. [`h2_Oakey()`] for calculating narrow sense heritability
#' - The upper case prefix `H2_` refers to the wrapper or subfunctions e.g. [`H2_Delta()`] for calculating broad sense heritability
#' @param model Model object of class `lmerMod/merMod` or `asreml`
#' @param method Character vector of name of method to calculate heritability. See details.
#' @param target The name of the random effect for which heritability is to be calculated.
#' @param options NULL by default, for internal checking of model object before calculations
#' @aliases H2
#' @usage
# h2(model, target, method = c("Oakey", "Delta"), options)
#' H2(model, target, method = c("Cullis", "Oakey", "Delta", "Piepho", "Standard"), options)
#' @returns A named numeric vector, length matching number of methods supplied
#' @details
#'
# The following methods are currently implemented for narrow-sense heritability `h2(method = "XX")`:
#
# - `"Oakey"`: \deqn{H^2_{Oakey} = \frac{\sum_{i = n_z+1}^{n_g} \lambda_i}{\sum_{n_g}^{\lambda_i\neq 0}}}
# - `"Delta"`: \deqn{H^2_{\Delta ..} = 1 - \frac{PEV^{BLUP}_{\overline\Delta ..}}{2\sigma^2_g}}
#'
#' The following methods are currently implemented for broad-sense heritability `H2(method = "XX")`:
#'
#' - `"Cullis"`: \deqn{H^2_{Cullis} = 1 - \frac{PEV^{BLUP}_{\overline\Delta ij}}{2\sigma^2_g}}
#' - `"Oakey"`: \deqn{H^2_{Oakey} = \frac{\sum_{i = n_z+1}^{n_g} \lambda_i}{\sum_{n_g}^{\lambda_i\neq 0}}}
#' - `"Delta"`: \deqn{H^2_{\Delta ..} = 1 - \frac{PEV^{BLUP}_{\overline\Delta ..}}{2\sigma^2_g}}
#' - `"Piepho"`: \deqn{H^2_{Piepho} = \frac{\sigma^2_g}{\sigma^2_g + \overline{PEV_{BLUE_g}} / 2}}
#' - `"Standard"`: \deqn{H^2_{Standard} = \frac{\sigma^2_g}{\sigma^2_g + \frac{1}{n_g}\sum_{n_g}^{i=1} \sigma^2_p / n_{gi}}}
#'
#' For further details of a specific method - take a look at help file for each subfunctions `?H2_Cullis`
#'
#' @references
#' - Cullis, B. R., Smith, A. B., & Coombes, N. E. (2006). On the design of early generation variety trials with correlated data. Journal of Agricultural, Biological, and Environmental Statistics, 11(4), 381–393. https://doi.org/10.1198/108571106X154443
#' - Oakey, H., Verbyla, A., Pitchford, W., Cullis, B., & Kuchel, H. (2006). Joint modeling of additive and non-additive genetic line effects in single field trials. Theoretical and Applied Genetics, 113(5), 809–819. https://doi.org/10.1007/s00122-006-0333-z
#' - Schmidt, P., Hartung, J., Rath, J., & Piepho, H.-P. (2019). Estimating Broad-Sense Heritability with Unbalanced Data from Agricultural Cultivar Trials. Crop Science, 59(2), 525–536. https://doi.org/10.2135/cropsci2018.06.0376
#' - Piepho, H.-P., & Möhring, J. (2007). Computing Heritability and Selection Response From Unbalanced Plant Breeding Trials. Genetics, 177(3), 1881–1888. https://doi.org/10.1534/genetics.107.074229
#' - Falconer, D. S., & Mackay, T. F. C. (1996). Introduction to quantitative genetics (4th ed.). Longman.
#' @seealso [H2_Cullis()], [H2_Oakey()], [H2_Delta()], [H2_Piepho()], [H2_Standard()]
# [`h2_Oakey()`], [`h2_Delta()`]
#' @export
H2 <- function(model,
               target = NULL,
               method = c("Cullis", "Oakey", "Delta", "Piepho", "Standard"),
               options = NULL) {
  UseMethod("H2")
}

#' @importFrom stats setNames
#' @noRd
#' @export
H2.default <- function(model,
                       target = NULL,
                       method = c("Cullis", "Oakey", "Piepho", "Delta", "Standard"),
                       options
                       ) {
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

#' Calculate standard heritability from model object
#' @description Compute standard heritability using the classic ratio method of
#' genotypic and phenotypic variance. See Falconer & Mackay (1996)
#' @usage H2_Standard(model, target, options)
#' @inheritParams H2
#' @return Single numeric value
#' @details
#' The equation used to calculate standard heritability is:
#' \deqn{H^2_{Standard} = \frac{\sigma^2_g}{\sigma^2_g + \frac{1}{n_g}\sum_{n_g}^{i=1} \sigma^2_p / n_{gi}}}
#' where:
#'
#' - \eqn{n_g} is the number of genotypes
#' - \eqn{n_{gi}} is the number of replicate for a given genotype i
#' - \eqn{\sigma_g} is the variance attributed to genotype differences
#' - \eqn{\sigma_p} is the variance attributed to phenotypic differences
#' @export
#' @references
#' Falconer, D. S., & Mackay, T. F. C. (1996). Introduction to quantitative genetics (4th ed.). Longman.
H2_Standard <- function(model, target, options) {
  UseMethod("H2_Standard")
}

#' Calculate Cullis' heritability from model object
#' @description Compute "generalised heritability" for unbalanced experimental designs.
#' See Cullis, Smith and Coombes (2006) for derivation.
#' @inheritParams H2
#' @usage H2_Cullis(model, target, options)
#' @return Single numeric value
#' @details The equation for Cullis heritability is as follow
#'
#' \deqn{H^2_{Cullis} = 1 - \frac{PEV^{BLUP}_{\overline\Delta ij}}{2\sigma^2_g}}
#'
#' where:
#' - \eqn{PEV} is the prediction error variance matrix of the pairwise differences among BLUPS
#' - \eqn{\sigma^2} is the variance attributed to differences between genotype
#' @references
#' Cullis, B. R., Smith, A. B., & Coombes, N. E. (2006). On the design of early generation variety trials with correlated data. Journal of Agricultural, Biological, and Environmental Statistics, 11(4), 381–393. https://doi.org/10.1198/108571106X154443
#' @export
H2_Cullis <- function(model, target, options) {
  UseMethod("H2_Cullis")
}

#' @noRd
#' @export
H2_Oakey <- function(model, target, options) {
  UseMethod("H2_Oakey")
}

#' Calculate Piepho's heritability from model object
#' Compute Piepho's heritability using variance differences between genotype BLUEs
#' @usage H2_Piepho(model, target, options)
#' @inheritParams H2
#' @details The equation for Piepho's heritability is as follows:
#'
#' \deqn{H^2_{Piepho} = \frac{\sigma^2_g}{\sigma^2_g + \overline{PEV_{BLUE_g}} / 2}}
#'
#' where:
#' - \eqn{\overline{PEV_{BLUE_g}}} is the prediction error variance matrix for genotype BLUEs
#' - \eqn{\sigma^2_g} is the variance attributed to differences between genotype
#'
#' See reference for full derivation and details.
#' @export
#' @references
#' Piepho, H.-P., & Möhring, J. (2007). Computing Heritability and Selection Response From Unbalanced Plant Breeding Trials. Genetics, 177(3), 1881–1888. https://doi.org/10.1534/genetics.107.074229
#'
H2_Piepho <- function(model, target, options) {
  UseMethod("H2_Piepho")
}

#' @noRd
#' @export
H2_Delta <- function(
    model,
    target,
    type = c("BLUP", "BLUE"),
    aggregate = c("arithmetic", "harmonic"),
    options) {
  UseMethod("H2_Delta")
}

#' @noRd
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

#' @noRd
#' @export
H2_Delta_by_genotype <- function(model, target, type = c("BLUE", "BLUP"), options) {
  UseMethod("H2_Delta_by_genotype")
}

#' @noRd
#' @export
H2_Delta_by_genotype.default <- function(model,
                                         target = NULL,
                                         type = NULL,
                                         options = NULL) {
  H2D_ij <- H2_Delta_pairwise(model, target, type, options)

  H2D_i <- as.matrix(H2D_ij) |>
    rowMeans(na.rm = TRUE) |>
    data.frame()

  H2D_i <- setNames(H2D_i, "H2D_i")

  H2D_i_list <- split(H2D_i, rownames(H2D_i))

  return(H2D_i_list)
}

#' @noRd
#' @export
H2_Delta_pairwise <- function(model, target, type = c("BLUE", "BLUP"), options) {
  UseMethod("H2_Delta_pairwise")
}


