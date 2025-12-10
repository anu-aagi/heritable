#' Calculate broad-sense or narrow sense heritability
#' @description
#' A case-specific wrapper for calculating broad / narrow sense heritability.
#'
#' - The lowercase prefix `h2_` refers to the wrapper or subfunctions e.g. [`h2_Oakey()`] for calculating narrow sense heritability
#' - The upper case prefix `H2_` refers to the wrapper or subfunctions e.g. [`H2_Delta()`] for calculating broad sense heritability
#' @param model Model object of class `lmerMod/merMod` or `asreml`
#' @param method Character vector of name of method to calculate heritability. See details.
#' @param target The name of the random effect for which heritability is to be calculated.
#' @param options NULL by default, for internal checking of model object before calculations
#' @aliases H2
#' @usage
#' h2(model, target, method = c("Oakey", "Delta"), options)
#' H2(model, target, method = c("Cullis", "Oakey", "Delta", "Piepho", "Standard"), options)
#' @return A named numeric vector, length matching number of methods supplied
#' @details
#'
#' The following methods are currently implemented for narrow-sense heritability `h2(method = "XX")`:
#'
#' - `"Oakey"`: \deqn{H^2_{Oakey} = \frac{\sum_{i = n_z+1}^{n_g} \lambda_i}{\sum_{n_g}^{\lambda_i\neq 0}}}
#' - `"Delta"`: \deqn{H^2_{\Delta ..} = 1 - \frac{PEV^{BLUP}_{\overline\Delta ..}}{2\sigma^2_g}}
#'
#' The following methods are currently implemented for broad-sense heritability `H2(method = "XX")`:
#'
#' - `"Cullis"`: \deqn{H^2_{Cullis} = 1 - \frac{PEV^{BLUP}_{\overline\Delta ij}}{2\sigma^2_g}}
#' - `"Oakey"`: \deqn{H^2_{Oakey} = \frac{\sum_{i = n_z+1}^{n_g} \lambda_i}{\sum_{n_g}^{\lambda_i\neq 0}}}
#' - `"Delta"`: \deqn{H^2_{\Delta ..} = 1 - \frac{PEV^{BLUP}_{\overline\Delta ..}}{2\sigma^2_g}}
#' - `"Piepho"`: \deqn{H^2_{Piepho} = \frac{\sigma^2_g}{\sigma^2_g + \overline{PEV_{BLUE_g}} / 2}}
#' - `"Standard"`: \deqn{H^2_{Standard} = \frac{\sigma^2_g}{\sigma^2_g + \frac{1}{n_g}\sum_{n_g}^{i=1} \sigma^2_p / n_{gi}}}
#'
#' For further details of a specific method - take a look at helpfile for each subfunctions `?H2_Cullis`
#'
#' @references
#' - Cullis, B. R., Smith, A. B., & Coombes, N. E. (2006). On the design of early generation variety trials with correlated data. Journal of Agricultural, Biological, and Environmental Statistics, 11(4), 381–393. https://doi.org/10.1198/108571106X154443
#' - Oakey, H., Verbyla, A., Pitchford, W., Cullis, B., & Kuchel, H. (2006). Joint modeling of additive and non-additive genetic line effects in single field trials. Theoretical and Applied Genetics, 113(5), 809–819. https://doi.org/10.1007/s00122-006-0333-z
#' - Schmidt, P., Hartung, J., Rath, J., & Piepho, H.-P. (2019). Estimating Broad-Sense Heritability with Unbalanced Data from Agricultural Cultivar Trials. Crop Science, 59(2), 525–536. https://doi.org/10.2135/cropsci2018.06.0376
#' - Piepho, H.-P., & Möhring, J. (2007). Computing Heritability and Selection Response From Unbalanced Plant Breeding Trials. Genetics, 177(3), 1881–1888. https://doi.org/10.1534/genetics.107.074229
#' - Falconer, D. S., & Mackay, T. F. C. (1996). Introduction to quantitative genetics (4th ed.). Longman.
#' @seealso [H2_Cullis()], [H2_Oakey()], [H2_Delta()], [H2_Piepho()], [H2_Standard()], [`h2_Oakey()`], [`h2_Delta()`]
#' @export

h2 <- function(model, target, method, options) {
  UseMethod("h2")
}

#' @noRd
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
           Delta = h2_Delta(model, target, options = list(check = FALSE)),
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

#' @title Calculate Oakey's heritability
#' @description
#' Compute heritability for genotype means using the variance–covariance matrix of the genotype BLUPs
#' as described by Oakey et al. (2006).
#' @inheritParams h2
#' @aliases H2_Oakey
#' @details
#' \deqn{H^2_{Oakey} = \frac{\sum_{i = n_z+1}^{n_g} \lambda_i}{\sum_{n_g}^{\lambda_i\neq 0}}}
#' where:
#' - \eqn{n_g} is the number of genotypes
#' - \eqn{n_z} is the number of zero eigenvalues
#' - \eqn{\lambda_i} is the ith eigenvalue of the matrix \eqn{I_{m} - G^{-1}C^{gg}}
#' - \eqn{\sigma^2} is the variance attributed to differences between genotype
#'
#' See pages 813 and 818 of the reference for full derivation and explanation for Oakey's heritability
#' @usage
#' h2_Oakey(model, target, method, options)
#' H2_Oakey(model, target, method, options)
#' @references
#' Oakey, H., Verbyla, A., Pitchford, W., Cullis, B., & Kuchel, H. (2006). Joint modeling of additive and non-additive genetic line effects in single field trials. Theoretical and Applied Genetics, 113(5), 809–819. https://doi.org/10.1007/s00122-006-0333-z
#'
#' @export
h2_Oakey <- function(model, target, options) {
  UseMethod("h2_Oakey")
}

#' Calculate average heritability of differences between genotypes
#' @inheritParams h2
#' @aliases H2_Delta
#' @param type character, whether heritability is calculated using BLUEs or BLUPs
#' @param aggregate character, when taking means in the calculation, should harmonic or arithmetic mean be used?
#' @param options NULL by default, for internal checking of model object before calculations
#' @usage
#' h2_Delta(model, target, type = c("BLUE", "BLUP"), aggregate = c("arithmetic", "harmonic"))
#' H2_Delta(model, target, type = c("BLUE", "BLUP"), aggregate = c("arithmetic", "harmonic"))
#' @details
#' The heritability of differences between genotypes is given by:
#'
#' \deqn{H^2_{\Delta ..} = 1 - \frac{PEV^{BLUP}_{\overline\Delta ..}}{2\sigma^2_g}}
#'
#' where:
#'
#' - \eqn{PEV} is the prediction error variance matrix of the pairwise differences among BLUPs (BLUEs if `method = "BLUE"`)
#' - \eqn{\sigma^2} is the variance attributed to differences between genotype
#'
#' See reference page 995 - 997 for full derivation of this heritability measure and related variants
#' @references
#' Schmidt, P., Hartung, J., Rath, J., & Piepho, H.-P. (2019). Estimating Broad-Sense Heritability with Unbalanced Data from Agricultural Cultivar Trials. Crop Science, 59(2), 525–536. https://doi.org/10.2135/cropsci2018.06.0376
#' @seealso [`h2_Delta_by_genotype()`], [`H2_Delta_by_genotype()`], [`h2_Delta_pairwise()`], [`H2_Delta_pairwise()`]
#' @export
h2_Delta <- function(model,
                     target = NULL,
                     type = c("BLUE", "BLUP"),
                     aggregate = c("arithmetic", "harmonic"), options) {
  UseMethod("h2_Delta")
}

#' @noRd
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

#' Calculate heritability of differences for a given genotype
#' @inheritParams h2_Delta
#' @aliases H2_Delta_by_genotype
#' @export
h2_Delta_by_genotype <- function(model, target, type = c("BLUE", "BLUP"), options) {
  UseMethod("h2_Delta_by_genotype")
}

#' @noRd
#' @export
h2_Delta_by_genotype.default <- function(model,
                                         target = NULL,
                                         type = c("BLUP", "BLUE"),
                                         options = NULL) {
  type <- match.arg(type)

  H2D_ij <- h2_Delta_pairwise(model, target, type = type, options)

  H2D_i <- as.matrix(H2D_ij) |>
    rowMeans(na.rm = TRUE) |>
    data.frame()

  H2D_i <- setNames(H2D_i, "H2D_i")

  H2D_i_list <- split(H2D_i, rownames(H2D_i))

  return(H2D_i_list)
}

#' Calculate pairwise heritability of differences between genotypes
#' @inheritParams h2_Delta
#' @aliases H2_Delta_pairwise
#' @export
h2_Delta_pairwise <- function(model, target, type = c("BLUE", "BLUP"), options) {
  UseMethod("h2_Delta_pairwise")
}
