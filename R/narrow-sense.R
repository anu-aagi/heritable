#' Calculate broad-sense or narrow sense heritability from model object
#' @description
#' A case-specific wrapper for calculating broad / narrow sense heritability.
#'
#' - The lowercase prefix `h2_` refers to the wrapper or subfunctions e.g. [`h2_Oakey()`] for calculating narrow sense heritability
#' - The upper case prefix `H2_` refers to the wrapper or subfunctions e.g. [`H2_Delta()`] for calculating broad sense heritability
#' @param model Model object of class `lmerMod/merMod` or `asreml`
#' @param method Character vector of name of method to calculate heritability. See details.
#' @param target The name of the random effect for which heritability is to be calculated.
#' @param options NULL by default, for internal checking of model object before calculations
#' @param marginal Logical; if `TRUE`, construct marginal (strata-averaged)
#'   mappings so that each genotype receives a single averaged effect per term.
#'   If `FALSE`, mappings will only consider the main genotype effect and ignore the
#'   iteracting terms.
#' @param stratification A one-row data frame defining the stratum in which
#'   genotype effects should be evaluated. The columns must correspond
#'   to model terms that interact with `target`.
#' @param source The known genomic relationship matrix (GRM) used in `model` fitted using `asreml::vm()`, provided as a named list.
#' When not provided (an empty list by default), the GRM variable used for `vm` calling will be searched in the global environment.
#' Ignored for broad-sense and `lmerMod` methods
#' @param vc A list of precomputed variance components. Should be in the same structure as the output of [`var_comp()`]
#' @param ... Additional arguments that specify heritability calculation when interactions with genotype effects are modelled
#' @usage
#' h2(model,
#'    target,
#'    method = c("Cullis", "Oakey", "Piepho", "Delta", "Standard"),
#'    options = NULL,
#'    marginal = TRUE,
#'    stratification = NULL,
#'    source = list(),
#'    vc = NULL,
#'    ...)
#'
#' H2(model,
#'    target,
#'    method = c("Cullis", "Oakey", "Piepho", "Delta", "Standard"),
#'    options = NULL,
#'    marginal = TRUE,
#'    stratification = NULL,
#'    source = list(),
#'    vc = NULL,
#'    ...
#'    )
#' @name H2
#' @aliases H2, h2
#' @returns A named numeric vector, length matching number of methods supplied
#' @details
#'
#' The following methods are currently implemented for narrow-sense heritability `h2(method = "XX")`:
#' - `"Cullis"`: \deqn{H^2_{Cullis} = 1 - \frac{PEV^{BLUP}_{\overline\Delta ..}}{2\sigma^2_g}}
#' - `"Oakey"`: \deqn{H^2_{Oakey} = \frac{\sum_{i = n_z+1}^{n_g} \lambda_i}{\sum_{n_g}^{\lambda_i\neq 0}}}
#' - `"Piepho"`: \deqn{H^2_{Piepho} = \frac{\sigma^2_g}{\sigma^2_g + \overline{PEV_{BLUE_g}} / 2}}
#' - `"Delta"`: \deqn{H^2_{\Delta ij} = 1 - \frac{PEV^{BLUP}_{\overline\Delta ij}}{\operatorname{Var}(g_i - g_j)}}
#' - `"Standard"`: \deqn{H^2_{Standard} = \frac{\operatorname{Var}(g_i - g_j)}{\operatorname{Var}(y_i.. - y_j..)}}
#' - `"Reliability"`: \deqn{\bar{r}^2 = \frac{1}{n_g}\sum_{i=1}^{n_g} \left(1 - \frac{var(\hat{g}^{BLUP}_i)}{var(g_i)}\right)}
#'
#' `"Reliability"` is not part of the default `method` set but can be requested
#' explicitly. Use [`h2_Reliability_by_genotype()`] / [`H2_Reliability_by_genotype()`]
#' to obtain the per-genotype values \eqn{r^2_i} instead of their mean.
#'
#' The following methods are currently implemented for broad-sense heritability `H2(method = "XX")`:
#' - `"Cullis"`: \deqn{H^2_{Cullis} = 1 - \frac{PEV^{BLUP}_{\overline\Delta ..}}{2\sigma^2_g}}
#' - `"Oakey"`: \deqn{H^2_{Oakey} = \frac{\sum_{i = n_z+1}^{n_g} \lambda_i}{\sum_{n_g}^{\lambda_i\neq 0}}}
#' - `"Piepho"`: \deqn{H^2_{Piepho} = \frac{\sigma^2_g}{\sigma^2_g + \overline{PEV_{BLUE_g}} / 2}}
#' - `"Delta"`: \deqn{H^2_{\Delta ij} = 1 - \frac{PEV^{BLUP}_{\overline\Delta ij}}{2\sigma^2_g}}
#' - `"Standard"`: \deqn{H^2_{Standard} = \frac{\sigma^2_g}{\sigma^2_g + \frac{1}{n_g}\sum_{n_g}^{i=1} \sigma^2_p / n_{gi}}}
#' - `"Reliability"`: \deqn{\bar{r}^2 = \frac{1}{n_g}\sum_{i=1}^{n_g} \left(1 - \frac{var(\hat{g}^{BLUP}_i)}{var(g_i)}\right)}
#'
#' For further details of a specific method - take a look at helpfile for each subfunctions `?H2_Cullis`
#'
#' @references
#' - Cullis, B. R., Smith, A. B., & Coombes, N. E. (2006). On the design of early generation variety trials with correlated data. Journal of Agricultural, Biological, and Environmental Statistics, 11(4), 381–393. https://doi.org/10.1198/108571106X154443
#' - Oakey, H., Verbyla, A., Pitchford, W., Cullis, B., & Kuchel, H. (2006). Joint modeling of additive and non-additive genetic line effects in single field trials. Theoretical and Applied Genetics, 113(5), 809–819. https://doi.org/10.1007/s00122-006-0333-z
#' - Schmidt, P., Hartung, J., Rath, J., & Piepho, H.-P. (2019). Estimating Broad-Sense Heritability with Unbalanced Data from Agricultural Cultivar Trials. Crop Science, 59(2), 525–536. https://doi.org/10.2135/cropsci2018.06.0376
#' - Piepho, H.-P., & Möhring, J. (2007). Computing Heritability and Selection Response From Unbalanced Plant Breeding Trials. Genetics, 177(3), 1881–1888. https://doi.org/10.1534/genetics.107.074229
#' - Falconer, D. S., & Mackay, T. F. C. (1996). Introduction to quantitative genetics (4th ed.). Longman.
#' @seealso [H2_Cullis()], [H2_Oakey()], [H2_Delta()], [H2_Piepho()], [H2_Standard()], [h2_Oakey()], [h2_Delta()], [h2_Standard()]
#' @examples
#' # lme4 model
#' lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
#' lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
#' H2(lettuce_lme4, target = "gen", method = c("Standard", "Delta"))
#'
#' # asreml model (Requires license)
#' \dontrun{
#' lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
#'                                  random = ~ gen,
#'                                  data = lettuce_subset,
#'                                  trace = FALSE
#'                                  )
#'
#' H2(lettuce_asreml, target = "gen", method = c("Standard", "Delta"))
#' }
#' @export
h2 <- function(model,
               target,
               method = c("Cullis", "Oakey", "Piepho", "Delta", "Standard"),
               options = NULL,
               marginal = TRUE,
               stratification = NULL,
               source = list(),
               vc = NULL,
               ...) {
  UseMethod("h2")
}

#' @noRd
#' @export
h2.default <- function(model,
                       target,
                       method = c("Cullis", "Oakey", "Piepho", "Delta", "Standard"),
                       options = NULL,
                       marginal = TRUE,
                       stratification = NULL,
                       source = list(),
                       vc = NULL,
                       ...) {
  method <- match.arg(
    method,
    c("Cullis", "Oakey", "Piepho", "Delta", "Standard", "Reliability"),
    several.ok = TRUE
  )

  initial_checks(model, target, options = options)

  # Check correct model specification.
  if(options$check %||% TRUE){
    check_model_specification(model, target, "narrow_sense", source)
  }

  # Check design exists
  model <- check_design_exists(model, source = source)

  if(is.null(vc)){
    # Build variance component
    vc <- var_comp(model,
                   target = target,
                   source = source,
                   marginal = marginal,
                   stratification = stratification,
                   ...)
  }

  h2_values <- sapply(method, function(m) {
    switch(m,
           Cullis = h2_Cullis(model, target, options = list(check = FALSE), vc = vc, ...),
           Oakey = h2_Oakey(model, target, options = list(check = FALSE), vc = vc, ...),
           Piepho = h2_Piepho(model, target, options = list(check = FALSE), vc = vc, ...),
           Delta = h2_Delta(model, target, options = list(check = FALSE), vc = vc, ...),
           Standard = h2_Standard(model, target, options = list(check = FALSE), vc = vc, ...),
           Reliability = h2_Reliability(model, target, options = list(check = FALSE), vc = vc, ...),
           cli::cli_abort(
             "{.fn h2} is not implemented for method {.value m} of class{?es} {.code {class(model)}}"
           )
    )
  })

  # Set names and class
  h2_values <- stats::setNames(h2_values, method)
  args <- list(
    model = model,
    target = target,
    method = method,
    source = source,
    options = options,
    marginal = marginal,
    stratification = stratification
  )
  dots <- list(...)
  if(length(dots) > 0){
    args <- c(args, dots)
  }

  structure(h2_values,
            class = c("heritable", class(h2_values)),
            type = "narrow_sense",
            args = args
  )
}

#' Calculate standard heritability from model object
#' @description Compute standard heritability using the classic ratio method of
#' genotypic and phenotypic variance. See Falconer & Mackay (1996)
#' @usage
#' h2_Standard(model,
#'             target,
#'             options = NULL,
#'             marginal = TRUE,
#'             stratification = NULL,
#'             vc = NULL,
#'             ...)
#' H2_Standard(model,
#'             target,
#'             options = NULL,
#'             marginal = TRUE,
#'             stratification = NULL,
#'             vc = NULL,
#'             ...)
#' @inheritParams H2
#' @name H2_Standard
#' @aliases H2_Standard, h2_Standard
#' @return Numeric value
#' @details
#' The equation used to calculate standard heritability (broad-sense) is:
#' \deqn{H^2_{Standard} = \frac{\sigma^2_g}{\sigma^2_g + \frac{1}{n_g}\sum_{n_g}^{i=1} \sigma^2_p / n_{gi}}}
#' where:
#' - \eqn{n_g} is the number of genotypes
#' - \eqn{n_{gi}} is the number of replicate for a given genotype i
#' - \eqn{\sigma_g} is the variance attributed to genotype differences
#' - \eqn{\sigma_p} is the variance attributed to phenotypic differences
#'
#' The equation used to calculate standard heritability (narrow-sense) is:
#' \deqn{h^2_{Standard} = \frac{\operatorname{Var}(g_i - g_j)}{\operatorname{Var}(y_i.. - y_j..)}}
#' where:
#' - \eqn{g_i} is the random effect of the \eqn{i^{th}} genotype
#' - \eqn{y_i..} is the sample average of the  \eqn{i^{th}} genotype
#'
#' @export
#' @references Falconer, D. S., & Mackay, T. F. C. (1996). Introduction to quantitative genetics (4th ed.). Longman.
#' @seealso [H2_Standard()], [h2_Standard()]
#' @examples
#' # lme4 model
#' lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
#' lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
#' H2_Standard(lettuce_lme4, target = "gen")
#'
#' # asreml model (Requires license)
#' \dontrun{
#' lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
#'                                  random = ~ gen,
#'                                  data = lettuce_subset,
#'                                  trace = FALSE
#'                                  )
#'
#' H2_Standard(lettuce_asreml, target = "gen")
#' }
#' @export
h2_Standard <- function(model,
                        target,
                        options = NULL,
                        marginal = TRUE,
                        stratification = NULL,
                        vc = NULL,
                        ...) {
  UseMethod("h2_Standard")
}

#' @title Calculate Oakey's heritability from model object
#' @description
#' Compute heritability for genotype means using the variance–covariance matrix of the genotype BLUPs
#' as described by Oakey et al. (2006).
#' @inheritParams H2
#' @name H2_Oakey
#' @aliases H2_Oakey, h2_Oakey
#' @details
#' \deqn{h^2_{Oakey} = \frac{\sum_{i = n_z+1}^{n_g} \lambda_i}{\sum_{n_g}^{\lambda_i\neq 0}}}
#' where:
#' - \eqn{n_g} is the number of genotypes
#' - \eqn{n_z} is the number of zero eigenvalues
#' - \eqn{\lambda_i} is the ith eigenvalue of the matrix \eqn{I_{m} - G^{-1}C^{gg}}
#' - \eqn{\sigma^2} is the variance attributed to differences between genotype
#'
#' See pages 813 and 818 of the reference for full derivation and explanation for Oakey's heritability
#' @usage
#' h2_Oakey(model,
#'             target,
#'             options = NULL,
#'             marginal = TRUE,
#'             stratification = NULL,
#'             vc = NULL,
#'             ...)
#' H2_Oakey(model,
#'             target,
#'             options = NULL,
#'             marginal = TRUE,
#'             stratification = NULL,
#'             vc = NULL,
#'             ...)
#' @returns Numeric
#' @references
#' Oakey, H., Verbyla, A., Pitchford, W., Cullis, B., & Kuchel, H. (2006). Joint modeling of additive and non-additive genetic line effects in single field trials. Theoretical and Applied Genetics, 113(5), 809–819. https://doi.org/10.1007/s00122-006-0333-z
#' @examples
#' # lme4 model
#' lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
#' lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
#' H2_Oakey(lettuce_lme4, target = "gen")
#' @seealso [H2_Oakey()], [h2_Oakey()]
#' @export
h2_Oakey <- function(model,
                     target,
                     options  = NULL,
                     marginal = TRUE,
                     stratification = NULL,
                     vc = NULL,
                     ...) {
  UseMethod("h2_Oakey")
}

#' @title Calculate reliability from model object
#' @description
#' Compute the reliability (\eqn{r^2}, also known as the coefficient of
#' determination) of the genotype BLUPs, as described by Schmidt et al. (2019).
#' `h2_Reliability()` / `H2_Reliability()` return the overall (mean) reliability
#' \eqn{\bar{r}^2}, while `h2_Reliability_by_genotype()` /
#' `H2_Reliability_by_genotype()` return the per-genotype values \eqn{r^2_i}.
#' @inheritParams H2
#' @name H2_Reliability
#' @aliases H2_Reliability h2_Reliability H2_Reliability_by_genotype h2_Reliability_by_genotype
#' @details
#' The reliability of the \eqn{i}th genotype is
#' \deqn{r^2_i = 1 - \frac{var(\hat{g}^{BLUP}_i)}{var(g_i)}}
#' where \eqn{var(\hat{g}^{BLUP}_i)} is the \eqn{i}th diagonal element of the
#' prediction error variance matrix \eqn{C_{22(g)}} and \eqn{var(g_i)} is the
#' \eqn{i}th diagonal element of the genotypic variance-covariance matrix
#' \eqn{G_{(g)}}. As only the diagonal elements are used, \eqn{r^2_i} ignores the
#' off-diagonal elements (covariances) of both matrices. The overall reliability
#' is the mean across genotypes,
#' \deqn{\bar{r}^2 = \frac{1}{n_g}\sum_{i=1}^{n_g} r^2_i.}
#' @returns `h2_Reliability()` / `H2_Reliability()` return a numeric value;
#'   `h2_Reliability_by_genotype()` / `H2_Reliability_by_genotype()` return a
#'   named numeric vector.
#' @references
#' - Schmidt, P., Hartung, J., Bennewitz, J., & Piepho, H.-P. (2019). Heritability in Plant Breeding on a Genotype-Difference Basis. Genetics, 212(4), 991–1008. https://doi.org/10.1534/genetics.119.302134
#' - Mrode, R. A. (2014). Linear Models for the Prediction of Animal Breeding Values (3rd ed.). CABI.
#' @seealso [H2_Reliability_parameters()]
#' @examples
#' # lme4 model
#' lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
#' lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
#' H2_Reliability(lettuce_lme4, target = "gen")
#' H2_Reliability_by_genotype(lettuce_lme4, target = "gen")
#' @export
h2_Reliability <- function(model,
                           target,
                           options = NULL,
                           marginal = TRUE,
                           stratification = NULL,
                           vc = NULL,
                           ...) {
  UseMethod("h2_Reliability")
}

#' @noRd
#' @export
h2_Reliability.default <- function(model,
                                   target,
                                   options = NULL,
                                   marginal = TRUE,
                                   stratification = NULL,
                                   vc = NULL,
                                   ...) {
  r2 <- h2_Reliability_by_genotype(model,
                                   target,
                                   options = options,
                                   marginal = marginal,
                                   stratification = stratification,
                                   vc = vc,
                                   ...)

  if (is.atomic(r2) && length(r2) == 1 && is.na(r2)) {
    return(NA)
  }

  mean(r2)
}

#' @rdname H2_Reliability
#' @export
h2_Reliability_by_genotype <- function(model,
                                       target,
                                       options = NULL,
                                       marginal = TRUE,
                                       stratification = NULL,
                                       vc = NULL,
                                       ...) {
  UseMethod("h2_Reliability_by_genotype")
}

#' Calculate average heritability of differences between genotypes from model object
#' @description
#' Instead of computing heritability on a "entry-mean" basis, this method
#' calculates heritability using "entry-differences". Entry here is
#' referring to the genotype, line or variety of interest. See
#' reference for origin and interpretation of `h2/H2_Delta` and it's variants
#' @inheritParams H2
#' @name H2_Delta
#' @aliases H2_Delta, h2_Delta
#' @param type character, whether heritability is calculated using BLUEs or BLUPs
#' @param options NULL by default, for internal checking of model object before calculations
#' @usage
#' h2_Delta(model,
#'          target,
#'          type = c("BLUP", "BLUE"),
#'          options = NULL,
#'          marginal = TRUE,
#'          stratification = NULL,
#'          vc = NULL,
#'          ...)
#'
#' H2_Delta(model,
#'          target,
#'          type = c("BLUP", "BLUE"),
#'          options = NULL,
#'          marginal = TRUE,
#'          stratification = NULL,
#'          vc = NULL,
#'          ...)
#' @returns Numeric
#' @details
#' The broad-sense heritability of differences between genotypes is given by:
#'
#' \deqn{H^2_{\Delta ..} = 1 - \frac{PEV^{BLUP}_{\overline\Delta ..}}{2\sigma^2_g}}
#'
#' where:
#' - \eqn{PEV^{BLUP}_{\overline\Delta ..}} is the mean of the prediction error variance matrix for the pairwise differences among BLUPs (BLUEs if `method = "BLUE"`) across all genotypes
#' - \eqn{\sigma^2} is the variance attributed to differences between genotype
#'
#' The narrow-sense heritability of differences between genotypes is given by:
#'
#' \deqn{h^2_{\Delta ij} = 1 - \frac{PEV^{BLUP}_{\overline\Delta ij}}{\operatorname{Var}(g_i - g_j)}}
#'
#' where:
#' - \eqn{g_i} is the random effect of the \eqn{i^{th}} genotype
#'
#' See reference page 995 - 997 for full derivation of this heritability measure and related variants
#' @references
#' Schmidt, P., Hartung, J., Rath, J., & Piepho, H.-P. (2019). Estimating Broad-Sense Heritability with Unbalanced Data from Agricultural Cultivar Trials. Crop Science, 59(2), 525–536. https://doi.org/10.2135/cropsci2018.06.0376
#' @seealso [`h2_Delta_by_genotype()`], [`H2_Delta_by_genotype()`], [`h2_Delta_pairwise()`], [`H2_Delta_pairwise()`]
#' @export
#' @examples
#' # lme4 model
#' lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
#' lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
#' H2_Delta(lettuce_lme4, target = "gen", type = "BLUP")
#'
#' # asreml model (Requires license)
#' \dontrun{
#' lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
#'                                  random = ~ gen,
#'                                  data = lettuce_subset,
#'                                  trace = FALSE
#'                                  )
#'
#' H2_Delta(lettuce_asreml, target = "gen", type = "BLUP")
#' }
h2_Delta <- function(model,
                     target,
                     type = c("BLUP", "BLUE"),
                     options = NULL,
                     marginal = TRUE,
                     stratification = NULL,
                     vc = NULL,
                     ...) {
  UseMethod("h2_Delta")
}

#' @noRd
#' @export
h2_Delta.default <- function(model,
                             target,
                             type = c("BLUP", "BLUE"),
                             options = NULL,
                             marginal = TRUE,
                             stratification = NULL,
                             vc = NULL,
                             ...) {
  type <- match.arg(type)

  h2D_ij <- h2_Delta_pairwise(model,
                              target,
                              type = type,
                              options = options,
                              marginal =  marginal,
                              stratification = stratification,
                              vc = vc,
                              ...)

  if (is.atomic(h2D_ij) && length(h2D_ij) == 1 && is.na(h2D_ij)) {
    return(NA)
  }

  delta_g <- mean(attr(h2D_ij, "delta_g")[upper.tri(h2D_ij)])
  delta_pev <- mean(attr(h2D_ij, "delta_pev")[upper.tri(h2D_ij)])

  switch(type,
         "BLUP" = 1 - delta_pev / delta_g,
         "BLUE" = delta_g / (delta_g + delta_pev)
  )
}

#' Calculate heritability of differences for a given genotype from model object
#' @description
#' Instead of computing heritability on a "entry-mean" basis, this method
#' calculates heritability using "entry-differences". Entry here is
#' referring to the genotype, line or variety of interest. See
#' reference for origin and interpretation of `h2/H2_Delta_by_genotype` and it's variants
#' @usage
#' h2_Delta_by_genotype(model,
#'                      target,
#'                      type = c("BLUP", "BLUE"),
#'                      options = NULL,
#'                      marginal = TRUE,
#'                      stratification = NULL,
#'                      vc = NULL,
#'                      ...)
#' H2_Delta_by_genotype(model,
#'                      target,
#'                      type = c("BLUP", "BLUE"),
#'                      options = NULL,
#'                      marginal = TRUE,
#'                      stratification = NULL,
#'                      vc = NULL,
#'                      ...)
#' @inheritParams H2_Delta
#' @name H2_Delta_by_genotype
#' @aliases H2_Delta_by_genotype, h2_Delta_by_genotype
#' @returns Numeric
#' @details
#' The broad-sense heritability of differences between genotypes is given by:
#'
#' \deqn{H^2_{\Delta ..} = 1 - \frac{PEV^{BLUP}_{\overline\Delta ..}}{2\sigma^2_g}}
#'
#' where:
#' - \eqn{PEV^{BLUP}_{\overline\Delta ..}} is the mean of the prediction error variance matrix for the pairwise differences among BLUPs (BLUEs if `method = "BLUE"`) across all genotypes
#' - \eqn{\sigma^2} is the variance attributed to differences between genotype
#'
#' The narrow-sense heritability of differences between genotypes is given by:
#'
#' \deqn{h^2_{\Delta ij} = 1 - \frac{PEV^{BLUP}_{\overline\Delta ij}}{\operatorname{Var}(g_i - g_j)}}
#'
#' where:
#' - \eqn{g_i} is the random effect of the \eqn{i^{th}} genotype
#'
#' See reference page 995 - 997 for full derivation of this heritability measure and related variants
#' @references
#' Schmidt, P., Hartung, J., Rath, J., & Piepho, H.-P. (2019). Estimating Broad-Sense Heritability with Unbalanced Data from Agricultural Cultivar Trials. Crop Science, 59(2), 525–536. https://doi.org/10.2135/cropsci2018.06.0376
#' @seealso [`h2_Delta()`], [`H2_Delta()`], [`h2_Delta_pairwise()`], [`H2_Delta_pairwise()`]
#' @returns Named list, with each element containing a named numeric vector
#' @export
#' @examples
#' # lme4 model
#' lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
#' lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
#' H2_Delta_by_genotype(lettuce_lme4, target = "gen", type = "BLUP")
#'
#' # asreml model (Requires license)
#' \dontrun{
#' lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
#'                                  random = ~ gen,
#'                                  data = lettuce_subset,
#'                                  trace = FALSE
#'                                  )
#'
#' H2_Delta_by_genotype(lettuce_asreml, target = "gen", type = "BLUP")
#' }
h2_Delta_by_genotype <- function(model,
                                 target,
                                 type = c("BLUP", "BLUE"),
                                 options = NULL,
                                 marginal = TRUE,
                                 stratification = NULL,
                                 vc = NULL,
                                 ...) {
  UseMethod("h2_Delta_by_genotype")
}

#' @noRd
#' @export
h2_Delta_by_genotype.default <- function(model,
                                         target,
                                         type = c("BLUP", "BLUE"),
                                         options = NULL,
                                         marginal = TRUE,
                                         stratification = NULL,
                                         vc = NULL,
                                         ...) {
  type <- match.arg(type)

  h2D_ij <- h2_Delta_pairwise(model,
                              target,
                              type = type,
                              options = options,
                              marginal =  marginal,
                              stratification = stratification,
                              vc = vc,
                              ...)

  delta_g <- attr(h2D_ij, "delta_g")
  delta_pev <- attr(h2D_ij, "delta_pev")
  delta_pev <- Matrix::rowSums(delta_pev)/(ncol(delta_pev)-1)

  if(type == "BLUP"){
    return(1 - delta_pev / delta_g)
  } else {
    return(delta_g / (delta_g + delta_pev))
  }
}

#' Calculate pairwise heritability of differences between genotypes from model object
#' @description
#' Instead of computing heritability on a "entry-mean" basis, this method
#' calculates heritability using "entry-differences". Entry here is
#' referring to the genotype, line or variety of interest. See
#' reference for origin and interpretation of `h2/H2_Delta_pairwise` and it's variants
#' @usage
#' h2_Delta_pairwise(model,
#'                   target,
#'                   type = c("BLUP", "BLUE"),
#'                   options = NULL,
#'                   marginal = TRUE,
#'                   stratification = NULL,
#'                   vc = NULL,
#'                   ...)
#' H2_Delta_pairwise(model,
#'                   target,
#'                   type = c("BLUP", "BLUE"),
#'                   options = NULL,
#'                   marginal = TRUE,
#'                   stratification = NULL,
#'                   vc = NULL,
#'                   ...)
#' @inheritParams H2_Delta
#' @name H2_Delta_pairwise
#' @aliases H2_Delta_pairwise, h2_Delta_pairwise
#' @returns A `dspMatrix`
#' @references
#' Schmidt, P., Hartung, J., Rath, J., & Piepho, H.-P. (2019). Estimating Broad-Sense Heritability with Unbalanced Data from Agricultural Cultivar Trials. Crop Science, 59(2), 525–536. https://doi.org/10.2135/cropsci2018.06.0376
#' @seealso [`h2_Delta_by_genotype()`], [`H2_Delta_by_genotype()`], [`h2_Delta()`], [`H2_Delta()`]
#' @export
#' @examples
#' # lme4 model
#' lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
#' lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
#' H2_Delta_pairwise(lettuce_lme4, target = "gen", type = "BLUP")
#'
#' # asreml model (Requires license)
#' \dontrun{
#' lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
#'                                  random = ~ gen,
#'                                  data = lettuce_subset,
#'                                  trace = FALSE
#'                                  )
#'
#' H2_Delta_pairwise(lettuce_asreml, target = "gen", type = "BLUP")
#' }
h2_Delta_pairwise <- function(model,
                              target,
                              type = c("BLUP", "BLUE"),
                              options = NULL,
                              marginal = TRUE,
                              stratification = NULL,
                              vc = NULL,
                              ...) {
  UseMethod("h2_Delta_pairwise")
}


#' Calculate Cullis' heritability from model object
#' @description Compute "generalised heritability" for unbalanced experimental designs.
#' See Cullis, Smith and Coombes (2006) for derivation.
#' @inheritParams H2
#' @name H2_Cullis
#' @aliases H2_Cullis, h2_Cullis
#' @usage
#' H2_Cullis(model,
#'           target,
#'           options = NULL,
#'           marginal = TRUE,
#'           stratification = NULL,
#'           vc = NULL,
#'           ...)
#' h2_Cullis(model,
#'           target,
#'           options = NULL,
#'           marginal = TRUE,
#'           stratification = NULL,
#'           vc = NULL,
#'           ...)
#' @return Numeric value
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
#' @examples
#' # lme4 model
#' lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
#' lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
#' H2_Cullis(lettuce_lme4, target = "gen")
#'
#' # asreml model (Requires license)
#' \dontrun{
#' lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
#'                                  random = ~ gen,
#'                                  data = lettuce_subset,
#'                                  trace = FALSE
#'                                  )
#'
#' H2_Cullis(lettuce_asreml, target = "gen")
#' }

h2_Cullis <- function(model,
                      target,
                      options = NULL,
                      marginal = TRUE,
                      stratification = NULL,
                      vc = NULL,
                      ...) {
  UseMethod("h2_Cullis")
}


#' Calculate Piepho's heritability from model object
#' Compute Piepho's heritability using variance differences between genotype BLUEs
#' @usage
#' H2_Piepho(model,
#'           target,
#'           options = NULL,
#'           marginal = TRUE,
#'           stratification = NULL,
#'           vc = NULL,
#'           ...)
#' h2_Piepho(model,
#'           target,
#'           options = NULL,
#'           marginal = TRUE,
#'           stratification = NULL,
#'           vc = NULL,
#'           ...)
#' @inheritParams H2
#' @name H2_Piepho
#' @aliases H2_Piepho, h2_Piepho
#' @details The equation for Piepho's heritability is as follows:
#'
#' \deqn{H^2_{Piepho} = \frac{\sigma^2_g}{\sigma^2_g + \overline{PEV_{BLUE_g}} / 2}}
#'
#' where:
#' - \eqn{\overline{PEV_{BLUE_g}}} is the prediction error variance matrix for genotype BLUEs
#' - \eqn{\sigma^2_g} is the variance attributed to differences between genotype
#'
#' See reference for full derivation and details.
#' @returns Numeric
#' @export
#' @references
#' Piepho, H.-P., & Möhring, J. (2007). Computing Heritability and Selection Response From Unbalanced Plant Breeding Trials. Genetics, 177(3), 1881–1888. https://doi.org/10.1534/genetics.107.074229
#' @examples
#' # lme4 model
#' lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
#' lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
#' H2_Piepho(lettuce_lme4, target = "gen")
#'
#' # asreml model (Requires license)
#' \dontrun{
#' lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
#'                                  random = ~ gen,
#'                                  data = lettuce_subset,
#'                                  trace = FALSE
#'                                  )
#'
#' H2_Piepho(lettuce_asreml, target = "gen")
#' }
h2_Piepho <- function(model,
                      target,
                      options = NULL,
                      marginal = TRUE,
                      stratification = NULL,
                      vc = NULL,
                      ...) {
  UseMethod("h2_Piepho")
}
