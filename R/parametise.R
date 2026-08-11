#' Calculate Standard heritability using variance parameters
#'
#' @description Compute Standard heritability for genotype means using the variance components of genotype and residuals.
#'
#' @details The equation for Standard heritability is as follows:
#'
#' \deqn{H^2_{Standard} = \frac{\sigma^2_g}{\sigma^2_g + \frac{1}{n_g}\sum_{n_g}^{i=1} \sigma^2_p / n_{gi}}}
#' where:
#' - \eqn{n_g} is the number of genotypes
#' - \eqn{n_{gi}} is the number of replicate for a given genotype i
#' - \eqn{\sigma_g} is the variance attributed to genotype differences
#' - \eqn{\sigma_p} is the variance attributed to phenotypic differences
#'
#' @param vc_g Numeric. Genotype variance component
#' @param vc_e Numeric. Residuals variance component
#' @param n_r A numeric vector of size n_g, the number of genotype replicates.
#' @return Numeric value
#' @examples
#' H2_Standard_parameters(vc_g = 0.25, vc_e = 0.8)
#'
#' @references
#' Falconer, D. S., & Mackay, T. F. C. (1996). Introduction to quantitative genetics (4th ed.). Longman.
#'
#' @export
H2_Standard_parameters <- function(vc_g, vc_e, n_r = 1) {
  H2_Standard <- vc_g / (vc_g +  mean(vc_e / n_r))
  return(H2_Standard)
}

#' Calculate Cullis heritability using variance parameters
#'
#' @description Compute the Cullis heritability for genotype means using the average
#' variance of pairwise differences of best linear unbiased predictors (BLUPs).
#'
#' @details The equation for Cullis heritability is as follow
#'
#' \deqn{H^2_{Cullis} = 1 - \frac{PEV^{BLUP}_{\overline\Delta ij}}{2\sigma^2_g}}
#'
#' where:
#' - \eqn{PEV} is the prediction error variance matrix of the pairwise differences among BLUPS
#' - \eqn{\sigma^2} is the variance attributed to differences between genotype
#'
#' @param vd_BLUP_avg Numeric. Average variance of pairwise differences among BLUPs
#' @param vc_g Numeric. Genotype variance component
#' @return Numeric value
#' @examples
#' H2_Cullis_parameters(vd_BLUP_avg = 0.25, vc_g = 0.8)
#'
#' @references
#' Cullis, B. R., Smith, A. B., & Coombes, N. E. (2006). On the design of early generation variety trials with correlated data. Journal of Agricultural, Biological, and Environmental Statistics, 11(4), 381–393. https://doi.org/10.1198/108571106X154443
#'
#' @export
H2_Cullis_parameters <- function(vd_BLUP_avg, vc_g){
  H2_Cullis <- 1 - (vd_BLUP_avg / 2 / vc_g)

  return(H2_Cullis)
}

#' Calculate Oakey's heritability using variance parameters
#'
#' @description Rather than providing a model object, supply the necessary components to compute
#' this heritability measure.
#'
#' @param Gg_inv The inverse of the genotypic variance-covariance matrix.
#' @param C22_g Prediction error variance matrix associated with the genotype effects.
#' @return Numeric value
#' @examples
#' Gg_inv = diag(1/0.15, 3, 3)
#' C22_g <- matrix(
#'   c(
#'     0.08, 0.01, 0.00,
#'     0.01, 0.07, 0.01,
#'     0.00, 0.01, 0.09
#'   ),
#'   nrow = 3, byrow = TRUE
#' )
#' H2_Oakey_parameters(Gg_inv, C22_g)
#' @export
H2_Oakey_parameters <- function(Gg_inv, C22_g) {
   n_g <- nrow(Gg_inv)
   svds <- svd(Gg_inv)
   Gg_inv_sqrt <- sweep(svds$u, 2, sqrt(svds$d), "*") %*% t(svds$v)
   M <- diag(n_g) - (Gg_inv_sqrt %*% C22_g %*% Gg_inv_sqrt)

   eM <- eigen(M, symmetric = TRUE)
   thres <- sqrt(.Machine$double.eps)
   H2_Oakey <- mean(eM$values[eM$values/max(eM$values) > thres])
   return(H2_Oakey)

}

#' Calculate reliability using variance parameters
#'
#' @description Compute the reliability (\eqn{r^2}, also known as the coefficient
#' of determination) for each genotype from the genotypic variance-covariance
#' matrix and the prediction error variance (PEV) matrix of the genotype BLUPs.
#'
#' @details The reliability of the \eqn{i}th genotype is
#'
#' \deqn{r^2_i = 1 - \frac{var(\hat{g}^{BLUP}_i)}{var(g_i)}}
#'
#' where:
#' - \eqn{var(\hat{g}^{BLUP}_i)} is the \eqn{i}th diagonal element of the PEV
#'   matrix \eqn{C_{22(g)}}
#' - \eqn{var(g_i)} is the \eqn{i}th diagonal element of the genotypic
#'   variance-covariance matrix \eqn{G_{(g)}}
#'
#' As only the diagonal elements are used, \eqn{r^2_i} does not account for the
#' off-diagonal elements (i.e. covariances) of either matrix. The overall
#' reliability is obtained as the mean across genotypes,
#' \eqn{\bar{r}^2 = \frac{1}{n_g}\sum_{i=1}^{n_g} r^2_i}.
#'
#' Genotypes whose genotypic variance (the corresponding diagonal element of
#' \eqn{G_{(g)}}) is not strictly positive yield an undefined reliability and are
#' returned as `NA` (with a warning) rather than `Inf`/`NaN`.
#'
#' @param G_g Genotypic variance-covariance matrix.
#' @param C22_g Prediction error variance matrix associated with the genotype effects.
#' @return A named numeric vector of per-genotype reliabilities. Entries for
#'   genotypes with a non-positive genotypic variance are `NA`.
#' @examples
#' G_g <- diag(0.15, 3, 3)
#' C22_g <- matrix(
#'   c(
#'     0.08, 0.01, 0.00,
#'     0.01, 0.07, 0.01,
#'     0.00, 0.01, 0.09
#'   ),
#'   nrow = 3, byrow = TRUE
#' )
#' H2_Reliability_parameters(G_g, C22_g)
#'
#' @references
#' - Schmidt, P., Hartung, J., Bennewitz, J., & Piepho, H.-P. (2019). Heritability in Plant Breeding on a Genotype-Difference Basis. Genetics, 212(4), 991–1008. https://doi.org/10.1534/genetics.119.302134
#' - Mrode, R. A. (2014). Linear Models for the Prediction of Animal Breeding Values (3rd ed.). CABI.
#'
#' @export
H2_Reliability_parameters <- function(G_g, C22_g) {
  var_g <- diag(as.matrix(G_g))
  pev <- diag(as.matrix(C22_g))

  r2 <- 1 - pev / var_g

  # Reliability r^2 = 1 - PEV / var(g) is undefined when the genotypic variance
  # is not strictly positive (e.g. a degenerate genotype in a relationship
  # matrix). Guard against silently propagating Inf/NaN into the overall (mean)
  # reliability by setting those genotypes to NA and warning.
  invalid <- !is.finite(var_g) | var_g <= 0
  if (any(invalid)) {
    r2[invalid] <- NA_real_
    n_invalid <- sum(invalid)
    cli::cli_warn(c(
      "Reliability is undefined for {n_invalid} genotype{?s} with a
       non-positive genotypic variance.",
      "i" = "Their reliability was set to {.val {NA}}."
    ))
  }

  r2
}

#' Calculate Piepho's heritability using variance parameters
#'
#' @description Compute Piepho's heritability using the variance of differences between two BLUES.
#'
#' @details The equation for Piepho's heritability is as follows:
#'
#' \deqn{H^2_{Piepho} = \frac{\sigma^2_g}{\sigma^2_g + \overline{PEV_{BLUE_g}} / 2}}
#'
#' where:
#' - \eqn{\overline{PEV_{BLUE_g}}} is the prediction error variance matrix for genotype BLUEs
#' - \eqn{\sigma^2_g} is the variance attributed to differences between genotype
#'
#' @param vc_g Numeric. Genotype variance component
#' @param vd_BLUE_avg Numeric. Mean variance of pairwise differences among BLUES
#' @return Numeric value
#' @examples
#' H2_Piepho_parameters(vc_g = 0.25, vd_BLUE_avg = 0.68)
#'
#' @references
#' Piepho, H.-P., & Möhring, J. (2007). Computing Heritability and Selection Response From Unbalanced Plant Breeding Trials. Genetics, 177(3), 1881–1888. https://doi.org/10.1534/genetics.107.074229
#'
#' @export
H2_Piepho_parameters <- function(vc_g, vd_BLUE_avg) {
  H2_Piepho <- vc_g / (vc_g + (vd_BLUE_avg / 2))
  return(H2_Piepho)
}

#' Calculate heritability of pairwise differences using variance parameters
#' @description Compute broad-sense heritability of differences
#' using the variance of differences between two BLUPs/BLUEs
#' @details See [H2_Delta()] and reference for full derivation
#'  and equation for heritability Delta
#' @param delta_g Numeric. Genotypic variance-covariance matrix.
#' @param delta_pev Matrix. Variance of pairwise differences among BLUES or BLUPs
#' @param type Character. Either BLUES or BLUPS used to compute the variance of pairwise differences.
#' @return Matrix of pairwise heritability of differences among BLUES or BLUPs
#' @examples
#'
#' H2_Delta_parameters(delta_g = diag(0.15, 2, 2),
#'                     delta_pev = matrix(c(NA,0.2,0.2,NA),2,2),
#'                     type = "BLUP"
#'                     )
#'
#' @references
#' Schmidt, P., Hartung, J., Rath, J., & Piepho, H.-P. (2019). Estimating
#' Broad-Sense Heritability with Unbalanced Data from Agricultural Cultivar
#' Trials. Crop Science, 59(2), 525–536. https://doi.org/10.2135/cropsci2018.06.0376
#' @export
H2_Delta_parameters <- function(delta_g, delta_pev, type = c("BLUP", "BLUE")) {
  type <- match.arg(type)

  if (type == "BLUP") {
    1 - delta_pev / delta_g
  } else {
    1 / (1 + delta_pev / delta_g)
  }
}

#' @noRd
#' @return Narrow sense heritability for the standard method.
#' @keywords internal
h2_Standard_parameters <- function(G_g, V, C, active_g) {

  # S = C' V C  (g x g)
  S <- crossprod(C, V %*% C)

  # Sample contrast var(y_i.. - y_j..) = den_ij
  # = (C_i - C_j)' V (C_i - C_j) = S_ii + S_jj - 2 S_ij
  delta_y <- var_diff(S)

  # genotype contrast var(g_i - g_j) = num_ij = G_ii + G_jj
  delta_g <- var_diff(G_g)
  # for narrow sense, use outer(diag(G_g), diag(G_g), "+") - 2 * G_g

  if(!is.null(active_g) && any(!active_g)){
    delta_y <- delta_y[active_g,active_g,drop=FALSE]
    delta_g <- delta_g[active_g,active_g,drop=FALSE]
  }

  mean(delta_g[upper.tri(delta_g)])/mean(delta_y[upper.tri(delta_y)])

}
