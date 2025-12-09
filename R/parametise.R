#' Estimate Cullis heritability
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
#' @return Single numeric value
#'
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

#' Estimate Oakey's heritability
#'
#' @description Compute heritability for genotype means using the variance–covariance matrix of the genotype BLUPs
#' as described by Oakey et al. (2006).
#'
#' @param Gg_inv The inverse of the genotypic variance-covariance matrix.
#' @param C_gg Prediction error variance matrix associated with the genotype effects.
#'
#' @details See pages 813 and 818 of reference for full derivation and explanation for Oakey's heritability
#'
#' \deqn{H^2_{Oakey} = \frac{\sum_{i = n_z+1}^{n_g} \lambda_i}{\sum_{n_g}^{\lambda_i\neq 0}}}
#' where:
#' - \eqn{n_g} is the number of genotypes
#' - \eqn{n_z} is the number of zero eigenvalues
#' - \eqn{\lambda_i} is the ith eigenvalue of the matrix \eqn{I_{m} - G^{-1}C^{gg}}
#' - \eqn{\sigma^2} is the variance attributed to differences between genotype
#'
#' @return Single numeric value
#' @examples
#' Gg_inv = diag(1/0.15, 3, 3)
#' C_gg <- matrix(
#'   c(
#'     0.08, 0.01, 0.00,
#'     0.01, 0.07, 0.01,
#'     0.00, 0.01, 0.09
#'   ),
#'   nrow = 3, byrow = TRUE
#' )
#' H2_Oakey_parameters(Gg_inv, C_gg)
#'
#' @references
#' Oakey, H., Verbyla, A., Pitchford, W., Cullis, B., & Kuchel, H. (2006). Joint modeling of additive and non-additive genetic line effects in single field trials. Theoretical and Applied Genetics, 113(5), 809–819. https://doi.org/10.1007/s00122-006-0333-z
#'
#' @export
H2_Oakey_parameters <- function(Gg_inv, C_gg) {
   n_g <- nrow(Gg_inv)
   M <- diag(n_g) - (Gg_inv %*% C_gg)
   eM <- eigen(M)

   thres <- 1e-5
   H2_Oakey <- mean(eM$values[eM$values > thres])
   return(H2_Oakey)
}

#' Estimate Standard heritability
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
#' @return Single numeric value
#'
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

#' Estimate Piepho's heritability
#'
#' @description Compute Piepho's heritability using the variance of differences between two BLUES.
#'
#' @details The equation for Piepho's heritability is as follows:
#'
#' \eqn{H^2 = 1 - (vc_g / (vc_g + vd_BLUE_avg / 2))}
#'
#' @param vc_g Numeric. Genotype variance component
#' @param vd_BLUE_avg Numeric. Mean variance of pairwise differences among BLUES
#' @return Single numeric value
#'
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

#' Estimate heritability of differences (Delta) for BLUEs or BLUPs
#'
#' @description Compute heritability of differences using the variance of differences between two BLUES.
#'
#' @details See reference for full derivation and equation for heritability Delta BLUES
#'
#' @param var1 Numeric. The variance of genotype 1
#' @param var2 Numeric. The variance of genotype 2
#' @param cov Numeric. Covariance between genotypes 1 and 2.
#' @param vd_matrix Matrix. Variance of pairwise differences among BLUES or BLUPs
#' @return Matrix of pairwise heritability of differences among BLUES or BLUPs
#'
#' @examples
#' H2_Delta_BLUE_parameters(vc_g = 0.25, vd_matrix = matrix(c(NA,0.2,0.2,NA),2,2))
#'
#' @references
#' Schmidt, P., Hartung, J., Rath, J., & Piepho, H.-P. (2019). Estimating Broad-Sense Heritability with Unbalanced Data from Agricultural Cultivar Trials. Crop Science, 59(2), 525–536. https://doi.org/10.2135/cropsci2018.06.0376
#'
#' @export
H2_Delta_BLUE_parameters <- function(vc_g, vd_matrix) {
   denom <- 2 * vc_g
   1 / (1 + vd_matrix / denom)
}


#' @export
H2_Delta_BLUP_parameters <- function(vc_g, vd_matrix) {
  denom <- 2 * vc_g
  1 - vd_matrix / denom
}


#' @export
h2_Delta_BLUE_parameters <- function(G_g, vd_matrix) {
  vd <- diag(G_g)
  n_g <- nrow(G_g)
  denom <- matrix(vd, n_g, n_g) + matrix(vd, n_g, n_g, byrow = TRUE) - 2 * G_g
  1 / (1 + vd_matrix / denom)
}


#' @export
h2_Delta_BLUP_parameters <- function(G_g, vd_matrix) {
  vd <- diag(G_g)
  n_g <- nrow(G_g)
  denom <- matrix(vd, n_g, n_g) + matrix(vd, n_g, n_g, byrow = TRUE) - 2 * G_g
  1 - vd_matrix / denom
}



