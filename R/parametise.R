#' Estimate Cullis heritability
#'
#' @description Compute the Cullis heritability for genotype means using the average
#' variance of pairwise differences of best linear unbiased predictors (BLUPs).
#' 
#' @details The equation for Cullis heritability is as follows:
#' 
#' \eqn{H^2 = 1 - (vd_BLUP_avg / (2 * vc_g))}
#'
#' @param vd_BLUP_avg Average variance of pairwise differences among BLUPs 
#' @param vc_g Genotype variance component
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
#' @param n_g Integer. Number of genotypes 
#' @param vc_g Numeric. Genotype variance component
#' @param vcov_g Numeric matrix. Variance–covariance matrix of the genotype BLUPs (prediction error variances and covariances). Expected to be an n_g by n_g symmetric matrix.
#'
#' #' @details The equation for Oakey heritability is as follows:
#' 
#' \eqn{H^2 = 1 - (vd_BLUP_avg / (2 * vc_g))}
#'
#' @return Single numeric value 
#'
#' @examples
#' # H2_Oakey_parameters(n_g = 50, vc_g = 1.2, vcov_g = diag(0.3, 50, 50))
#'
#' @references
#' Oakey, H., Verbyla, A., Pitchford, W., Cullis, B., & Kuchel, H. (2006). Joint modeling of additive and non-additive genetic line effects in single field trials. Theoretical and Applied Genetics, 113(5), 809–819. https://doi.org/10.1007/s00122-006-0333-z
#'
#' @export
H2_Oakey_parameters <- function(n_g, vc_g, vcov_g) {

  Gg_inv <- diag(1 / vc_g, nrow = n_g, ncol = n_g)
  M <- diag(n_g) - (Gg_inv %*% vcov_g)
  eM <- eigen(M)

  H2_Oakey <- sum(eM$values) / (n_g - 1)
  return(H2_Oakey)
}