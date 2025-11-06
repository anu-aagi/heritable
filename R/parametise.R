#' Estimate Cullis heritability
#'
#' @description Compute the Cullis heritability for genotype means using the average
#' variance of pairwise differences of best linear unbiased predictors (BLUPs).
#' 
#' @details The equation for Cullis heritability is as follows:
#' 
#' \eqn{H^2 = 1 - (vd_BLUP_avg / (2 * vc_g))}
#'
#' @param vd_BLUP_avg Average variance of pairwise differences
#'   among BLUPs 
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
    H2_Cullis <- 1 - (vdBLUP_avg / 2 / vc_g)
    
    return(H2_Cullis)
 }