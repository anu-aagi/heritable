
H2_Delta <- function(model, target) {

  # BLUPs for genotype main effect
  g_pred  <- asreml::predict.asreml(model,
                                    classify = target,
                                    only = target,
                                    sed = TRUE,
                                    vcov = TRUE)

  genotype_names <- levels(model$mf[[target]]) # list of genotype names
  ngeno    <- length(genotype_names)  # number of genotypes


  vc_g <- asreml::summary.asreml(model)$varcomp[target, 'component'] # varcomp of geno
  # TODO: Change below to use elements of Kinship/relationship as necessary for narrowsense
  var1 <- vc_g
  var2 <- vc_g
  cov <- 0
  Vd_g <- g_pred$sed^2 # cov
  dimnames(Vd_g) <- list(genotype_names, genotype_names)

  v <- var1 + var2 - 2*cov
  H2D_ij <- 1 - Vd_g/v
  H2D_ij
}
