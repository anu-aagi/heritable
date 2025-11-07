H2.lmer <- function(model, method = c("Cullis", "Oakey", "Piepho", "Delta", "Naive"), target = NULL) {
  method <- match.arg(method)

  #TODO: If model has not converged, warn
  

  #TODO: Check if target is in model, if not throw error

  H2 <- switch(method,
    Cullis = H2_Cullis.lmer(model, target),
    Oakey = H2_Oakey.lmer(model, target),
    Piepho = H2_Piepho.lmer(model, target),
    Delta = H2_Delta.lmer(model, target),
    Naive = H2_Naive.lmer(model, target),
    H2.default(model)
  )

  structure(H2, class = c("heritable", class(H2)))

  return(stats::setNames(H2, method))
}

H2_Naive.lmerMod <- function(model, target = NULL) {
  # Get genotype variance
  vc <- model |> lme4::VarCorr() |> as.data.frame()
  vc_g <- subset(vc, grp==target)$vcov

  # Get residual variance
  vc_e <- subset(vc, grp=="Residual")$vcov

  H2_Naive <- H2_Naive_parameters(vc_g, vc_e)

  return(H2_Naive)
}