#' @importFrom stats setNames
#' @export
H2.lmerMod <- function(model, method = c("Cullis", "Oakey", "Piepho", "Delta", "Naive"), target = NULL) {
  method <- match.arg(method)

  #TODO: If model has not converged, warn


  #TODO: Check if target is in model, if not throw error

  H2 <- switch(method,
    Cullis = H2_Cullis.lmerMod(model, target),
    Oakey = H2_Oakey.lmerMod(model, target),
    Piepho = H2_Piepho.lmerMod(model, target),
    Delta = H2_Delta.lmerMod(model, target),
    Naive = H2_Naive.lmerMod(model, target),
    H2.default(model)
  )

  structure(H2, class = c("heritable", class(H2)))

  return(stats::setNames(H2, method))
}

#' @export
H2_Naive.lmerMod <- function(model, target = NULL) {
  # TODO: We need to import lme4 if I am going to use VarCorr
  # Get genotype variance
  vc <- model |> lme4::VarCorr() |> as.data.frame()
  vc_g <- subset(vc, grp==target)$vcov

  # Get residual variance
  vc_e <- subset(vc, grp=="Residual")$vcov

  H2_Naive <- H2_Naive_parameters(vc_g, vc_e)

  return(H2_Naive)
}

#' @export
H2_Cullis.lmerMod <- function(model, target = NULL) {
  vc <- lme4::VarCorr(model)
  ngrps <- lme4::ngrps(model)
  # Note the index and kronecker order needs to be followed careful downstream
  Glist <- lapply(names(vc), function(agrp) {
    Matrix::kronecker(vc[[agrp]], diag(ngrps[[agrp]]))
  })
  G <- do.call(Matrix::bdiag, Glist)

  n <- nrow(model@frame)
  R <- diag(n) * sigma(model)^2

  X <- as.matrix(lme4::getME(model, "X"))
  Z <- as.matrix(lme4::getME(model, "Z"))

  C11 <- t(X) %*% solve(R) %*% X
  C12 <- t(X) %*% solve(R) %*% Z
  C21 <- t(Z) %*% solve(R) %*% X
  C22 <- t(Z) %*% solve(R) %*% Z + solve(G)

  C <- rbind(cbind(C11, C12),
             cbind(C21, C22))
  C_inv <- solve(C)
  gnames <- levels(model@flist[[target]])
  C22_g <- C_inv[gnames, gnames]
  n_g <- ngrps[[target]]
  vc_g <- vc[[target]][1]

  one <- matrix(1, nrow = n_g, ncol = 1)
  P_mu <- diag(n_g, n_g) - one %*% t(one)
  vdBLUP_sum <- sum(diag(P_mu %*% C22_g))
  vdBLUP_avg <- vdBLUP_sum * (2 / (n_g * (n_g - 1)))

  1 - (vdBLUP_avg / 2 / vc_g)
}
