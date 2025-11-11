#' @importFrom stats setNames
#' @export
H2.lmerMod <- function(model, target = NULL, method = c("Cullis", "Oakey", "Piepho", "Delta", "Naive")) {
  method <- match.arg(method)

  # If model has not converged, warn
  check_model_convergence(model)

  # If target is not in model, error
  check_target_exists(model, target)

  # If there is more than one target, error
  check_target_single(target)

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
  # If model has not converged, warn
  check_model_convergence(model)

  # If target is not in model, error
  check_target_exists(model, target)

  # If there is more than one target, error
  check_target_single(target)

  # Get genotype variance
  vc <- model |>
    lme4::VarCorr() |>
    as.data.frame()
  vc_g <- subset(vc, grp == target)$vcov

  # Get residual variance
  vc_e <- subset(vc, grp == "Residual")$vcov

  H2_Naive <- H2_Naive_parameters(vc_g, vc_e)

  return(H2_Naive)
}

#' @export
H2_Cullis.lmerMod <- function(model, target = NULL) {
  # If model has not converged, warn
  check_model_convergence(model)

  # If target is not in model, error
  check_target_exists(model, target)

  # If there is more than one target, error
  check_target_single(target)

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

  C <- rbind(
    cbind(C11, C12),
    cbind(C21, C22)
  )
  C_inv <- solve(C)
  gnames <- levels(model@flist[[target]])
  C22_g <- C_inv[gnames, gnames]
  n_g <- ngrps[[target]]
  vc_g <- vc[[target]][1]

  one <- matrix(1, nrow = n_g, ncol = 1)
  P_mu <- diag(n_g, n_g) - one %*% t(one)
  vdBLUP_sum <- sum(diag(P_mu %*% C22_g))
  vdBLUP_avg <- vdBLUP_sum * (2 / (n_g * (n_g - 1)))

  H2_Cullis_parameters(vdBLUP_avg, vc_g)
}

#' @export
H2_Oakey.lmerMod <- function(model, target = NULL) {
  # If model has not converged, warn
  check_model_convergence(model)

  # If target is not in model, error
  check_target_exists(model, target)

  # If there is more than one target, error
  check_target_single(target)

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

  C <- rbind(
    cbind(C11, C12),
    cbind(C21, C22)
  )
  C_inv <- solve(C)
  gnames <- levels(model@flist[[target]])
  C22_g <- C_inv[gnames, gnames]
  n_g <- ngrps[[target]]
  vc_g <- vc[[target]][1]

  vcov_g <- C22_g

  return(H2_Oakey_parameters(n_g, vc_g, C22_g))
}

H2_Piepho.lmerMod <- function(model, target = NULL) {
  # If model has not converged, warn
  check_model_convergence(model)

  # If target is not in model, error
  check_target_exists(model, target)

  # If there is more than one target, error
  check_target_single(target)

  # Check if target is fixed or random
  # If random, refit fixed model
  if(check_target_random(model, target)){
    model_ran <- model
    model_fixed <- fit_counterpart_model(model, target)
  } 
  # If fixed, refit random model
  else if(!check_target_random(model, target)){
    model_fixed <- model
    model_ran <- fit_counterpart_model(model, target)
  }

  # Get genotype variance
  vc <- lme4::VarCorr(model_ran)
  vc_g <- vc[[target]][1]

  # Get mean variance of a difference between genotypes
  # Get the fixed effects design matrix and coefficients
  X <- as.matrix(lme4::getME(model, "X"))
  beta <- lme4::fixef(model_fix)

  # Get predictions for each level of target
  target_levels <- levels(model_fix@frame[[target]])
  n_levels <- length(target_levels)

  # Calculate all pairwise differences
  diffs <- outer(beta[grep(target, names(beta))], 
                beta[grep(target, names(beta))], 
                "-")

  # Get variance-covariance matrix of fixed effects
  vcov_beta <- vcov(model_fix)
  vcov_target <- vcov_beta[grep(target, rownames(vcov_beta)), 
                            grep(target, colnames(vcov_beta))]

  # Calculate variance of differences: Var(β_i - β_j) = Var(β_i) + Var(β_j) - 2*Cov(β_i, β_j)
  vd_matrix <- outer(diag(vcov_target), diag(vcov_target), "+") - 
              2 * vcov_target

  # Average variance of differences
  vdBLUE_avg <- mean(vd_matrix[upper.tri(vd_matrix)])

  # vc_g / (vc_g + vdBLUE_avg / 2)
  H2_Piepho <- H2_Piepho_parameters(vc_g, vdBLUE_avg)
  
  return(H2_Piepho)
}
