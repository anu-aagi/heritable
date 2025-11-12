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
  R <- diag(n) * lme4::sigma(model)^2

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

#' @export
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
    model_fix <- fit_counterpart_model(model, target)
  } 
  # If fixed, refit random model
  else if(!check_target_random(model, target)){
    model_fix <- model
    model_ran <- fit_counterpart_model(model, target)
  }
  # browser()
  # Get genotype variance
  vc <- lme4::VarCorr(model_ran)
  vc_g <- vc[[target]][1]

 # Get mean variance of a difference between genotypes
  d_BLUE <- emmeans::emmeans(model_fix, specs = as.formula(paste("pairwise ~", target)))$contrasts |> as.data.frame()
  d_BLUE$var <- d_BLUE$SE^2 # Get variance

  vd_BLUE_avg <- mean(d_BLUE$var)

  # vc_g / (vc_g + vdBLUE_avg / 2)
  H2_Piepho <- H2_Piepho_parameters(vc_g, vd_BLUE_avg)
  
  return(H2_Piepho)
}

H2_Delta.lmerMod <- function(model, target = NULL) {
  # If model has not converged, warn
  check_model_convergence(model)

  # If target is not in model, error
  check_target_exists(model, target)

  # If there is more than one target, error
  check_target_single(target)

  # Check if target is random or fixed
  # TODO: How to handle if both fixed and random?
  if(check_target_both(model, target)) {
    cli::cli_abort("The target {.var {target}} is fitted as both fixed and random effect")
  } 

  browser()

  # If fixed, compute H2 delta with BLUES
  if(!check_target_random(model, target)) { 
    # Obtain BLUES
    blues <- emmeans::emmeans(model, specs = target) |> as.data.frame()
    
  }

  # Get genotype variance
  vc <- lme4::VarCorr(model)
  vc_g <- vc[[target]][1]

  # Get residual variance
  vc_e <- subset(as.data.frame(vc), grp == "Residual")$vcov

  n_g <- lme4::ngrps(model)[[target]]

  H2_Delta <- H2_Delta_parameters(n_g, vc_g, vc_e)

  return(H2_Delta)
}