#' @noRd
#' @export
h2_Standard.lmerMod <- function(model,
                               target,
                               options = NULL,
                               marginal = TRUE,
                               stratification = NULL,
                               vc = NULL,
                               ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "narrow_sense")
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  if(!is.null(vc)) stratification <- vc$stratification

  if(!is.null(stratification)){
    cli::cli_warn("Stratified heritability is not defined for the standard method, retuning {.value {NA}}.")
    return(NA)
  }

  # Get genotype variance
  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = FALSE, calc_V = TRUE, calc_C11 = FALSE,
                   marginal = marginal, stratification = stratification, ...)
  }
  V <- vc$V
  Z <- vc$Z
  G <- vc$G
  idx <- vc$idx

  if(vc$marginal || !is.null(vc$stratification)){
    X <- vc$X
    y <- vc$y
    m <- vc$m

    X_tilde <- cbind(X, Z[, idx, drop=FALSE])
    C <- ginv_sym_sparse(crossprod(X_tilde)) %*% t(X_tilde)
    C <- C[-seq_len(ncol(X)),]
    C <- crossprod(C, m)
    W <-  t(C) %*% Z[,idx]
    G_g <- W %*% G[idx,idx,drop=FALSE] %*% t(W)

    h2_Standard <-  h2_Standard_parameters(G_g, V, C)

    # Check estimability
    # P <- t(X_tilde) %*% ginv_sym_sparse(tcrossprod(X_tilde)) %*% X_tilde
    # c <- c(rep(0, ), m[,1] - m[,3])
    # plot(c - P %*% c)

  } else {
    g <- stats::model.frame(model)[[target]]
    gnames <- levels(g)
    Z_g <- Matrix::sparse.model.matrix(~ 0 + g)
    C <- Z_g %*% Matrix::Diagonal(x = 1 / as.numeric(Matrix::colSums(Z_g)))
    W <- t(C) %*% Z[,idx] # C^T Z G Z C
    G_g <- W %*% G[idx,idx,drop=FALSE] %*% t(W)

    h2_Standard <- h2_Standard_parameters(G_g, V, C)

  }

  return(h2_Standard)
}


#' @noRd
#' @export
h2_Oakey.lmerMod <- function(model,
                            target,
                            options = NULL,
                            marginal = TRUE,
                            stratification = NULL,
                            vc = NULL,
                            ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "narrow_sense")
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = TRUE, calc_V = FALSE, calc_C11 = FALSE,
                   marginal = marginal, stratification = stratification, ...)
  }
  G_g_inv <- ginv_sym_sparse(vc$G_g)

  # Parametrise is the same for broad sense and narrow sense.
  return(H2_Oakey_parameters(G_g_inv, vc$C22_g))
}

#' @noRd
#' @export
h2_Delta_pairwise.lmerMod <- function(model,
                                     target,
                                     type = c("BLUP", "BLUE"),
                                     options = NULL,
                                     marginal = TRUE,
                                     stratification = NULL,
                                     vc = NULL,
                                     ...) {
  initial_checks(model, target, options)
  type <- match.arg(type)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "narrow_sense")
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  # Check if target is random or fixed
  if (type == "BLUE") {
    h2_Delta <- h2_Delta_BLUE_pairwise.lmerMod(model, target, options,
                                              marginal, stratification, vc, ...)
  } else if (type == "BLUP") {
    h2_Delta <- h2_Delta_BLUP_pairwise.lmerMod(model, target, options,
                                              marginal, stratification, vc, ...)
  }

  return(h2_Delta)
}

#' @keywords internal
h2_Delta_BLUP_pairwise.lmerMod<- function(model,
                                         target,
                                         options = NULL,
                                         marginal = TRUE,
                                         stratification = NULL,
                                         vc = NULL,
                                         ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "narrow_sense")
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = TRUE, calc_V = FALSE, calc_C11 = FALSE,
                   marginal = marginal, stratification = stratification, ...)
  }
  G_g <- vc$G_g
  C22_g <- vc$C22_g

  # Compute variance of difference from G
  delta_g <- var_diff(G_g)

  # Compute variance of difference from PEV
  delta_pev <- var_diff(C22_g)
  diag(delta_pev) <- NA

  # H2 Delta BLUP, same parameterisation as the broad sense.
  h2_Delta_BLUP <- H2_Delta_parameters(delta_g, delta_pev, "BLUP")

  dimnames(h2_Delta_BLUP) <- dimnames(delta_g)

  h2_Delta_BLUP
}

#' @keywords internal
h2_Delta_BLUE_pairwise.lmerMod<- function(model,
                                         target,
                                         options = NULL,
                                         marginal = TRUE,
                                         stratification = NULL,
                                         vc = NULL, ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "narrow_sense")
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  # Get genotype variance
  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = FALSE, calc_V = FALSE, calc_C11 = TRUE,
                   marginal = marginal, stratification = stratification, ...)
  }
  gnames <- vc$gnames

  delta <- var_diff(vc$C11_g)
  dimnames(delta) <- list(gnames, gnames)
  diag(delta) <- NA

  delta_g <- var_diff(vc$G_g)
  h2_Delta_BLUE <- H2_Delta_parameters(delta_g, delta, "BLUE")

  # # Extract vc_g and vc_e
  # if(is.null(vc)){
  #   vc <- var_comp(model, target, calc_C22 = FALSE, calc_V = FALSE,
  #                  marginal, stratification, ...)
  # }
  # if(vc$marginal || !is.null(vc$stratification)){
  #   cli::cli_warn("Delta (BLUE) heritability can only be computed on main genetic effects, retuning {.value {NA}}.")
  #   return(NA)
  # }
  # G_g <- vc$G_g
  # gnames <- vc$gnames
  # n_g <- vc$n_g
  #
  # # Get delta
  # delta_g <- var_diff(G_g)
  #
  # # Calculate mean variance of a difference between genotypes
  # conterpart <- fit_counterpart_model(model, target)
  #
  # # Take pairwise differences and turn into variance-covariance matrix
  # frm <- as.formula(paste("pairwise ~", target))
  # EMM_fit <- emmeans::emmeans(conterpart, specs = frm)$contrasts
  # EMM_fit <- data.frame(EMM_fit) # Get variance
  # EMM_fit$var <- EMM_fit$SE^2
  #
  # # Start with empty variance matrix for differences
  # delta <- matrix(0, nrow = n_g, ncol = n_g)
  # dimnames(delta) <- list(gnames, gnames)
  #
  # # Fill in the pairwise variance of differences
  # for (i in 1:nrow(EMM_fit)) {
  #   # Extract genotype names from contrast column
  #   pair <- strsplit(as.character(EMM_fit$contrast[i]), " - ")[[1]]
  #   g1 <- pair[1]
  #   g2 <- pair[2]
  #
  #   delta[g1, g2] <- EMM_fit$var[i]
  #   delta[g2, g1] <- EMM_fit$var[i] # symmetric
  # }
  #
  # delta <- Matrix::Matrix(delta)
  # diag(delta) <- NA
  #
  # # H2 Delta BLUE
  # h2_Delta_BLUE <- H2_Delta_parameters(delta_g, delta, "BLUE")
  #
  # dimnames(h2_Delta_BLUE) <- dimnames(delta_g)
  #
  # h2_Delta_BLUE

  dimnames(h2_Delta_BLUE) <- dimnames(delta_g)

  h2_Delta_BLUE
}

#' @noRd
#' @export
h2_Cullis.lmerMod <- function(model,
                              target = NULL,
                              options = NULL,
                              marginal = TRUE,
                              stratification = NULL,
                              vc = NULL,
                              ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense")
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = TRUE, calc_V = FALSE, calc_C11 = FALSE,
                   marginal = marginal, stratification = stratification, ...)
  }
  s2_g <- mean(diag(vc$G_g))
  n <- vc$n_g
  C22_g <- vc$C22_g

  # This is equivalent to delta <- var_diff(C22_g); delta_avg = mean(delta[lower.tri(delta)])
  delta_avg <- (2 / (n * (n - 1))) * (n * sum(diag(C22_g)) - sum(C22_g))

  return(H2_Cullis_parameters(delta_avg, s2_g))
}


#' @noRd
#' @export
h2_Piepho.lmerMod <- function(model,
                              target = NULL,
                              options = NULL,
                              marginal = TRUE,
                              stratification = NULL,
                              vc = NULL,
                              ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense")
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  # Get genotype variance
  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = FALSE, calc_V = FALSE, calc_C11 = TRUE,
                   marginal = marginal, stratification = stratification, ...)
  }
  gnames <- vc$gnames

  delta <- var_diff(vc$C11_g)
  dimnames(delta) <- list(gnames, gnames)
  diag(delta) <- NA
  delta_avg <- mean(delta[upper.tri(delta)])

  # s2_g / (s2_g + delta_avg / 2)
  s2_g <- mean(diag(vc$G_g))
  H2_Piepho <- H2_Piepho_parameters(s2_g, delta_avg)

  ## Conterpart method

  # Get genotype variance
  # if(is.null(vc)){
  #   vc <- var_comp(model, target, calc_C22 = FALSE, calc_V = FALSE, marginal, stratification, ...)
  # }
  # if(vc$marginal || !is.null(vc$stratification)){
  #   cli::cli_warn("Piepho heritability can only be computed on main genetic effects, retuning {.value {NA}}.")
  #   return(NA)
  # }
  # G_g <- vc$G_g
  # s2_g <- mean(diag(G_g))
  #
  # conterpart <- fit_counterpart_model(model, target)
  #
  # # Get mean variance of a difference between genotypes
  # frm <- as.formula(paste("pairwise ~", target))
  # EMM_fit <- emmeans::emmeans(conterpart, specs = frm)$contrasts
  # delta_avg <- mean(data.frame(EMM_fit)$SE^2) # Get variance
  #
  # # s2_g / (s2_g + delta_avg / 2)
  # H2_Piepho <- H2_Piepho_parameters(s2_g, delta_avg)

  return(H2_Piepho)
}


#' @noRd
#' @export
H2_Standard.lmerMod <- function(model,
                                target = NULL,
                                options = NULL,
                                marginal = TRUE,
                                stratification = NULL,
                                vc = NULL,
                                ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense")
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  # Check if all random terms contain the target.
  trms <- pull_terms_without_specials(model)$random

  # Only use the simple definition when there is only main genetic effect.
  if(all(trms == target)){
    # Get genotype variance
    if(is.null(vc)){
      G_g <- var_comp(model, target, calc_C22 = FALSE, calc_V = FALSE, calc_C11 = FALSE,
                      marginal = marginal, stratification = stratification, ...)$G_g
    } else {
      G_g <- vc$G_g
    }
    s2_g <- mean(diag(G_g))

    # Get residual variance
    s2_eps <- stats::sigma(model)^2

    n_r <- table(model@flist[[target]])

    H2_Standard <- H2_Standard_parameters(s2_g, s2_eps, n_r)
  } else {

    # Get genotype variance
    if(is.null(vc)){
      vc <- var_comp(model, target, calc_C22 = FALSE, calc_V = TRUE, calc_C11 = FALSE,
                     marginal = marginal, stratification = stratification, ...)
    }

    V <- vc$V
    Z <- vc$Z
    G <- vc$G
    idx <- vc$idx

    if(vc$marginal || !is.null(vc$stratification)){
      X <- vc$X
      m <- vc$m

      X_tilde <- cbind(X, Z[, idx, drop=FALSE])
      C <- ginv_sym_sparse(crossprod(X_tilde)) %*% t(X_tilde)
      C <- C[-seq_len(ncol(X)),]
      C <- crossprod(C, m)
      W <-  t(C) %*% Z[,idx]
      G_g <- W %*% G[idx,idx,drop=FALSE] %*% t(W)

      H2_Standard <-  h2_Standard_parameters(G_g, V, C)

      # Check estimability
      # P <- t(X_tilde) %*% ginv_sym_sparse(tcrossprod(X_tilde)) %*% X_tilde
      # c <- c(rep(0, ), m[,1] - m[,3])
      # plot(c - P %*% c)

      # Validate stratified heritability using the following model
      # model <-  lme4::lmer(y ~  rep + (1 | gen * loc), data = lettuce_phenotypes)
      #
      # g <- stats::model.frame(model)[[target]]
      # gnames <- levels(g)
      # Z_g <- Matrix::sparse.model.matrix(~ 0 + g)
      # Z_g <- Z_g * (stats::model.frame(model)[["loc"]] == "L3")
      # C <- Z_g %*% Diagonal(x = 1 / as.numeric(Matrix::colSums(Z_g)))
      # C[is.na(C)] <- 0
      # W <- t(C) %*% Z[,idx] # C^T Z G Z C
      # G_g <- W %*% G[idx,idx,drop=FALSE] %*% t(W)
      # S <- crossprod(C, V %*% C)
      # delta_y <- var_diff(S)
      # delta_g <- var_diff(G_g)
      # 1/mean((delta_y/delta_g)[lower.tri(delta_y)], na.rm=TRUE)
      #
      # this yields the same result as using the above code.


    } else {
      g <- stats::model.frame(model)[[target]]
      gnames <- levels(g)
      Z_g <- Matrix::sparse.model.matrix(~ 0 + g)
      C <- Z_g %*% Matrix::Diagonal(x = 1 / as.numeric(Matrix::colSums(Z_g)))
      W <- t(C) %*% Z[,idx] # C^T Z G Z C
      G_g <- W %*% G[idx,idx,drop=FALSE] %*% t(W)

      H2_Standard <- h2_Standard_parameters(G_g, V, C)

    }

  }

  return(H2_Standard)
}


#' @noRd
#' @export
H2_Cullis.lmerMod <- function(model,
                              target = NULL,
                              options = NULL,
                              marginal = TRUE,
                              stratification = NULL,
                              vc = NULL,
                              ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense")
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = TRUE, calc_V = FALSE, calc_C11 = FALSE,
                   marginal = marginal, stratification = stratification, ...)
  }
  s2_g <- mean(diag(vc$G_g))
  n <- vc$n_g
  C22_g <- vc$C22_g

  # This is equivalent to delta <- var_diff(C22_g); delta_avg = mean(delta[lower.tri(delta)])
  delta_avg <- (2 / (n * (n - 1))) * (n * sum(diag(C22_g)) - sum(C22_g))

  return(H2_Cullis_parameters(delta_avg, s2_g))
}

#' @noRd
#' @export
H2_Oakey.lmerMod <- function(model,
                             target = NULL,
                             options = NULL,
                             marginal = TRUE,
                             stratification = NULL,
                             vc = NULL,
                             ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense")
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = TRUE, calc_V = FALSE, calc_C11 = FALSE,
                   marginal = marginal, stratification = stratification, ...)
  }
  G_g_inv <- ginv_sym_sparse(vc$G_g)

  return(H2_Oakey_parameters(G_g_inv, vc$C22_g))
}

#' @noRd
#' @export
H2_Piepho.lmerMod <- function(model,
                              target = NULL,
                              options = NULL,
                              marginal = TRUE,
                              stratification = NULL,
                              vc = NULL,
                              ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense")
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  # Get genotype variance
  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = FALSE, calc_V = FALSE, calc_C11 = TRUE,
                   marginal = marginal, stratification = stratification, ...)
  }
  gnames <- vc$gnames

  delta <- var_diff(vc$C11_g)
  dimnames(delta) <- list(gnames, gnames)
  diag(delta) <- NA
  delta_avg <- mean(delta[upper.tri(delta)])

  # s2_g / (s2_g + delta_avg / 2)
  s2_g <- mean(diag(vc$G_g))
  H2_Piepho <- H2_Piepho_parameters(s2_g, delta_avg)

  ## Conterpart method

  # Get genotype variance
  # if(is.null(vc)){
  #   vc <- var_comp(model, target, calc_C22 = FALSE, calc_V = FALSE, marginal, stratification, ...)
  # }
  # if(vc$marginal || !is.null(vc$stratification)){
  #   cli::cli_warn("Piepho heritability can only be computed on main genetic effects, retuning {.value {NA}}.")
  #   return(NA)
  # }
  # G_g <- vc$G_g
  # s2_g <- mean(diag(G_g))
  #
  # conterpart <- fit_counterpart_model(model, target)
  #
  # # Get mean variance of a difference between genotypes
  # frm <- as.formula(paste("pairwise ~", target))
  # EMM_fit <- emmeans::emmeans(conterpart, specs = frm)$contrasts
  # delta_avg <- mean(data.frame(EMM_fit)$SE^2) # Get variance
  #
  # # s2_g / (s2_g + delta_avg / 2)
  # H2_Piepho <- H2_Piepho_parameters(s2_g, delta_avg)

  return(H2_Piepho)
}

#' @noRd
#' @export
H2_Delta_pairwise.lmerMod <- function(model,
                                      target = NULL,
                                      type = c("BLUP", "BLUE"),
                                      options = NULL,
                                      marginal = TRUE,
                                      stratification = NULL,
                                      vc = NULL,
                                      ...) {
  initial_checks(model, target, options)
  type <- match.arg(type)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense")
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  # Check if target is random or fixed
  if (type == "BLUE") {
    H2_Delta <- H2_Delta_BLUE_pairwise.lmerMod(model, target, options,
                                               marginal, stratification, vc, ...)
  } else if (type == "BLUP") {
    H2_Delta <- H2_Delta_BLUP_pairwise.lmerMod(model, target, options,
                                               marginal, stratification, vc, ...)
  }

  return(H2_Delta)
}

#' @keywords internal
H2_Delta_BLUE_pairwise.lmerMod <- function(model,
                                           target = NULL,
                                           options = NULL,
                                           marginal = TRUE,
                                           stratification = NULL,
                                           vc = NULL, ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense")
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  # Get genotype variance
  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = FALSE, calc_V = FALSE, calc_C11 = TRUE,
                   marginal = marginal, stratification = stratification, ...)
  }
  gnames <- vc$gnames

  delta <- var_diff(vc$C11_g)
  dimnames(delta) <- list(gnames, gnames)
  diag(delta) <- NA

  s2_g <- mean(diag(vc$G_g))
  H2_Delta_BLUE <- H2_Delta_parameters(2*s2_g, delta, "BLUE")

  # Counterpart approach
  #
  # # Extract vc_g and vc_e
  # if(is.null(vc)){
  #   vc <- var_comp(model, target, calc_C22 = FALSE, calc_V = FALSE, marginal, stratification, ...)
  # }
  # if(vc$marginal || !is.null(vc$stratification)){
  #   cli::cli_warn("Delta (BLUE) heritability can only be computed on main genetic effects, retuning {.value {NA}}.")
  #   return(NA)
  # }
  # s2_g <- mean(diag(vc$G_g))
  #
  # conterpart <- fit_counterpart_model(model, target)
  #
  # # Calculate mean variance of a difference between genotypes
  # frm <- as.formula(paste("pairwise ~", target))
  # EMM_fit <- emmeans::emmeans(conterpart, specs = frm)$contrasts
  # EMM_fit <- data.frame(EMM_fit) # Get variance
  # EMM_fit$var <- EMM_fit$SE^2
  #
  # # Take pairwise differences and turn into variance-covariance matrix
  # gnames <- vc$gnames
  # n_g <- vc$n_g
  #
  # # Start with empty variance matrix for differences
  # delta <- matrix(0, nrow = n_g, ncol = n_g)
  # dimnames(delta) <- list(gnames, gnames)
  #
  # # Fill in the pairwise variance of differences
  # for (i in 1:nrow(EMM_fit)) {
  #   # Extract genotype names from contrast column
  #   pair <- strsplit(as.character(EMM_fit$contrast[i]), " - ")[[1]]
  #   g1 <- pair[1]
  #   g2 <- pair[2]
  #
  #   delta[g1, g2] <- EMM_fit$var[i]
  #   delta[g2, g1] <- EMM_fit$var[i] # symmetric
  # }
  #
  # delta <- Matrix::Matrix(delta)
  # diag(delta) <- NA

  # H2 Delta BLUE
  # H2_Delta_BLUE <- H2_Delta_parameters(2*s2_g, delta, "BLUE")

  dimnames(H2_Delta_BLUE) <- dimnames(delta)

  return(H2_Delta_BLUE)
}

#' @keywords internal
H2_Delta_BLUP_pairwise.lmerMod <- function(model,
                                           target = NULL,
                                           options = NULL,
                                           marginal = TRUE,
                                           stratification = NULL,
                                           vc = NULL, ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense")
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = TRUE, calc_V = FALSE, calc_C11 = FALSE,
                   marginal = marginal, stratification = stratification, ...)
  }
  s2_g <- mean(diag(vc$G_g))
  C22_g <- vc$C22_g

  # Compute variance of difference from PEV
  delta <- var_diff(C22_g)
  diag(delta) <- NA
  dimnames(delta) <- list(vc$gnames, vc$gnames)

  # H2 Delta BLUP
  H2_Delta_BLUP <- H2_Delta_parameters(2*s2_g, delta, "BLUP")

  dimnames(H2_Delta_BLUP) <- dimnames(delta)

  H2_Delta_BLUP
}
