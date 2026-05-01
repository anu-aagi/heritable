#' Calculate standard heritability from asreml model
#' @export
#' @noRd
#' @return Numeric
#' @examples
#' # asreml model (Requires license)
#' \dontrun{
#' lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
#'                                  random = ~ vm(gen, lettuce_GRM),
#'                                  data = lettuce_subset,
#'                                  trace = FALSE
#'                                  )
#'
#' h2_Standard.asreml(lettuce_asreml, target = "gen", source = list(lettuce_GRM = lettuce_GRM))
#' }
h2_Standard.asreml <- function(model,
                               target,
                               options = NULL,
                               marginal = TRUE,
                               stratification = NULL,
                               vc = NULL,
                               source = list(),
                               ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "narrow_sense", source)
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  model <- check_deisgn_exsits(model, build_design = FALSE, source = source)
  mf <- model$mf

  # Get genotype variance
  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = FALSE, calc_V = TRUE, calc_C11 = FALSE,
                   marginal = marginal, stratification = stratification, source = source, ...)
  }
  V <- vc$V
  Z <- vc$Z
  G <- vc$G
  idx <- vc$idx


  if(vc$marginal || !is.null(vc$stratification)){
    X <- vc$X
    W <- vc$W

    active_g <- attr(W, "active")

    X_tilde <- cbind(X, Z[, idx, drop=FALSE])
    C <- ginv_sym_sparse(crossprod(X_tilde)) %*% t(X_tilde)
    C <- C[-seq_len(ncol(X)),]
    C <- crossprod(C, W)
    W <-  t(C) %*% Z[,idx]
    G_g <- W %*% G[idx,idx,drop=FALSE] %*% t(W)

    h2_Standard <-  h2_Standard_parameters(G_g, V, C, active_g)

    # Check estimability
    # P <- t(X_tilde) %*% ginv_sym_sparse(tcrossprod(X_tilde)) %*% X_tilde
    # c <- c(rep(0, ), W[,1] - W[,3])
    # plot(c - P %*% c)

  } else {
    g <- mf[[target]]
    gnames <- levels(g)
    Z_g <- Matrix::sparse.model.matrix(~ 0 + g)
    C <- Z_g %*% Matrix::Diagonal(x = 1 / as.numeric(Matrix::colSums(Z_g)))
    W <- t(C) %*% Z[,idx]
    G_g <- W %*% G[idx,idx,drop=FALSE] %*% t(W)

    h2_Standard <- h2_Standard_parameters(G_g, V, C, NULL)
  }

  return(h2_Standard)
}


#' Calculate Oakey's heritability from asreml model
#' @export
#' @noRd
#' @return Numeric
#' @examples
#' # asreml model (Requires license)
#' \dontrun{
#' lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
#'                                  random = ~ vm(gen, lettuce_GRM),
#'                                  data = lettuce_subset,
#'                                  trace = FALSE
#'                                  )
#'
#' h2_Oakey.asreml(lettuce_asreml, target = "gen", source = list(lettuce_GRM = lettuce_GRM))
#' }
h2_Oakey.asreml <- function(model,
                            target,
                            options = NULL,
                            marginal = TRUE,
                            stratification = NULL,
                            vc = NULL,
                            source = list(),
                            ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "narrow_sense", source)
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = TRUE, calc_V = FALSE, calc_C11 = FALSE,
                   marginal = marginal, stratification = stratification, source = source, ...)

  }
  G_g_inv <- ginv_sym_sparse(vc$G_g)

  # Parametrise is the same for broad sense and narrow sense.
  return(H2_Oakey_parameters(G_g_inv, vc$C22_g))
}

#' Calculate pairwise heritability from asreml model
#' @export
#' @noRd
#' @return Numeric
#' @examples
#' \dontrun{
#' lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
#'                                  random = ~ vm(gen, lettuce_GRM),
#'                                  data = lettuce_subset,
#'                                  trace = FALSE
#'                                  )
#'
#' h2_Delta_pairwise.asreml(lettuce_asreml, target = "gen", type = "BLUP", source = list(lettuce_GRM = lettuce_GRM))
#' }
h2_Delta_pairwise.asreml <- function(model,
                                     target,
                                     type = c("BLUP", "BLUE"),
                                     options = NULL,
                                     marginal = TRUE,
                                     stratification = NULL,
                                     vc = NULL,
                                     source = list(),
                                     ...) {
  initial_checks(model, target, options)
  type <- match.arg(type)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "narrow_sense", source)
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  # Check if target is random or fixed
  if (type == "BLUE") {
    h2_Delta <- h2_Delta_BLUE_pairwise.asreml(model, target, options,
                                              marginal, stratification, vc, source, ...)
  } else if (type == "BLUP") {
    h2_Delta <- h2_Delta_BLUP_pairwise.asreml(model, target, options,
                                              marginal, stratification, vc, source, ...)
  }

  return(h2_Delta)
}

#' @keywords internal
h2_Delta_BLUP_pairwise.asreml<- function(model,
                                         target,
                                         options = NULL,
                                         marginal = TRUE,
                                         stratification = NULL,
                                         vc = NULL,
                                         source = list(),
                                         ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "narrow_sense", source)
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = TRUE, calc_V = FALSE, calc_C11 = FALSE,
                   marginal = marginal, stratification = stratification, source = source, ...)
  }
  gnames <- vc$gnames
  G_g <- vc$G_g
  C22_g <- vc$C22_g

  # Compute variance of difference from G
  delta_g <- var_diff(G_g)

  # Compute variance of difference from PEV
  delta_pev <- var_diff(C22_g)

  # H2 Delta BLUP, same parameterisation as the broad sense.
  h2_Delta_BLUP <- H2_Delta_parameters(delta_g, delta_pev, "BLUP")

  dimnames(h2_Delta_BLUP) <- list(gnames, gnames)
  diag(h2_Delta_BLUP) <- NA

  attr(h2_Delta_BLUP, "delta_g") <- delta_g
  attr(h2_Delta_BLUP, "delta_pev") <- delta_pev

  h2_Delta_BLUP

}

#' @keywords internal
h2_Delta_BLUE_pairwise.asreml<- function(model,
                                         target,
                                         options = NULL,
                                         marginal = TRUE,
                                         stratification = NULL,
                                         vc = NULL,
                                         source = list(),
                                         ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "narrow_sense", source)
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  # Get genotype variance
  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = FALSE, calc_V = FALSE, calc_C11 = TRUE,
                   marginal = marginal, stratification = stratification, source = source, ...)
  }
  gnames <- vc$gnames

  delta <- var_diff(vc$C11_g)

  delta_g <- var_diff(vc$G_g_tilde)

  h2_Delta_BLUE <- H2_Delta_parameters(delta_g, delta, "BLUE")

  dimnames(h2_Delta_BLUE) <- list(gnames, gnames)
  diag(h2_Delta_BLUE) <- NA

  attr(h2_Delta_BLUE, "delta_g") <- delta_g
  attr(h2_Delta_BLUE, "delta_pev") <- delta

  h2_Delta_BLUE

  # # Extract vc_g and vc_e
  # if(is.null(vc)){
  #   vc <- var_comp(model, target, calc_C22 = FALSE, calc_V = FALSE, marginal, stratification, source = source, ...)
  # }
  # if(vc$marginal || !is.null(vc$stratification)){
  #   cli::cli_warn("Delta (BLUE) heritability can only be computed on main genetic effects, retuning {.value {NA}}.")
  #   return(NA)
  # }
  # G_g <- vc$G_g
  # gnames <- vc$gnames
  #
  # conterpart <- fit_counterpart_model(model, target)
  #
  # # Get delta
  # delta_g <- var_diff(G_g)
  #
  # delta <- predict(conterpart, classify = target, sed = TRUE, trace = FALSE)$sed^2
  # diag(delta) <- NA
  #
  # # H2 Delta BLUE
  # h2_Delta_BLUE <- H2_Delta_parameters(delta_g, delta, "BLUE")
}

#' Calculate Cullis's heritability from asreml model
#' @export
#' @noRd
#' @importFrom stats predict
#' @return Numeric
#' @examples
#' # asreml model (Requires license)
#' \dontrun{
#' lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
#'                                  random = ~ vm(gen, lettuce_GRM),
#'                                  data = lettuce_subset,
#'                                  trace = FALSE
#'                                  )
#'
#' h2_Cullis.asreml(lettuce_asreml, target = "gen", source = list(lettuce_GRM = lettuce_GRM))
#' }
h2_Cullis.asreml <- function(model,
                             target = NULL,
                             options = NULL,
                             marginal = TRUE,
                             stratification = NULL,
                             vc = NULL,
                             source = list(),
                             ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "narrow_sense", source)
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = TRUE, calc_V = FALSE, calc_C11 = FALSE,
                   marginal = marginal, stratification = stratification, source = source, ...)
  }
  s2_g <- mean(diag(vc$G_g))
  n <- vc$n_g
  C22_g <- vc$C22_g

  # This is equivalent to delta <- var_diff(C22_g); delta_avg = mean(delta[lower.tri(delta)])
  delta_avg <- (2 / (n * (n - 1))) * (n * sum(diag(C22_g)) - sum(C22_g))

  return(H2_Cullis_parameters(delta_avg, s2_g))
}


#' Calculate Piepho's heritability from asreml model
#' @export
#' @noRd
#' @return Numeric
#' @examples
#' # asreml model (Requires license)
#' \dontrun{
#' lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
#'                                  random = ~ vm(gen, lettuce_GRM),
#'                                  data = lettuce_subset,
#'                                  trace = FALSE
#'                                  )
#'
#' h2_Piepho.asreml(lettuce_asreml, target = "gen", source = list(lettuce_GRM = lettuce_GRM))
#' }
h2_Piepho.asreml <- function(model,
                             target = NULL,
                             options = NULL,
                             marginal = TRUE,
                             stratification = NULL,
                             vc = NULL,
                             source = list(),
                             ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "narrow_sense", source)
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  # Get genotype variance
  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = FALSE, calc_V = FALSE, calc_C11 = TRUE,
                   marginal = marginal, stratification = stratification, source = source, ...)
  }
  gnames <- vc$gnames

  delta <- var_diff(vc$C11_g)
  dimnames(delta) <- list(gnames, gnames)
  diag(delta) <- NA
  delta_avg <- mean(delta[upper.tri(delta)])

  # s2_g / (s2_g + delta_avg / 2)
  G_g_tilde <- vc$G_g_tilde_no_cov
  s2_g_tilde <- mean(var_diff(G_g_tilde)[upper.tri(G_g_tilde)]) / 2
  return(H2_Piepho_parameters(s2_g_tilde, delta_avg))

  # # Get genotype variance
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
  # delta <- predict(conterpart,
  #   classify = target,
  #   sed = TRUE,
  #   trace = FALSE
  # )$sed^2
  #
  # delta_avg <- mean(delta[upper.tri(delta, diag = FALSE)])
  #
  # # Calculate Piepho's H2
  # H2_Piepho <- H2_Piepho_parameters(s2_g, delta_avg)
}

#' Calculate Cullis's heritability from asreml model
#' @export
#' @noRd
#' @importFrom stats predict
#' @return Numeric
#' @examples
#' # asreml model (Requires license)
#' \dontrun{
#' lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
#'                                  random = ~ gen,
#'                                  data = lettuce_subset,
#'                                  trace = FALSE
#'                                  )
#'
#' H2_Cullis.asreml(lettuce_asreml, target = "gen")
#' }
H2_Cullis.asreml <- function(model,
                              target = NULL,
                              options = NULL,
                              marginal = TRUE,
                              stratification = NULL,
                              vc = NULL,
                              source = list(),
                              ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense", source)
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = TRUE, calc_V = FALSE, calc_C11 = FALSE,
                   marginal = marginal, stratification = stratification, source = source, ...)
  }
  s2_g <- mean(diag(vc$G_g))
  n <- vc$n_g
  C22_g <- vc$C22_g

  # This is equivalent to delta <- var_diff(C22_g); delta_avg = mean(delta[lower.tri(delta)])
  delta_avg <- (2 / (n * (n - 1))) * (n * sum(diag(C22_g)) - sum(C22_g))

  return(H2_Cullis_parameters(delta_avg, s2_g))
}

#' Calculate Oakey's heritability from asreml model
#' @export
#' @noRd
#' @return Numeric
#' @examples
#' # asreml model (Requires license)
#' \dontrun{
#' lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
#'                                  random = ~ gen,
#'                                  data = lettuce_subset,
#'                                  trace = FALSE
#'                                  )
#'
#' H2_Oakey.asreml(lettuce_asreml, target = "gen")
#' }
H2_Oakey.asreml <- function(model,
                             target = NULL,
                             options = NULL,
                             marginal = TRUE,
                             stratification = NULL,
                             vc = NULL,
                             source = list(),
                             ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense", source)
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = TRUE, calc_V = FALSE, calc_C11 = FALSE,
                   marginal = marginal, stratification = stratification, source = source, ...)
  }
  G_g_inv <- ginv_sym_sparse(vc$G_g)

  return(H2_Oakey_parameters(G_g_inv, vc$C22_g))
}

#' Calculate Piepho's heritability from asreml model
#' @export
#' @noRd
#' @return Numeric
#' @examples
#' # asreml model (Requires license)
#' \dontrun{
#' lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
#'                                  random = ~ gen,
#'                                  data = lettuce_subset,
#'                                  trace = FALSE
#'                                  )
#'
#' H2_Piepho.asreml(lettuce_asreml, target = "gen")
#' }
H2_Piepho.asreml <- function(model,
                              target = NULL,
                              options = NULL,
                              marginal = TRUE,
                              stratification = NULL,
                              vc = NULL,
                              source = list(),
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
                   marginal = marginal, stratification = stratification, source = source, ...)
  }
  gnames <- vc$gnames

  delta <- var_diff(vc$C11_g)
  dimnames(delta) <- list(gnames, gnames)
  diag(delta) <- NA
  delta_avg <- mean(delta[upper.tri(delta)])

  # s2_g / (s2_g + delta_avg / 2)
  G_g_tilde <- vc$G_g_tilde_no_cov
  s2_g_tilde <- mean(var_diff(G_g_tilde)[upper.tri(G_g_tilde)]) / 2
  return(H2_Piepho_parameters(s2_g_tilde, delta_avg))

  # # Get genotype variance
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
  # delta <- predict(conterpart,
  #   classify = target,
  #   sed = TRUE,
  #   trace = FALSE
  # )$sed^2
  #
  # delta_avg <- mean(delta[upper.tri(delta, diag = FALSE)])
  #
  # # Calculate Piepho's H2
  # H2_Piepho <- H2_Piepho_parameters(s2_g, delta_avg)
}

#' Calculate pairwise heritability from asreml model
#' @export
#' @noRd
#' @return Numeric
#' @examples
#' \dontrun{
#' lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
#'                                  random = ~ gen,
#'                                  data = lettuce_subset,
#'                                  trace = FALSE
#'                                  )
#'
#' H2_Delta_pairwise.asreml(lettuce_asreml, target = "gen", type = "BLUP")
#' }
H2_Delta_pairwise.asreml <- function(model,
                                      target = NULL,
                                      type = c("BLUP", "BLUE"),
                                      options = NULL,
                                      marginal = TRUE,
                                      stratification = NULL,
                                      vc = NULL,
                                      source = list(),
                                     ...) {
  initial_checks(model, target, options)
  type <- match.arg(type)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense", source)
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  # Check if target is random or fixed
  if (type == "BLUE") {
    H2_Delta <- H2_Delta_BLUE_pairwise.asreml(model, target, options,
                                               marginal, stratification, vc, source, ...)
  } else if (type == "BLUP") {
    H2_Delta <- H2_Delta_BLUP_pairwise.asreml(model, target, options,
                                               marginal, stratification, vc, source, ...)
  }

  return(H2_Delta)
}

#' @keywords internal
H2_Delta_BLUP_pairwise.asreml<- function(model,
                                           target = NULL,
                                           options = NULL,
                                           marginal = TRUE,
                                           stratification = NULL,
                                           vc = NULL,
                                           source = list(),
                                         ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense", source)
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = TRUE, calc_V = FALSE, calc_C11 = FALSE,
                   marginal = marginal, stratification = stratification, source = source, ...)
  }
  gnames <- vc$gnames
  s2_g <- mean(diag(vc$G_g))
  C22_g <- vc$C22_g

  # Compute variance of difference from PEV
  delta <- var_diff(C22_g)

  # H2 Delta BLUP
  H2_Delta_BLUP <- H2_Delta_parameters(2*s2_g, delta, "BLUP")

  dimnames(H2_Delta_BLUP) <- list(vc$gnames, vc$gnames)
  diag(H2_Delta_BLUP) <- NA

  attr(H2_Delta_BLUP, "delta_g") <- 2*s2_g
  attr(H2_Delta_BLUP, "delta_pev") <- delta

  H2_Delta_BLUP

}

#' @keywords internal
H2_Delta_BLUE_pairwise.asreml<- function(model,
                                         target = NULL,
                                         options = NULL,
                                         marginal = TRUE,
                                         stratification = NULL,
                                         vc = NULL,
                                         source = list(),
                                         ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense", source)
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }


  # Get genotype variance
  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = FALSE, calc_V = FALSE, calc_C11 = TRUE,
                   marginal = marginal, stratification = stratification, source = source, ...)
  }
  gnames <- vc$gnames

  delta <- var_diff(vc$C11_g)

  G_g_tilde <- vc$G_g_tilde_no_cov
  s2_g_tilde <- mean(var_diff(G_g_tilde)[upper.tri(G_g_tilde)]) / 2

  H2_Delta_BLUE <- H2_Delta_parameters(2*s2_g_tilde, delta, "BLUE")

  dimnames(H2_Delta_BLUE) <- list(gnames, gnames)
  diag(H2_Delta_BLUE) <- NA

  attr(H2_Delta_BLUE, "delta_g") <- 2*s2_g_tilde
  attr(H2_Delta_BLUE, "delta_pev") <- delta

  H2_Delta_BLUE

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
  # # Get delta
  # delta <- predict(conterpart, classify = target, sed = TRUE, trace = FALSE)$sed^2
  # diag(delta) <- NA
  # dimnames(delta) <- list(vc$gnames, vc$gnames)
  #
  # # H2 Delta BLUE
  # H2_Delta_BLUE <- H2_Delta_parameters(2*s2_g, delta, "BLUE")
}



#' Calculate standard heritability from asreml model
#' @export
#' @noRd
#' @return Numeric
#' @examples
#' # asreml model (Requires license)
#' \dontrun{
#' lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
#'                                  random = ~ gen,
#'                                  data = lettuce_subset,
#'                                  trace = FALSE
#'                                  )
#'
#' H2_Standard.asreml(lettuce_asreml, target = "gen")
#' }
H2_Standard.asreml <- function(model,
                                target = NULL,
                                options = NULL,
                                marginal = TRUE,
                                stratification = NULL,
                                vc = NULL,
                                source = list(),
                                ...) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense", source)
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  model <- check_deisgn_exsits(model, build_design = FALSE, source = source)
  mf <- model$mf

  # Check if all random terms contain the target.
  trms <- pull_terms_without_specials(model)$random

  # Only use the simple definition when there is only main genetic effect.
  if(all(trms == target)){
    # Get genotype variance
    if(is.null(vc)){
      G_g <- var_comp(model, target, calc_C22 = FALSE, calc_V = FALSE, calc_C11 = FALSE,
                      marginal = marginal, stratification = stratification, source = source, ...)$G_g
    } else {
      G_g <- vc$G_g
    }
    s2_g <- mean(diag(G_g))

    # Get residual variance
    s2_eps <- model$sigma2

    n_r <- table(mf[[target]])

    H2_Standard <- H2_Standard_parameters(s2_g, s2_eps, n_r)
  } else {

    # Get genotype variance
    if(is.null(vc)){
      vc <- var_comp(model, target, calc_C22 = FALSE, calc_V = TRUE, calc_C11 = FALSE,
                     marginal = marginal, stratification = stratification, source = source, ...)
    }

    V <- vc$V
    Z <- vc$Z
    G <- vc$G
    idx <- vc$idx

    if(vc$marginal || !is.null(vc$stratification)){
      X <- vc$X
      W <- vc$W

      active_g <- attr(W, "active")

      X_tilde <- cbind(X, Z[, idx, drop=FALSE])
      C <- ginv_sym_sparse(crossprod(X_tilde)) %*% t(X_tilde)
      C <- C[-seq_len(ncol(X)),]
      C <- crossprod(C, W)
      W <-  t(C) %*% Z[,idx]
      G_g <- W %*% G[idx,idx,drop=FALSE] %*% t(W)

      # Check estimability
      # P <- t(X_tilde) %*% ginv_sym_sparse(tcrossprod(X_tilde)) %*% X_tilde
      # c <- c(rep(0, ), W[,1] - W[,3])
      # plot(c - P %*% c)

      H2_Standard <-  h2_Standard_parameters(G_g, V, C, active_g)

    } else {
      g <- mf[[target]]
      gnames <- levels(g)
      Z_g <- Matrix::sparse.model.matrix(~ 0 + g)
      C <- Z_g %*% Matrix::Diagonal(x = 1 / as.numeric(Matrix::colSums(Z_g)))
      W <- t(C) %*% Z[,idx]
      G_g <- W %*% G[idx,idx,drop=FALSE] %*% t(W)

      H2_Standard <-  h2_Standard_parameters(G_g, V, C, NULL)
    }

  }

  return(H2_Standard)
}
