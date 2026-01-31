#' @noRd
#' @keywords internal
h2_Cullis.asreml <- function(model, target = NULL, options = NULL) {

  initial_checks(model, target, options)

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  vm <- target_vm_term_asreml(model, target)

  n_g <-  model$noeff[[vm$target_vm]]
  vc_g <- model$vparameters[[vm$target_vm]] * model$sigma2 * semivariance(vm$GRM)


  vdBLUP_mat <- predict(model,
                        classify = target,
                        only = target,
                        sed = TRUE,
                        trace = FALSE,
  )$sed^2

  vd_BLUP_avg <- mean(vdBLUP_mat[upper.tri(vdBLUP_mat, diag = FALSE)])

  H2_Cullis_parameters(vd_BLUP_avg, vc_g)
}

#' @keywords internal
h2_Oakey.asreml <- function(model, target = NULL, source = NULL, options = NULL) {
  initial_checks(model, target, options)

  if(check_GRM_exists(model, target, source)){
    vm <- target_vm_term_asreml(model, target)
    n_g <- model$noeff[[vm$target_vm]]
    Gg_inv <- 1 / (model$vparameters[[vm$target_vm]] * model$sigma2) * vm$GRMinv
    vcov_g <- predict(model,
                      classify = target,
                      only = target,
                      vcov = TRUE,
                      trace = FALSE
    )$vcov

    names_clean <- stringr::str_remove(
      rownames(model$coefficients$random),
      "vm\\(.+\\)_"
    )

    dimnames(vcov_g) <- list(names_clean, names_clean)

    vcov_g <- vcov_g[rownames(Gg_inv),colnames(Gg_inv)]

    H2_Oakey_parameters(Gg_inv, vcov_g)
  }
}

#' @keywords internal
h2_Delta_pairwise.asreml <- function(model, target = NULL, source = NULL, type = NULL, options = NULL) {
  initial_checks(model, target, options)

  if(check_GRM_exists(model, target, source)){

  vm <- target_vm_term_asreml(model, target)
  n_g <- model$noeff[[vm$target_vm]]
  Gg <- model$vparameters[[vm$target_vm]] * model$sigma2 * solve(vm$GRMinv)

  if (type == "BLUP") {
    gpred <- predict(model, classify = target, sed = TRUE, trace = FALSE)
    Vd_g <- gpred$sed^2 # Variance of difference
    genotype_names <- gpred$pvals[[target]] # list of genotype names
    dimnames(Vd_g) <- list(genotype_names, genotype_names) # name the covariance matrix
    h2_Delta_BLUP_parameters(Gg, Vd_g)
  } else if (type == "BLUE") {
    model_fix <- fit_counterpart_model.asreml(model, target)
    gpred <- predict(model_fix, classify = target, sed = TRUE, trace = FALSE)
    Vd_g <- gpred$sed^2 # Variance of difference
    genotype_names <- gpred$pvals[[target]] # list of genotype names
    dimnames(Vd_g) <- list(genotype_names, genotype_names) # name the covariance matrix
    h2_Delta_BLUE_parameters(Gg, Vd_g)
  }
  }
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
    vc <- var_comp(model, target, calc_C22 = TRUE, marginal, stratification)
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
    vc <- var_comp(model, target, calc_C22 = TRUE, marginal, stratification)
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
    vc <- var_comp(model, target, calc_C22 = FALSE, marginal, stratification)
  }
  if(!vc$main){
    cli::cli_warn("Piepho heritability can only be computed on main genetic effects, retuning {.value {NA}}.")
    return(NA)
  }
  G_g <- vc$G_g
  s2_g <- mean(diag(G_g))

  conterpart <- fit_counterpart_model(model, target)

  delta <- predict(conterpart,
    classify = target,
    sed = TRUE,
    trace = FALSE
  )$sed^2

  delta_avg <- mean(delta[upper.tri(delta, diag = FALSE)])

  # Calculate Piepho's H2
  H2_Piepho <- H2_Piepho_parameters(s2_g, delta_avg)

  return(H2_Piepho)
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
    H2_Delta <- H2_Delta_BLUE_pairwise.asreml(model, target, options,
                                               marginal, stratification, vc)
  } else if (type == "BLUP") {
    H2_Delta <- H2_Delta_BLUP_pairwise.asreml(model, target, options,
                                               marginal, stratification, vc)
  }

  return(H2_Delta)
}

#' @keywords internal
H2_Delta_BLUP_pairwise.asreml<- function(model,
                                           target = NULL,
                                           options = NULL,
                                           marginal = TRUE,
                                           stratification = NULL,
                                           vc = NULL) {
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
    vc <- var_comp(model, target, calc_C22 = TRUE, marginal, stratification)
  }
  s2_g <- mean(diag(vc$G_g))
  C22_g <- vc$C22_g

  # Compute variance of difference from PEV
  delta <- var_diff(C22_g)
  diag(delta) <- NA
  dimnames(delta) <- list(vc$gnames, vc$gnames)

  # H2 Delta BLUP
  H2_Delta_BLUP <- H2_Delta_parameters(s2_g, delta, "BLUP")

  dimnames(H2_Delta_BLUP) <- dimnames(delta)

  H2_Delta_BLUP
}

#' @keywords internal
H2_Delta_BLUE_pairwise.asreml<- function(model,
                                         target = NULL,
                                         options = NULL,
                                         marginal = TRUE,
                                         stratification = NULL,
                                         vc = NULL) {
  initial_checks(model, target, options)

  if (options$check %||% TRUE) {
    # Check correct model specification.
    check_model_specification(model, target, "broad_sense")
  }

  # Check if target is random or fixed
  if (!check_target_random(model, target)) {
    return(NA)
  }

  # Extract vc_g and vc_e
  if(is.null(vc)){
    vc <- var_comp(model, target, calc_C22 = FALSE, marginal, stratification)
  }
  if(!vc$main){
    cli::cli_warn("Delta (BLUE) heritability can only be computed on main genetic effects, retuning {.value {NA}}.")
    return(NA)
  }
  s2_g <- mean(diag(vc$G_g))

  conterpart <- fit_counterpart_model(model, target)

  # Get delta
  delta <- predict(conterpart, classify = target, sed = TRUE, trace = FALSE)$sed^2
  diag(delta) <- NA
  dimnames(delta) <- list(vc$gnames, vc$gnames)

  # H2 Delta BLUE
  H2_Delta_BLUE <- H2_Delta_parameters(s2_g, delta, "BLUE")

  dimnames(H2_Delta_BLUE) <- dimnames(delta)

  H2_Delta_BLUE
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
    model <- check_deisgn_exsits(model)
    mf <- model$mf
    G_g <- var_comp(model, target, calc_C22 = TRUE, marginal, stratification)$G_g
  } else {
    model <- check_deisgn_exsits(model, build_design = FALSE)
    mf <- model$mf
    G_g <- vc$G_g
  }
  s2_g <- mean(diag(G_g))

  # Get residual variance
  s2_eps <- model$sigma2

  n_r <- table(model$mf[[target]])

  H2_Standard <- H2_Standard_parameters(s2_g, s2_eps, n_r)

  return(H2_Standard)
}
