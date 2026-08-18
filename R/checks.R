#' Null coalescing operator
#'
#' Returns the left-hand side if it is not NULL, otherwise returns the right-hand side.
#'
#' @param x Left-hand side value
#' @param y Right-hand side value (default if x is NULL)
#'
#' @return x if x is not NULL, otherwise y
#' @noRd
#' @keywords internal
#' @examples
#' NULL heritable:::`%||%` "default" # returns "default"
#' "value" heritable:::`%||%` "default" # returns "value"
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' @keywords internal
initial_checks <- function(model, target, options) {
  if (options$check %||% TRUE) {
    # If there is more than one target, error
    check_target_single(target)

    # If the model is not Gaussian with an identity link, error
    check_model_family(model)

    # If model has not converged, warn
    check_model_convergence(model)

    # Check if target appears in the model
    check_target_exists(model, target)

    # Check if target is random
    if(!check_target_random(model, target)){
      cli::cli_abort(
        "The target {.var {target}} is fitted as both fixed and random effect"
      )
    }

    # Check if the target is specified as both fixed and random
    check_target_both(model, target)

  }
}

########################### Check all ##############################

# Check if model converged
#' @keywords internal
check_model_convergence.asreml <- function(model) {
  if (!model$converge) {
    cli::cli_warn(
      "The input model has not converged, estimates may be unreliable"
    )
  }
}

#' @keywords internal
check_model_convergence.lmerMod <- function(model) {
  if (model@optinfo$conv$opt != 0) {
    warning("The input model has not converged")
  }
}

#' @keywords internal
check_model_convergence.glmmTMB <- function(model) {
  conv <- model$fit$convergence
  pdHess <- isTRUE(model$sdr$pdHess)
  if ((!is.null(conv) && conv != 0) || !pdHess) {
    warning("The input model has not converged")
  }
}

#' @keywords internal
check_model_convergence <- function(model) {
  UseMethod("check_model_convergence")
}
.S3method("check_model_convergence", "asreml", check_model_convergence.asreml)
.S3method("check_model_convergence", "lmerMod", check_model_convergence.lmerMod)
.S3method("check_model_convergence", "glmmTMB", check_model_convergence.glmmTMB)


# Check the response distribution and link function
# Only Gaussian models with an identity link are currently supported for the
# `lme4` backend. Catches `glmer()` / `glmer.nb()` (class `glmerMod`) models,
# which otherwise fail later with an uninformative dispatch error.
#' @keywords internal
check_model_family.merMod <- function(model) {
  fam <- stats::family(model)
  if (!identical(fam$family, "gaussian") || !identical(fam$link, "identity")) {
    cli::cli_abort(
      c(
        "Only Gaussian models with an identity link are currently supported for {.pkg lme4} models.",
        "x" = "The supplied model has family {.val {fam$family}} with link {.val {fam$link}}.",
        "i" = "Heritability for generalized linear mixed models is not currently implemented."
      )
    )
  }
  invisible(model)
}

#' @keywords internal
check_model_family.default <- function(model) {
  invisible(model)
}

# Gate the `glmmTMB` backend to the model structures it can currently handle.
# The variance-component machinery assumes a plain Gaussian (identity-link)
# model with a homoscedastic residual `R = sigma^2 * I` and simple correlated
# random-effects blocks. Model features that break those assumptions are
# rejected here with an informative error instead of failing later with a
# cryptic message (or, worse, returning a silently wrong result).
#' @keywords internal
check_model_family.glmmTMB <- function(model) {
  fam <- model$modelInfo$family
  if (!identical(fam$family, "gaussian") || !identical(fam$link, "identity")) {
    cli::cli_abort(
      c(
        "Only Gaussian models with an identity link are currently supported for {.pkg glmmTMB} models.",
        "x" = "The supplied model has family {.val {fam$family}} with link {.val {fam$link}}.",
        "i" = "Heritability for generalized linear mixed models is not currently implemented."
      )
    )
  }

  # Zero-inflation: var_comp() only extracts the conditional model, so a
  # zero-inflation component would be silently dropped.
  zi <- stats::terms(model$modelInfo$allForm$ziformula)
  if (length(labels(zi)) > 0 || attr(zi, "intercept") == 1) {
    cli::cli_abort(
      c(
        "Zero-inflated {.pkg glmmTMB} models are not currently supported.",
        "i" = "Heritability is derived from the conditional model only; a zero-inflation component would be ignored."
      )
    )
  }

  # Non-constant dispersion: the backend assumes a homoscedastic residual
  # variance `R = sigma^2 * I`. A modelled `dispformula` also makes
  # `sigma()` return `NA`.
  dsp <- stats::terms(model$modelInfo$allForm$dispformula)
  if (length(labels(dsp)) > 0 || !is.finite(stats::sigma(model))) {
    cli::cli_abort(
      c(
        "{.pkg glmmTMB} models with a non-constant {.code dispformula} are not currently supported.",
        "i" = "The residual variance is assumed homoscedastic ({.code R = sigma^2 * I})."
      )
    )
  }

  # A weighted Gaussian fit has residual covariance `sigma^2 * diag(1 / w)`, but
  # the shared variance-component core rebuilds it as `sigma^2 * I`, dropping the
  # per-observation weighting, so the resulting heritability does not account for
  # the weights. How prior weights should enter a heritability computation is
  # itself unresolved and affects every backend alike, so a consistent fix
  # belongs in a dedicated cross-backend change rather than here. For now we warn
  # (rather than silently returning a weight-unaware result, or hard-blocking
  # only this backend while lme4 accepts weights without comment).
  # TODO(anu-aagi/heritable#52): handle prior weights consistently across backends.
  w <- model$frame[["(weights)"]]
  if (!is.null(w) && !isTRUE(all(w == 1))) {
    cli::cli_warn(
      c(
        "Prior {.code weights} are not accounted for in {.pkg glmmTMB} heritability estimates.",
        "i" = "The residual variance is treated as homoscedastic ({.code R = sigma^2 * I}), not the weighted fit's {.code sigma^2 * diag(1 / w)}, so the estimate ignores the weights and may be biased.",
        "i" = "Consistent handling of prior weights is tracked in anu-aagi/heritable#52."
      )
    )
  }

  # This first glmmTMB backend supports only plain random-effect *intercepts*:
  # single-column blocks whose one column is `(Intercept)`, i.e. `(1 | group)`
  # terms (including interactions such as `(1 | gen:block)`). That is the
  # canonical heritability setup and the model class lme4 and glmmTMB reproduce
  # identically. Everything else glmmTMB can fit is deferred to a follow-up,
  # because it packs more into a covariance block than the downstream term
  # mapping (one formula bar <-> one RE block) can currently align:
  #   * random slopes -- `(1 + x | gen)` (a multi-column "us" block) and
  #     `(1 + x || gen)` (a single `diag` block that findbars() splits into
  #     several bars, so bars and blocks no longer line up);
  #   * structured G-side covariance -- `cs()`, `ar1()`, `toep()`, `us(...)`, ...
  # TODO(anu-aagi/heritable #19 follow-up): support random slopes and structured
  # G-/R-side covariance (in particular on non-target effects), then relax this.
  #
  # Gating on "single column AND that column is the intercept" (rather than on
  # `blockCode == "us" && single column`) is deliberate and load-bearing:
  #   * a lone random *slope* `(0 + x | g)` is ALSO a single-column "us" block,
  #     and the uncorrelated `(1 + x || g)` can reach us pre-expanded as two
  #     single-column blocks on the same factor, `(1 | g) + (0 + x | g)`. Both
  #     must be rejected here, otherwise they slip through and either crash or
  #     mis-map downstream (the very failure this guard exists to prevent).
  #   * for a single intercept column the covariance code is necessarily plain
  #     ("us"/"diag"/"homdiag" all coincide), so `diag(1 | g)` -- identical to
  #     `(1 | g)` -- is correctly accepted without a special case.
  vcc  <- glmmTMB::VarCorr(model)$cond
  cnms <- model$modelInfo$reTrms$cond$cnms
  code <- vapply(vcc, function(v) {
    bc <- attr(v, "blockCode")
    if (is.null(bc)) "us" else names(bc)
  }, character(1))
  is_plain_intercept <- vapply(cnms, function(cc) {
    length(cc) == 1L && identical(unname(cc), "(Intercept)")
  }, logical(1))
  if (any(!is_plain_intercept)) {
    found <- unique(ifelse(code != "us", code, "random slope")[!is_plain_intercept])
    cli::cli_abort(
      c(
        "Only plain random-effect intercepts {.code (1 | group)} are currently supported for {.pkg glmmTMB} models.",
        "x" = "Found an unsupported random-effect structure: {.val {found}}.",
        "i" = "Random slopes and structured covariance (e.g. {.code cs()}, {.code ar1()}, {.code toep()}, {.code ||}) are planned for a future release."
      )
    )
  }

  invisible(model)
}

#' @keywords internal
check_model_family <- function(model) {
  UseMethod("check_model_family")
}
.S3method("check_model_family", "merMod", check_model_family.merMod)
.S3method("check_model_family", "default", check_model_family.default)
.S3method("check_model_family", "glmmTMB", check_model_family.glmmTMB)


#' Check whether the fitted model contains random terms not grouped by `target`
#' @keywords internal
check_all_random_terms_match_target <- function(model, target) {
  matched_grp <- pull_terms_without_specials(model)$random == target
  if(!all(matched_grp)){
    FALSE
  } else {
    TRUE
  }
}


# Target level checks
# Check if only one target has been supplied
#' @keywords internal
check_target_single <- function(target) {
  if (is.null(target)) {
    cli::cli_abort(
      "The target is {.value NULL}. Please specify the target name."
    )
  }
  if (length(target) > 1) {
    cli::cli_abort("Only one target can be supplied to calculate heritability")
  }
}

# Check if target is in model
#' @keywords internal
check_target_exists <- function(model, target) {
  model_terms <- pull_terms_without_specials(model)
  if (!target %in% c(model_terms$fixed, model_terms$random)) {
    cli::cli_abort(
      "The specified target {.code {target}} is not found in the model"
    )
  }
}

# Check if target is in fixed or random
#' @keywords internal
check_target_random <- function(model, target) {
  model_terms <- pull_terms_without_specials(model)
  if (target %in% model_terms$random) {
    TRUE
  } else {
    FALSE
  }
}


# Check if target is in both fixed and random
#' @keywords internal
check_target_both <- function(model, target) {
  model_terms <- pull_terms_without_specials(model)
  # Use regex to check for presence in both fixed and random
  if (
    any(grepl(target, model_terms$fixed, fixed = TRUE)) &&
      any(grepl(target, model_terms$random, fixed = TRUE))
  ) {
    cli::cli_abort("The target {.code {target}} is fitted as both fixed and random effect")
  }
}

# Check if the design matrix exists, otherwise builds one.
#' @keywords internal
check_deisgn_exsits <- function(model, build_design = TRUE, build_mf = TRUE, source = list()){

  if(!inherits(model, "asreml")){
    return(model)
  }

  if(build_design){
    # Get the design matrix
    design <- model$design
    if (is.null(design)) {
      cli::cli_inform("A design matrix was not found in the asreml object. Building a design matrix.")
      build_design <- TRUE
    } else {
      build_design <- FALSE
    }
  }

  if(build_mf){
    # Get the design matrix
    mf <- model$mf
    if (is.null(mf)) {
      cli::cli_inform("A model frame was not found in the asreml object. Building a model frame.")
      build_mf <- TRUE
    } else {
      build_mf <- FALSE
    }
  }

  if(build_design && build_mf){
    design_default <- asreml::asreml.options()$design
    asreml::asreml.options(design = TRUE)
    model <- asreml::update.asreml(model, model.frame = TRUE)
    asreml::asreml.options(design = design_default)
  }

  if(!build_design && build_mf){
    model <- asreml::update.asreml(model, model.frame = TRUE)
  }

  if(build_design && !build_mf){
    design_default <- asreml::asreml.options()$design
    asreml::asreml.options(design = TRUE)
    model <- asreml::update.asreml(model, model.frame = TRUE)
    asreml::asreml.options(design = design_default)
  }

  model
}


########################### Method specific check ##############################
# Helper function to check if GRM exists in environment
#' @keywords internal
check_GRM_exists <- function(model, target = NULL, source = list(), return = FALSE) {

  if(!is.list(source)){
    cli::cli_abort("source must be a named list.")
  }

  trms <- pull_terms(model)$random
  trms_no_special <- pull_terms_without_specials(model)$random

  if(!is.null(target)){
    contain_target <- sapply(trms_no_special, function(trm){
      target %in% stringr::str_split(trm, ":")[[1]]
    }, USE.NAMES = FALSE)
  } else {
    contain_target <- rep(TRUE, length(model$G.param))
  }

  trms <- lapply(model$G.param[contain_target], function(x) {
    x <- lapply(x[-1], function(y) {
      if(y[["model"]] == "vm") y[["facnam"]] else NULL
    }) |> unname()
    do.call(c, x)
  }) |> unname()
  trms <- do.call(c, trms)

  if (length(trms) > 0) {
    name_GRM <- stringr::str_match(
      trms,
      "source\\s*=\\s*([^\\s,\\)]+)"
    )[,2]
    na_idx <- is.na(name_GRM)

    name_GRM[na_idx] <- stringr::str_match(
      trms[na_idx],
      "vm\\([^,]+,\\s*([^\\s,\\)]+)"
    )[,2]

    name_GRM <- unique(name_GRM)

    for(x in name_GRM){

      if(!is.null(source[[x]])){

        if(!return) source[[x]] <- TRUE

      } else if (exists(x, envir = .GlobalEnv, inherits = FALSE)) {

        source[[x]] <- if(return) get(x, envir = .GlobalEnv, inherits = FALSE) else TRUE

      } else {
        # Source doesn't exist and not supplied
        cli::cli_abort("Cannot find the source {.code {x}}.")
      }
    }

  }
  source
}

#' Check target term specification for borad-sense heritability
#' For lme4, target as a random effect can only be specified once as (1|target)
#' For asreml, target as a random effect can only be specified once as target
#' @keywords internal
#' @noRd
check_model_specification.asreml <- function(model, target,
                                             type = c("broad_sense", "narrow_sense"),
                                             source = list(),
                                             ...){
  type <- match.arg(type)
  ran_trms <- pull_terms_without_specials(model)$random
  ran_trms_with_special <- pull_terms(model)$random
  spec <- attr(ran_trms_with_special, "spec")

  # Check GRM exists
  check_GRM_exists(model, target, source)

  if(target %in% ran_trms){
    if(type == "broad_sense"){
      if(sum(ran_trms == target)!=1){
        cli::cli_warn("The target {.code {target}} as a grouping variable should be specified once. Heritability calculation can be misleading.")
      } else {
        simple_model <- !any(ran_trms == target & spec != "id")
        if(!simple_model){
          cli::cli_warn("The target {.code {target}} should be modelled as a random term without special: {.code id({target})}. Heritability calculation can be misleading.")
        }
      }
    }
  }
}

#' @keywords internal
#' @noRd
check_model_specification.lmerMod <- function(model, target,
                                              type = c("broad_sense", "narrow_sense"),
                                              ...){
  type <- match.arg(type)
  ran_trms <- pull_terms_without_specials(model)$random
  ran_trms_with_special <- pull_terms(model)$random

  if(target %in% ran_trms){
    if(type == "broad_sense"){
      if(sum(ran_trms == target)!=1){
        cli::cli_warn("Duplicated random intercept detected for the target {.code {target}}. Heritability calculation can be misleading.")
      } else {
        grp_name <- attr(ran_trms_with_special, "grouping_variable")
        simple_intercept  <- attr(ran_trms_with_special, "simple_intercept")
        contain_intercept  <- attr(ran_trms_with_special, "contain_intercept")
        simple_model <- all(simple_intercept[contain_intercept & grp_name== target])
        if(!simple_model){
          cli::cli_warn("The target {.code {target}} should be modelled as a random term without special: {.code (1 | {target})}. Heritability calculation can be misleading.")
        }
      }
    }
  }

}

#' @keywords internal
#' @noRd
check_model_specification.glmmTMB <- function(model, target,
                                              type = c("broad_sense", "narrow_sense"),
                                              ...){
  # Body is purely attribute-based (pull_terms), so the lme4 logic is shared.
  check_model_specification.lmerMod(model, target, type = match.arg(type), ...)
}

#' @keywords internal
#' @noRd
check_model_specification <- function(model, target,
                                      type = c("broad_sense", "narrow_sense"),
                                      ...) {
  UseMethod("check_model_specification")
}
.S3method("check_model_specification", "asreml", check_model_specification.asreml)
.S3method("check_model_specification", "lmerMod", check_model_specification.lmerMod)
.S3method("check_model_specification", "glmmTMB", check_model_specification.glmmTMB)
