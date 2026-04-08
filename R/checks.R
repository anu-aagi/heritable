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
check_model_convergence <- function(model) {
  UseMethod("check_model_convergence")
}
.S3method("check_model_convergence", "asreml", check_model_convergence.asreml)
.S3method("check_model_convergence", "lmerMod", check_model_convergence.lmerMod)


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

  if(build_design ||  build_mf){
    source <- check_GRM_exists(model = model, source = source, return = TRUE)
    if(length(source) > 0){
      list2env(source, environment())
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
    cli::cli_abort("souce must be a named list.")
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
check_model_specification <- function(model, target,
                                      type = c("broad_sense", "narrow_sense"),
                                      ...) {
  UseMethod("check_model_specification")
}
.S3method("check_model_specification", "asreml", check_model_specification.asreml)
.S3method("check_model_specification", "lmerMod", check_model_specification.lmerMod)

#' Check if asreml is installed and load it
#'
#' Checks whether asreml is installed, errors if not, and loads it if present.
#'
#' @keywords internal
#' @noRd

check_and_load_asreml <- function() {
  if (!requireNamespace("asreml", quietly = TRUE)) {
    cli::cli_abort(
      "The {.pkg asreml} package is required for this function.
       Please install it before proceeding."
    )
  }
  invisible(library("asreml", character.only = TRUE))
}
