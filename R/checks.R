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

    # Check if G x E is fitted
    check_G_by_E_exists(model, target)

    if (check_target_both(model, target)) {
      cli::cli_abort("The target {.var {target}} is fitted as both fixed and random effect")
    }

    if (options$target_once %||% FALSE) {
      # Check if target only appears exactly once in the model
      check_target_appears_once(model, target)
    }

    # Check if target is random or fixed
    if (!check_target_random(model, target)) {
      cli::cli_abort("Heritability can only be calculated if {.value target} is a random effect.")
    }
  }
}

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

check_target_appears_once <- function(model, target) {
  model_terms <- pull_terms_without_specials(model)
  form <- as.formula(paste("~", paste(
    c(
      model_terms$fixed,
      model_terms$random
    ),
    collapse = " + "
  )))
  form <- update(form, paste0("~ . - ", target))
  fcts <- rownames(attr(terms(form), "factors"))
  if (target %in% fcts) {
    cli::cli_abort(
      "The specified target {.code {target}} is found in multiple terms. Please specify only once."
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
    TRUE
  } else {
    FALSE
  }
}

# Check length of random effects
#' @keywords internal
check_single_random_effect <- function(terms) {
  if (length(terms$random) == 1) { #
    TRUE
  } else {
    FALSE
  }
}

# Helper function to check if GRM exists in environment
#' @keywords internal
check_GRM_in_environment <- function(model, target) {
  vpars <- names(model$vparameters)
  env <- attr(model$formulae$random, ".Environment")
  w <- grepl(paste0("^vm\\(", target), vpars)
  if (sum(w) == 1) {
    target_vm <- vpars[w]
    #name_GRM <- stringr::str_extract(vpars[w], paste0("vm\\(", target, ", (.+)\\)"), group = 1)
    name_GRM <- stringr::str_match(
      vpars[w],
      paste0("vm\\(", target, "\\s*,\\s*([^,\\)]+)")
    )[,2]
    if (exists(name_GRM, envir = env, inherits = FALSE)) {
      return(TRUE)
    }
  }
  FALSE
}

#' Check if GRM is supplied if not search environment
#' @keywords internal
check_GRM_exists <- function(model, target, source = NULL){
  # Is source supplied?
  if(!is.null(source)){
    TRUE
  } else if (check_GRM_in_environment(model, target)) {
    # Is it in the environment?
    TRUE
  } else {
    # Source doesn't exist and not supplied
    cli::cli_abort("Cannot find the source for {.code vm({target}, ...)}.")
  }
}

#' Check if G X E is fitted
#' @keywords internal
#' @noRd
check_G_by_E_exists.asreml <- function(model, target){
  nms <- names(model$vparameters)
  pattern <- paste0("(^|:)", target, "($|:)")

  if(sum(grepl(pattern, nms)) > 1){
    cli::cli_abort("G x E models are currently not supported")
  } else {
    FALSE
  }
}

#' @keywords internal
#' @noRd
check_G_by_E_exists.lmerMod <- function(model, target) {
  bars <- reformulas::findbars(formula(model))
  groups <- vapply(bars, function(b) paste(deparse(b[[3]]), collapse = ""), "")
  pattern <- paste0("(^|:)", target, "($|:)")

    if(sum(grepl(pattern, groups)) > 1){
      cli::cli_abort("G x E models are currently not supported")
    }
}

#' @keywords internal
#' @noRd
check_G_by_E_exists <- function(model, target) {
  UseMethod("check_G_by_E_exists")
}
.S3method("check_G_by_E_exists", "asreml", check_G_by_E_exists.asreml)
.S3method("check_G_by_E_exists", "lmerMod", check_G_by_E_exists.lmerMod)



