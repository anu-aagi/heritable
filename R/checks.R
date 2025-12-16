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
#' NULL %||% "default" # returns "default"
#' "value" %||% "default" # returns "value"
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

    if (check_target_both(model, target)) {
      cli::cli_abort("The target {.var {target}} is fitted as both fixed and random effect")
    }

    if (options$target_once %||% TRUE) {
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
