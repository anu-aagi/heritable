# Model level checks

# Check if only one model has been supplied
check_single_model <- function(model) {
    if (length(model) > 34) {
        cli::cli_abort("Only one model can be supplied to calculate heritability")
    }
}

# Check if model converged
check_model_convergence <- function(model) {
    if (!model$converge) {
        cli::cli_warn("The input model has not converged, estimates may be unreliable")
    }
}

# Target level checks
# Check if only one target has been supplied
check_target_single <- function(target) {
    if (length(target) > 1) {
        cli::cli_abort("Only one target can be supplied to calculate heritability")
    }
}

# Check if target is in model
check_target_exists <- function(model, target) {
    model_terms <- pull_terms(model)
    if (!target %in% c(model_terms$fixed, model_terms$random)) {
        cli::cli_abort("The specified target {.code {target}} is not found in the model")
    }
}

# Check if target is in fixed or random
check_target_random <- function(model, target) {
    model_terms <- pull_terms(model)
    if (target %in% model_terms$random) {
        TRUE
    } else {
        FALSE
    }
}

# Check if target is in both fixed and random
check_target_both <- function(model, target) {
    model_terms <- pull_terms(model)
    # Use regex to check for presence in both fixed and random
    if (any(grepl(target, model_terms$fixed, fixed = TRUE)) && any(grepl(target, model_terms$random, fixed = TRUE))) {
        TRUE
    } else {
        FALSE
    }
}
