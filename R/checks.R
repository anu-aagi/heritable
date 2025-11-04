# Model level checks

# Check if only one model has been supplied
check_single_model <- function(model) {
    if (length(model) > 1) {
        cli::cli_abort("Only one model can be supplied to calculate heritability")
    }
}

# Check if model converged
check_model_convergence <- function(model) {
    if (!model$converge) {
        cli::cli_warn("The input model has not converged, estimates may be unreliable")
    }
}

# Check if model is of a supported
## Check if model is of class asreml

## Check if model is of class lme4




# Target level checks
# Check if only one target has been supplied

# Check if target is in model

# Check if target is in fixed or random
