#' Pull fixed and random terms from a model formula
#'
#' Extract the labels of fixed and random terms from a model object that exposes
#' a formula with `fixed` and `random` components (for example objects produced
#' by asreml::asreml). The function returns a named list containing two character
#' vectors: `fixed` and `random`.
#'
#' @param model A fitted model object with a `formula` method that returns a
#'   list containing `fixed` and `random` formula components.
#' @return A named list with components:
#'   \item{fixed}{Character vector of labels for fixed-effect terms.}
#'   \item{random}{Character vector of labels for random-effect terms.}
#' @importFrom stats terms formula
#' @keywords internal

pull_terms.asreml <- function(model) {
    fixed_trms <- terms(formula(model)$fixed) |> labels()
    ran_trms <- terms(formula(model)$random) |> labels()
    return(list(fixed = fixed_trms, random = ran_trms))
}

#' @keywords internal
pull_terms.lmerMod <- function(model) {
    model_formula <- formula(model)
    term_labels <- attr(terms(model_formula), "term.labels")

    ran_trms <- names(lme4::ranef(model))
    fixed_trms <- setdiff(term_labels, paste0("1 | ", ran_trms))

    return(list(fixed = fixed_trms, random = ran_trms))
}

#' @keywords internal
pull_terms <- function(model) {
    UseMethod("pull_terms")
}
.S3method("pull_terms", "asreml", pull_terms.asreml)
.S3method("pull_terms", "lmerMod", pull_terms.lmerMod)

#' Fit the counterpart of an asreml model by swapping a specified term between fixed and random
#'
#' Fit a "counterpart" model to an existing asreml model by moving a specified
#' term from the random effects to the fixed effects or vice
#' versa. This is useful for procedures that require the same term to be fitted
#' as the opposite effect (for example, calculating Piepho's heritability).
#'
#'
#' @param model An existing fitted asreml model object. The function expects
#'   that formulas can be retrieved via formula(model)$fixed and
#'   formula(model)$random and that the model frame is available as model$mf.
#' @param target Character(1). Name of the term (e.g. a factor variable) to be
#'   switched between fixed and random effects. Must match one of the terms
#'   present in either the model's fixed or random formulas.
#'
#' @return A fitted asreml model object identical to the input model except
#'   that `target` has been moved from the random effects to the fixed effects
#'   or from the fixed effects to the random effects, depending on where it
#'   appeared in the original model. The returned object is produced by
#'   asreml::asreml and can be used for subsequent model comparisons or
#'   variance-component calculations.
#' @importFrom stats as.formula update
#' @keywords internal

fit_counterpart_model.asreml <- function(model, target = NULL) {
    # get the terms from model object
    fixed_trms <- pull_terms.asreml(model)$fixed
    ran_trms <- pull_terms.asreml(model)$random

    # when target is in random
    if (target %in% ran_trms) {
        #cli::cli_inform("{.var {target}} was fitted as a random effect. We will fit {.var {target}} as a fixed effect to calculate heritability.")
        # fit model with target as fixed effect
        model_counter <- asreml::update.asreml(model,
                                               fixed = as.formula(paste(". ~ . +", target)),
                                               random =  as.formula(paste("~ . -", target)))
    } else if (target %in% fixed_trms) { # when target is in fixed
        #cli::cli_inform("{.var {target}} was fitted as a fixed effect. We will fit {.var {target}} as a random effect to calculate heritability.")
        # fit model with target as random effect
        model_counter <- asreml::update.asreml(model,
                                               fixed = as.formula(paste(". ~ . -", target)),
                                               random =  as.formula(paste("~ . +", target)))
    } else {
        cli::cli_abort("{.var {target}} not found in either fixed or random effects of the model.")
    }

    return(model_counter)
}

#' @keywords internal
fit_counterpart_model.lmerMod <- function(model, target = NULL) {
    # get the terms from model object
    fixed_trms <- pull_terms.lmerMod(model)$fixed
    ran_trms <- pull_terms.lmerMod(model)$random

    current_formula <- formula(model)
    random_sym <- reformulas::findbars(current_formula)

    # Non-target random effects
    other_re <- random_sym[[which(sapply(random_sym, function(x) !any(grepl(target, deparse(x)))))]]

    # If target is in random effects
    if (target %in% ran_trms) {
        #cli::cli_inform("{.var {target}} was fitted as a random effect. We will fit {.var {target}} as a fixed effect to calculate heritability.")

        updated_formula <-
            reformulas::nobars_(current_formula) |> # Remove random effect terms
            update(paste(". ~ . +", target)) # Add target as a fixed effect
    } else if (target %in% fixed_trms) { # If target is in fixed effects
        #cli::cli_inform("{.var {target}} was fitted as a fixed effect. We will fit {.var {target}} as a random effect to calculate heritability.")
        # Create new formula with target as random effect
        updated_formula <-
            reformulas::nobars_(current_formula) |> # Remove random effect terms
            update(as.formula(paste(". ~ . -", target))) |> # Remove target from fixed effects
            update(as.formula(paste(". ~ . + (1 |", target, ")"))) # Add target as random effect
    } else {
        cli::cli_abort("{.var {target}} not found in either fixed or random effects of the model.")
    }

    if (!is.null(other_re)) { # If there are other random effects, add them back in
        updated_formula <-
            update(updated_formula, as.formula(paste(". ~ . + (", deparse(other_re), ")")))
    }

    # Refit the model
    refit_model <- update(model, formula = updated_formula)
    return(refit_model)
}


#' @keywords internal
fit_counterpart_model <- function(model, target = NULL) {
    UseMethod("fit_counterpart_model")
}
.S3method("fit_counterpart_model", "asreml", fit_counterpart_model.asreml)
.S3method("fit_counterpart_model", "lmerMod", fit_counterpart_model.lmerMod)

#' Print method for heritable objects
#'
#' @param x An object of class "heritable"
#' @param digits Number of digits to print
#' @param ... Additional arguments passed to print
#'
#' @export
print.heritable <- function(x, digits = getOption("digits"), ...) {
    print(unclass(x))
    # # Format all values to specified digits
    # x_rounded <- round(x, digits)
    #
    # # If multiple methods, print with method names
    # if (length(x) > 1) {
    #     cli::cli_h2("Heritability estimates:")
    #     for (method in names(x_rounded)) {
    #         cli::cli_text("{.var {method}}: {x_rounded[[method]]}")
    #     }
    # } else {
    #     # For single method, print simple value
    #     cli::cli_text("Heritability ({.var {names(x_rounded)}}): {x_rounded}")
    # }
    #
    # invisible(x)
}
