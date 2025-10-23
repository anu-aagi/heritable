
#' Pull fixed and random terms from a model formula
#'
#' Extract the labels of fixed and random terms from a model object that exposes
#' a formula with `fixed` and `random` components (for example objects produced
#' by asreml::asreml). The function returns a named list containing two character
#' vectors: `fixed` and `random`.
#'
#' @param model A fitted model object with a `formula` method that returns a
#'   list containing `fixed` and `random` formula components. Typically this is
#'   an `nlme::lme` object or any model providing the same structure.
#' @return A named list with components:
#'   \item{fixed}{Character vector of labels for fixed-effect terms.}
#'   \item{random}{Character vector of labels for random-effect terms.}
#' @importFrom stats terms formula
#' @keywords internal

pull_terms <- function(model){
    fixed_trms <- terms(formula(model)$fixed) |> labels() 
    ran_trms <- terms(formula(model)$random) |> labels() 
    return(list(fixed = fixed_trms, random = ran_trms))
}   

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

fit_counterpart_model.asreml <- function(model, target = NULL){
    # get the terms from model object
    fixed_trms <- pull_terms(model)$fixed
    ran_trms <- pull_terms(model)$random
    
    # TODO: when target is in both random and fixed

    # when target is in random
    if(target %in% ran_trms){
    cli::cli_inform("{.var {target}} was fitted as a random effect. We will fit {.var {target}} as a fixed effect to calculate Piepho's heritability.")
    # fit model with target as fixed effect
    model_counter <- asreml::asreml(
        fixed = update(formula(model)$fixed, as.formula(paste(". ~ . +", target))),
        random = update(formula(model)$random, as.formula(paste("~ . -", target))), 
        data = model$mf,
        trace = FALSE
    )
    } else if (target %in% fixed_trms){ # when target is in fixed
    cli::cli_inform("{.var {target}} was fitted as a fixed effect. We will fit {.var {target}} as a random effect to calculate Piepho's heritability.")
    # fit model with target as random effect
    model_counter <- asreml::asreml(
        fixed = update(formula(model)$fixed, as.formula(paste(". ~ . -", target))),
        random = update(formula(model)$random, as.formula(paste("~ . +", target))),
        data = model$mf,
        trace = FALSE
    )
    } else {
    cli::cli_abort("{.var {target}} not found in either fixed or random effects of the model.")
    }
    
    return(model_counter)

}
