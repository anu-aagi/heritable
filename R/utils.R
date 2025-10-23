
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
#'
#' @details
#' - If `target` is found among the model's random terms, the function will:
#'     1. Inform the user that the term was random and will be fitted as fixed.
#'     2. Refit the model with `target` added to the fixed formula and removed
#'        from the random formula.
#' - If `target` is found among the model's fixed terms, the function will:
#'     1. Inform the user that the term was fixed and will be fitted as random.
#'     2. Refit the model with `target` removed from the fixed formula and
#'        added to the random formula.
#' - If `target` is not found in either set of terms, the function aborts with
#'   an error via cli::cli_abort.
#'
#' @note
#' - The function relies on asreml formula extraction and on the data stored in
#'   the original model object; ensure those components are present before
#'   calling this helper.
#' - Re-fitting models may take time and may change convergence/status
#'   information; inspect the returned model for warnings or fit issues.
#'
#' @examples
#' \dontrun{
#' # assume `m` is a previously fitted asreml model containing a term "geno"
#' # fitted as random; refit with "geno" as fixed to compute fixed-effect
#' # summaries or heritability estimates that require a fixed genotype.
#' m_counter <- fit_counterpart_model.asreml(m, target = "geno")
#' }
#'
#' @seealso asreml::asreml, cli::cli_inform, cli::cli_abort
#' @export
#' @importFrom asreml asreml
#' @importFrom cli cli_inform cli_abort
fit_counterpart_model.asreml <- function(model, target = NULL){
    # get the terms from model object
    fixed_trms <- pull_terms(model)$fixed
    ran_trms <- pull_terms(model)$random
    
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
