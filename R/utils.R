
pull_terms <- function(model){
    fixed_trms <- terms(formula(model)$fixed) |> labels() 
    ran_trms <- terms(formula(model)$random) |> labels() 
    return(list(fixed = fixed_trms, random = ran_trms))
}   

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
