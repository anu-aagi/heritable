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
    fixed_trms <- attr(model$formulae$fixed, "term.labels")
    ran_trms <- attr(model$formulae$random, "term.labels")
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


#' @keywords internal
pull_terms_without_specials <- function(model) {
  UseMethod("pull_terms_without_specials")
}

#' @keywords internal
pull_terms_without_specials.lmerMod <- function(model) {
  model_terms <- pull_terms(model)
  model_terms
}

semivariance <- function(X) {
  n <- nrow(X)
  1/(n - 1) * (sum(diag(X)) - 1 / n * sum(X))
}

#' @keywords internal
pull_terms_without_specials.asreml <- function(model) {
  model_terms <- pull_terms(model)
  pattern <- paste0("^(",
                    paste0(asreml_Spcls, collapse = "|"),
                    ")\\(([^,]+),?.*\\)")
  clean_which <- stringr::str_detect(model_terms$fixed, pattern)
  model_terms$fixed[clean_which] <- stringr::str_extract(model_terms$fixed[clean_which],
                                                         pattern,
                                                         group = 2)
  clean_which <- stringr::str_detect(model_terms$random, pattern)
  model_terms$random[clean_which] <- stringr::str_extract(model_terms$random[clean_which],
                                                          pattern,
                                                          group = 2)
  model_terms
}


asreml_Spcls <- c("con", "C", "lin", "pow", "pol", "leg", "spl", "dev", "ped",
                  "ide", "giv", "vm", "ma", "at", "dsum", "and", "grp", "mbf",
                  "sbs", "gpf", "uni", "id", "idv", "idh", "ar1", "ar1v", "ar1h",
                  "ar2", "ar2v", "ar2h", "ar3", "ar3v", "ar3h", "sar", "sarv",
                  "sarh", "sar2", "sar2v", "sar2h", "ma1", "ma1v", "ma1h", "ma2",
                  "ma2v", "ma2h", "arma", "armav", "armah", "cor", "corv", "corh",
                  "corb", "corbv", "corbh", "corg", "corgv", "corgh", "diag", "us",
                  "sfa", "chol", "cholc", "ante", "exp", "expv", "exph", "iexp",
                  "iexpv", "iexph", "aexp", "aexpv", "bexpv", "aexph", "gau", "gauv",
                  "gauh", "lvr", "lvrv", "lvrh", "igau", "igauv", "igauh", "agau",
                  "agauv", "agauh", "ieuc", "ieucv", "ieuch", "ilv", "ilvv", "ilvh",
                  "sph", "sphv", "sphh", "cir", "cirv", "cirh", "mtrn", "mtrnv",
                  "mtrnh", "mthr", "facv", "fa", "rr", "str", "own")



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
    trms <- pull_terms.lmerMod(model)

    # If target is in random effects
    if (target %in% trms$random) {
      refit_model <- update(model, as.formula(paste(". ~ . - (1|", target, ") + ", target)))
    } else if (target %in% trms$fixed) { # If target is in fixed effects
      refit_model <- update(model, as.formula(paste(". ~ . + (1|", target, ") - ", target)))
    } else {
        cli::cli_abort("{.var {target}} not found in either fixed or random effects of the model.")
    }
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
