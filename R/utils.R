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

#' @keywords internal
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

#' @keywords internal
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

#' @keywords internal
target_vm_term_asreml <- function(model, target) {
  vpars <- names(model$vparameters)
  env <- attr(model$formulae$random, ".Environment")
  w <- grepl(paste0("^vm\\(", target), vpars)
  if(sum(w) == 1) {
    target_vm <- vpars[w]
    name_GRM <- stringr::str_extract(vpars[w], paste0("vm\\(", target, ", (.+)\\)"), group = 1)
    if(exists(name_GRM, envir = env)) {
      GRM_source <- get(name_GRM, envir = env)
      if(is.data.frame(GRM_source) & ncol(GRM_source) == 3) {
        GRMinv <- solve(sp2Matrix(GRM_source))
      } else {
        GRMinv <- solve(GRM_source)
      }
      if(inherits(GRM_source, "ginv") || isTRUE(attr(GRM_source, "INVERSE"))) {
        GRMinv <- GRM_source
      }
    } else {
      cli::cli_abort("Cannot get the source {.value target_vm} for vm().")
    }
    return(list(target_vm = vpars[w],
                GRMinv = GRMinv))
  } else {
    cli::cli_abort("The {.value target} should be wrapped with vm() in the model with a known relationship matrix.")
  }
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

fit_counterpart_model.asreml <- function(model, target = NULL) {
    # get the terms from model object
    fixed_trms <- pull_terms.asreml(model)$fixed
    ran_trms_without_specials <- pull_terms_without_specials.asreml(model)$random
    ran_trms <- pull_terms.asreml(model)$random

    # when target is in random
    if (target %in% ran_trms_without_specials) {
      if(target %in% ran_trms) {
        model_counter <- update(model,
                                               fixed = as.formula(paste(". ~ . +", target)),
                                               random =  as.formula(paste("~ . -", target)),
                                               trace = FALSE)
      } else {
        target_spcl <- ran_trms[which(ran_trms_without_specials == target)]
        model_counter <- update(model,
                                               fixed = as.formula(paste(". ~ . +", target)),
                                               random =  as.formula(paste("~ . -", target_spcl)),
                                               trace = FALSE)
      }
    } else if (target %in% fixed_trms) { # when target is in fixed
        model_counter <- update(model,
                                               fixed = as.formula(paste(". ~ . -", target)),
                                               random =  as.formula(paste("~ . +", target)),
                                               trace = FALSE)
    } else {
        cli::cli_abort("{.var {target}} not found in either fixed or random effects of the model.")
    }
    check_model_convergence(model_counter)

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
    check_model_convergence(refit_model)
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
#' @noRd
#' @export
print.heritable <- function(x, digits = getOption("digits"), ...) {
    attr(x, "model") <- NULL
    attr(x, "target") <- NULL
    print(unclass(x))
}

#' @keywords internal
build_precompiled_vignette <- function(){
  knitr::knit(here::here("vignettes/heritable.Rmd.orig"),
              output = here::here("vignettes/heritable.Rmd"))
}

#' @keywords internal
#' @importFrom methods canCoerce hasMethod as
sp2Matrix <- function (x, dense = FALSE, triplet = FALSE)
{
  triplet <- ifelse(triplet, yes = "T", no = "C")
  A <- Matrix::sparseMatrix(x[, 1], x[, 2], x = x[, 3], repr = triplet,
                    symmetric = TRUE)
  if (dense) {
    if (canCoerce(A, "packedMatrix")) {
      A <- as(A, "packedMatrix")
    }
    else if (canCoerce(A, "dsyMatrix") && hasMethod("coerce",
                                                    c("dsyMatrix", "dspMatrix"))) {
      A <- as(as(A, "dsyMatrix"), "dspMatrix")
    }
    else {
      stop("Unable to return a dense matrix")
    }
  }
  if (inherits(x, "ginv")) {
    dimnames(A) <- list(attr(x, "rowNames"), attr(x, "rowNames"))
    attr(A, "INVERSE") <- TRUE
    for (i in c("inbreeding", "logdet", "geneticGroups")) attr(A,
                                                               i) <- attr(x, i)
  }
  else {
    att <- c("rowNames", "INVERSE")
    w <- which(is.element(att, names(attributes(x))))
    if (length(w) > 0) {
      for (i in w) {
        if (i == 1)
          dimnames(A) <- list(attr(x, "rowNames"), attr(x,
                                                        "rowNames"))
        else attr(A, att[i]) <- attr(x, att[i])
      }
    }
  }
  return(A)
}
