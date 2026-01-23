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
  ran_trms_formula <-
    stringr::str_extract_all(deparse1(model_formula), "(?<=\\()[^|()]+\\|+[^|()]+(?=\\))")[[1]]
  fixed_trms <- setdiff(term_labels, ran_trms_formula)
  ran_trms <- sapply(
    reformulas::findbars(model_formula), deparse
  )

  return(list(fixed = fixed_trms, random = ran_trms))
}

#' @keywords internal
pull_terms <- function(model) {
  UseMethod("pull_terms")
}
.S3method("pull_terms", "asreml", pull_terms.asreml)
.S3method("pull_terms", "lmerMod", pull_terms.lmerMod)


#' @keywords internal
pull_terms_without_specials.lmerMod <- function(model) {
  model_terms <- pull_terms(model)
  model_terms$random <- sapply(
    model_terms$random,
    function(frm){
      stringr::str_extract(frm, "(?<=\\|\\ ).+")
    }
  ) |> unname()
  model_terms
}

#' @keywords internal
semivariance <- function(X) {
  n <- nrow(X)
  1 / (n - 1) * (sum(diag(X)) - 1 / n * sum(X))
}

#' @keywords internal
pull_terms_without_specials.asreml <- function(model) {
  model_terms <- pull_terms(model)
  pattern <- paste0(
    "^(",
    paste0(asreml_Spcls, collapse = "|"),
    ")\\(([^,]+),?.*\\)"
  )
  clean_which <- stringr::str_which(model_terms$fixed, pattern)
  model_terms$fixed[clean_which] <- stringr::str_extract(model_terms$fixed[clean_which],
    pattern,
    group = 2
  )
  clean_which <- stringr::str_which(model_terms$random, pattern)
  model_terms$random[clean_which] <- stringr::str_extract(model_terms$random[clean_which],
    pattern,
    group = 2
  )
  model_terms
}

#' @keywords internal
asreml_Spcls <- c(
  "con", "C", "lin", "pow", "pol", "leg", "spl", "dev", "ped",
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
  "mtrnh", "mthr", "facv", "fa", "rr", "str", "own"
)

#' @keywords internal
pull_terms_without_specials <- function(model) {
  UseMethod("pull_terms_without_specials")
}
.S3method("pull_terms_without_specials", "asreml", pull_terms_without_specials.asreml)
.S3method("pull_terms_without_specials", "lmerMod", pull_terms_without_specials.lmerMod)


#' @keywords internal
target_vm_term_asreml <- function(model, target) {
  vpars <- names(model$vparameters)
  env <- attr(model$formulae$random, ".Environment")
  w <- grepl(paste0("^vm\\(", target), vpars)
  if (sum(w) == 1) {
    target_vm <- vpars[w]
    # name_GRM <- stringr::str_extract(vpars[w], paste0("vm\\(", target, ", (.+)\\)"), group = 1)
    name_GRM <- stringr::str_match(
      vpars[w],
      paste0("vm\\(", target, "\\s*,\\s*([^,\\)]+)")
    )[,2]
    if (exists(name_GRM, envir = env, inherits = FALSE)) {
      GRM_source <- get(name_GRM, envir = env)
      if (is.data.frame(GRM_source) & ncol(GRM_source) == 3) {
        #TODO: Change in future because GRM_source could be singular
        GRMinv <- solve(sp2Matrix(GRM_source))
      }
      if (inherits(GRM_source, "ginv") || isTRUE(attr(GRM_source, "INVERSE"))) {
        GRMinv <- GRM_source
      } else {
        GRMinv <- MASS::ginv(GRM_source)
      }
    } else {
      cli::cli_abort("Cannot get the source {.value target_vm} for vm().")
    }
    # Assign names
    dimnames(GRMinv) <- dimnames(GRM_source)

    return(list(
      target_vm = vpars[w],
      GRM_source = GRM_source,
      GRMinv = GRMinv
    ))
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
#' @param model An existing fitted asreml model object.
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
#' @noRd
fit_counterpart_model.asreml <- function(model, target = NULL) {
  # get the terms from model object
  fixed_trms <- pull_terms.asreml(model)$fixed
  ran_trms_without_specials <- pull_terms_without_specials.asreml(model)$random
  ran_trms <- pull_terms.asreml(model)$random

  # when target is in random
  if (target %in% ran_trms_without_specials) {
    if (target %in% ran_trms) {
      model_counter <- update(model,
        fixed = as.formula(paste(". ~ . +", target)),
        random =  as.formula(paste("~ . -", target)),
        trace = FALSE
      )
    } else {
      target_spcl <- ran_trms[which(ran_trms_without_specials == target)]
      model_counter <- update(model,
        fixed = as.formula(paste(". ~ . +", target)),
        random =  as.formula(paste("~ . -", target_spcl)),
        trace = FALSE
      )
    }
  } else if (target %in% fixed_trms) { # when target is in fixed
    model_counter <- update(model,
      fixed = as.formula(paste(". ~ . -", target)),
      random =  as.formula(paste("~ . +", target)),
      trace = FALSE
    )
  } else {
    cli::cli_abort("{.var {target}} not found in either fixed or random effects of the model.")
  }
  check_model_convergence(model_counter)

  return(model_counter)
}

#' @importFrom stats lm model.frame
#' @keywords internal
fit_counterpart_model.lmerMod <- function(model, target = NULL) {
  # get the terms from model object
  trms <- pull_terms_without_specials(model)
  trms <- lapply(trms, unique)
  fixed_trms <- trms$fixed
  ran_trms <- trms$random

  # If target is in random effects
  if(target %in% ran_trms){
    # check whether there is only a single RE
    if(length(ran_trms) == 1){
      # Fit a lm instead
      fixed_formula <- reformulas::nobars(formula(model))
      fixed_formula <- update(fixed_formula, paste(". ~ . +", ran_trms))
      # Pull out data
      model_data <- model@frame %||% model.frame(model)
      refit_model <- lm(fixed_formula, data = model_data)
    } else {
      ran_frms <- reformulas::findbars(formula(model))
      contains_target <- sapply(ran_frms, function(frm){
          frm <- stringr::str_split(deparse1(frm), " \\| ")[[1]] |> utils::tail(n=1)
        }) == target
      target_ran_frms <- sapply(ran_frms[contains_target], function(frm){
        paste0("(", deparse1(frm), ")")
      })
      target_ran_frms <- paste(target_ran_frms, collapse = "-")
      refit_model <- update(model, as.formula(paste(". ~ . - ", target_ran_frms, " + ", target)))
      check_model_convergence(refit_model)
    }
  } else if (target %in% fixed_trms) { # If target is in fixed effects
    refit_model <- update(model, as.formula(paste(". ~ . + (1|", target, ") - ", target)))
    check_model_convergence(refit_model)
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
#' @noRd
#' @export
print.heritable <- function(x, digits = getOption("digits"), ...) {
  attr(x, "model") <- NULL
  attr(x, "target") <- NULL
  attr(x, "type") <- NULL
  print(unclass(x))
}

#' @keywords internal
build_precompiled_vignette <- function(input = "vignettes/heritable.Rmd.orig",
                                       output = "vignettes/heritable.Rmd") {
  # Build the precompiled vignette
  knitr::knit(here::here(input),
    output = here::here(output)
  )

  # Fix the paths
  lines <- readLines(here::here(output), warn = FALSE)

  # Update the img src and remove any alt attribute
  lines <- gsub(
    'src="vignettes/figs/unnamed-chunk-3-1.png"',
    'src="figs/unnamed-chunk-3-1.png"',
    lines,
    fixed = TRUE
  )
  lines <- gsub(' alt="[^"]+"', "", lines)

  # # Drop caption paragraph lines
  # lines <- lines[!grepl('<p class="caption">', lines, fixed = TRUE)]

  writeLines(lines, here::here(output), useBytes = TRUE)
  invisible(here::here(output))
}

#' @keywords internal
#' @importFrom methods canCoerce hasMethod as
sp2Matrix <- function(x, dense = FALSE, triplet = FALSE) {
  triplet <- ifelse(triplet, yes = "T", no = "C")
  A <- Matrix::sparseMatrix(x[, 1], x[, 2],
    x = x[, 3], repr = triplet,
    symmetric = TRUE
  )
  if (dense) {
    if (canCoerce(A, "packedMatrix")) {
      A <- as(A, "packedMatrix")
    } else if (canCoerce(A, "dsyMatrix") && hasMethod(
      "coerce",
      c("dsyMatrix", "dspMatrix")
    )) {
      A <- as(as(A, "dsyMatrix"), "dspMatrix")
    } else {
      stop("Unable to return a dense matrix")
    }
  }
  if (inherits(x, "ginv")) {
    dimnames(A) <- list(attr(x, "rowNames"), attr(x, "rowNames"))
    attr(A, "INVERSE") <- TRUE
    for (i in c("inbreeding", "logdet", "geneticGroups")) {
      attr(
        A,
        i
      ) <- attr(x, i)
    }
  } else {
    att <- c("rowNames", "INVERSE")
    w <- which(is.element(att, names(attributes(x))))
    if (length(w) > 0) {
      for (i in w) {
        if (i == 1) {
          dimnames(A) <- list(attr(x, "rowNames"), attr(
            x,
            "rowNames"
          ))
        } else {
          attr(A, att[i]) <- attr(x, att[i])
        }
      }
    }
  }
  return(A)
}

#' @noRd
#' @keywords internal
#' @importFrom Matrix Diagonal Matrix t diag
var_diff <- function(V) {
  d <- diag(V)
  delta <- - 2 * V
  delta <- delta + d
  delta <- sweep(delta, 2, d, "+")
  delta
}


#' @noRd
#' @keywords internal
var_comp.lmerMod <- function(model, target, calc_C22 = TRUE,
                             marginal = TRUE, stratification = NULL) {
  X <- lme4::getME(model, "X")
  Z <- lme4::getME(model, "Z")

  sigma2 <- stats::sigma(model)^2
  Lambda <- lme4::getME(model, "Lambda")
  G <- tcrossprod(Lambda) * sigma2
  dimnames(G) <- list(colnames(Z), colnames(Z))

  # Get BLUP weight
  mapper <- map_target_terms(model, target, marginal)
  g <- mapper$idx

  if(is.null(stratification)){
    m <- mapper$m
    intercept <- mapper$intercept
    if(sum(intercept) != 0 && !marginal){
      g <- g[intercept]
      m <- m[intercept, , drop=FALSE]
    }
  } else {
    m <- build_new_Z(model, target, stratification)
  }

  gnames <- levels(model@flist[[target]])
  G_g <- crossprod(m,G[g, g, drop=FALSE]) %*% m
  dimnames(G_g) <- list(gnames, gnames)
  n_g <- length(gnames)

  if(calc_C22){
    R <- diag(nrow(X)) * sigma2
    V <- R + Z %*% G %*% t(Z)
    Vinv <- solve(V)
    P <- Vinv - Vinv %*% X %*% solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv
    C22 <- G - G %*% t(Z) %*% P %*% Z %*% G
    dimnames(C22) <- list(colnames(Z), colnames(Z))
    C22_g <- crossprod(m, C22[g, g, drop=FALSE]) %*% m
    dimnames(C22_g) <- list(gnames, gnames)
  } else {
    C22_g <- NULL
  }

  list(n_g = n_g, G_g = G_g, C22_g = C22_g, gnames = gnames)
}

# To Do, asreml

#' @keywords internal
#' @noRd
var_comp <- function(model, target, calc_C22 = TRUE,
                     marginal = TRUE, stratification = NULL) {
  UseMethod("var_comp")
}
.S3method("var_comp", "lmerMod", var_comp.lmerMod)

#' @noRd
#' @keywords internal
#' @importFrom Matrix colMeans
map_target_terms.lmerMod <- function(model, target, marginal = TRUE){
  mmlist   <- lme4::getME(model, "mmList")
  grp_list <- lme4::getME(model, "flist")
  grp_names <- names(lme4::getME(model, "cnms"))
  Z <- lme4::getME(model, "Z")
  Gp <- lme4::getME(model, "Gp")

  pattern <- paste0("(^|:)", target, "(:|$)")
  matched_grp <- which(grepl(pattern, grp_names))
  idx <- do.call(c,
                 lapply(matched_grp, function(x) (Gp[x]+1):Gp[x+1])
  )
  terms <- colnames(Z[,idx])
  target_grp <- levels(grp_list[[target]])
  n_tg <- length(target_grp)


  w_list <- list()
  m_list <- list()
  intercept_idx <- c()

  for(g in matched_grp){
    mm <- Matrix::Matrix(mmlist[[g]])
    grp <- levels(grp_list[[grp_names[g]]])

    n <- nrow(mm)
    p <- ncol(mm)
    q <- length(grp)

    # Get weighting matrix
    if(grp_names[g] != target){
      target_order <- which(stringr::str_split(grp_names[g], ":", simplify = TRUE) == target)
      grp_no_target <- stringr::str_split(grp, ":", simplify = TRUE)[, -target_order, drop=FALSE]
      grp_no_target <- apply(grp_no_target, 1, function(g) paste0(g, collapse = ":"))
      w <- as.numeric((table(grp_no_target)/q)[grp_no_target])

      # Get intercept terms
      intercept_idx <- c(intercept_idx, rep(0, p*q))
    } else {
      w <- rep(1, p*q)

      # Get intercept terms
      pi <- which("(Intercept)" %in% colnames(mm))
      z <- rep(0, p*q)
      if(length(pi) == 1) z[pi + p * (seq_len(q) - 1)] <- 1
      intercept_idx <- c(intercept_idx, z)
    }
    # BLUP weight
    w <- w * rep(colMeans(mm), q)
    if(!marginal) w <- rep(1, p * q)
    m <- Matrix::Matrix(0, nrow = p * q, ncol = n_tg)

    # BLUP weight matrix
    colnames(m) <- target_grp
    names(w) <- rownames(m) <- rep(grp, each = p)
    for(tg in target_grp){
      pattern <- paste0("(^|:)", tg, "(:|$)")
      m[grepl(pattern, names(w)),tg] <- 1
    }
    m <- m * w
    m_list[[g]] <- m
    w_list[[g]] <- w
  }

  m <- do.call(rbind, m_list)
  w <- do.call(c, w_list)

  intercept <- intercept_idx==1
  list(m = m,
       w = w,
       idx = setNames(idx, terms),
       intercept = setNames(intercept, terms))
}

# To Do, asreml

#' @keywords internal
#' @noRd
map_target_terms <- function(model, target, marginal = TRUE) {
  UseMethod("map_target_terms")
}
.S3method("map_target_terms", "lmerMod", map_target_terms.lmerMod)


#' @keywords internal
#' @noRd
build_new_Z.lmerMod <- function(model, target, new_dat){
  trms <- names(new_dat)
  g <- lme4::getME(model,"flist")[[target]]
  gnames <- levels(g)
  n_g <- nlevels(g)
  new_dat <- matrix(rep(new_dat, n_g), ncol = n_g) |> t() |>
    data.frame()
  colnames(new_dat) <- trms
  new_dat[[target]] <- gnames

  lme4::mkNewReTrms(
    object = model,
    newdata = new_dat,
    re.form = NULL,
    allow.new.levels = TRUE
  )$Z
}

# To Do, asreml

#' @keywords internal
#' @noRd
build_new_Z <- function(model, target, new_dat) {
  UseMethod("build_new_Z")
}
.S3method("build_new_Z", "lmerMod", build_new_Z.lmerMod)

