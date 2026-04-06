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
  ran_trms <- sapply(model$G.param, function(x) {
    facnam <- lapply(x[-1], function(y) y[["facnam"]])
    paste0(do.call(c, facnam), collapse = ":")
  }) |> unname()
  attr(ran_trms, "spec") <- sapply(model$G.param, function(x) {
    facnam <- lapply(x[-1], function(y) y[["model"]])
    paste0(do.call(c, facnam), collapse = ":")
  }) |> unname()

  return(list(fixed = fixed_trms, random = ran_trms))
}

#' @keywords internal
pull_terms.lmerMod <- function(model) {

  fixed_trms <- reformulas::nobars(formula(model))
  if(inherits(fixed_trms,"formula")){
    fixed_trms <- terms(fixed_trms) |>
      attr("term.labels")
  } else {
    fixed_trms <- character()
  }

  mmList <- lme4::getME(model, "mmList")
  ran_trms <- names(mmList)

  attr(ran_trms, "grouping_variable") <- sapply(ran_trms, function(trm) {
    stringr::str_split(trm, " \\| ")[[1]] |> utils::tail(n = 1)
  }, USE.NAMES = FALSE)

  design_trms <- lapply(ran_trms, function(trm) {
    trm <- stringr::str_split(trm, " \\| ")[[1]] |> utils::head(n = 1)
    as.formula(paste("~ ", trm)) |> terms()
  })

  attr(ran_trms, "simple_intercept") <- sapply(design_trms, function(trm) {
    length(attr(trm, "term.labels"))  == 0
  })

  attr(ran_trms, "contain_intercept") <- sapply(design_trms, function(trm) {
    attr(trm, "intercept")  == 1
  })

  return(list(fixed = fixed_trms, random = ran_trms))
}

#' @keywords internal
pull_terms.glmmTMB <- function(model) {
  fixed_trms <- terms(model) |>
    attr("term.labels")

  ran_trms <- glmmTMB::VarCorr(model) |>
    _$cond |>
    names()

  return(list(fixed = fixed_trms, random = ran_trms))
}

#' @keywords internal
pull_terms <- function(model) {
  UseMethod("pull_terms")
}
.S3method("pull_terms", "asreml", pull_terms.asreml)
.S3method("pull_terms", "lmerMod", pull_terms.lmerMod)
.S3method("pull_terms", "glmmTMB", pull_terms.glmmTMB)


#' @keywords internal
pull_terms_without_specials.lmerMod <- function(model) {
  model_terms <- pull_terms(model)
  contain_intercept <- attr(model_terms$random, "contain_intercept")
  model_terms$random <- lapply(
    model_terms$random,
    function(frm) {
      trms <- stringr::str_split(frm, " \\| ")[[1]]
      grp <- utils::tail(trms, 1)
      design <- utils::head(trms, 1)
      design <- as.formula(paste("~ ", design)) |> terms()
      design_terms <- attr(design, "term.labels")
      contain_intercept <- attr(design, "intercept") == 1

      if(length(design_terms) == 0){
        grp
      } else if(contain_intercept){
        c(grp, paste(design_terms, grp, sep = " | "))
      } else {
        paste(design_terms, grp, sep = " | ")
      }
    }
  )
  model_terms$random <- do.call(c,model_terms$random)
  attr(model_terms$random, "grouping_variable") <-
    sapply(model_terms$random, function(trm) {
      stringr::str_split(trm, " \\| ")[[1]] |> utils::tail(n = 1)
    }, USE.NAMES = FALSE)

  model_terms
}

#' @keywords internal
semivariance <- function(X) {
  n <- nrow(X)
  1 / (n - 1) * (sum(diag(X)) - 1 / n * sum(X))
}

#' @keywords internal
pull_terms_without_specials.asreml <- function(model) {
  fixed_trms <- attr(model$formulae$fixed, "term.labels")
  ran_trms <- sapply(model$G.param, function(x) {
    facnam <- lapply(x[-1], function(y) names(y[["model"]]))
    paste0(do.call(c, facnam), collapse = ":")
  }) |> unname()
  attr(ran_trms, "spec") <- sapply(model$G.param, function(x) {
    facnam <- lapply(x[-1], function(y) y[["model"]])
    paste0(do.call(c, facnam), collapse = ":")
  }) |> unname()

  return(list(fixed = fixed_trms, random = ran_trms))
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
pull_terms_without_specials.glmmTMB <- function(model) {
  pull_terms.glmmTMB(model)
}

#' @keywords internal
pull_terms_without_specials <- function(model) {
  UseMethod("pull_terms_without_specials")
}
.S3method("pull_terms_without_specials", "asreml", pull_terms_without_specials.asreml)
.S3method("pull_terms_without_specials", "lmerMod", pull_terms_without_specials.lmerMod)
.S3method("pull_terms_without_specials", "glmmTMB", pull_terms_without_specials.glmmTMB)


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
    )[, 2]
    if (exists(name_GRM, envir = env, inherits = FALSE)) {
      GRM_source <- get(name_GRM, envir = env)
      if (is.data.frame(GRM_source) & ncol(GRM_source) == 3) {
        # TODO: Change in future because GRM_source could be singular
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
fit_counterpart_model.asreml <- function(model, target) {
  # get the terms from model object
  fixed_trms <- pull_terms.asreml(model)$fixed
  ran_trms_without_specials <- pull_terms_without_specials.asreml(model)$random
  ran_trms <- pull_terms.asreml(model)$random

  # when target is in random
  if (target %in% ran_trms_without_specials) {

    target_spcl <- ran_trms[which(ran_trms_without_specials == target)]
    target_ran_frms <- paste(target_spcl, collapse = "-")
    model_counter <- update(model,
      fixed = as.formula(paste(". ~ . +", target)),
      random =  as.formula(paste("~ . -", target_ran_frms)),
      trace = FALSE
    )

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
fit_counterpart_model.lmerMod <- function(model, target) {
  # get the terms from model object
  trms <- pull_terms_without_specials(model)
  fixed_trms <- trms$fixed
  ran_trms <- trms$random

  # If target is in random effects
  if (target %in% ran_trms) {
    # check whether there is only a single RE
    if (length(unique(ran_trms)) == 1) {
      # Fit a lm instead
      fixed_formula <- reformulas::nobars(formula(model))
      fixed_formula <- update(fixed_formula, paste(". ~ . +", ran_trms))
      # Pull out data
      model_data <- model@frame %||% model.frame(model)
      refit_model <- lm(fixed_formula, data = model_data)
    } else {
      ran_trms <- pull_terms(model)$random
      grp_name <- attr(ran_trms, "grouping_variable")
      simple_intercept  <- attr(ran_trms, "simple_intercept")
      contain_intercept <- attr(ran_trms, "contain_intercept")

      # Get interaction terms that contain intercept and use the target as
      # the grouping variable.
      interaction_trm <- !simple_intercept & grp_name == target & attr(ran_trms, "contain_intercept")
      keep_trm <- !(simple_intercept & grp_name == target)
      if(any(interaction_trm)){
        design_trms <- lapply(ran_trms[interaction_trm], function(trm) {
          trm <- stringr::str_split(trm, " \\| ")[[1]] |> utils::head(n = 1)
          as.formula(paste("~ ", trm)) |> terms() |> attr("term.labels")
        })
        new_trm <- sapply(design_trms, function(trm) {
          paste("0 +", paste(trm, collapse = "+"), "|", target)
        })
        ran_trms[interaction_trm] <- new_trm
      }
      frm <- paste(". ~", paste(trms$fixed, collapse = "+"), "+",
             target, "+",
             paste( paste0("(", ran_trms[keep_trm], ")"), collapse = " + ")
      )

      refit_model <- update(model, as.formula(frm))
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
  attr(x, "type") <- NULL
  attr(x, "args") <- NULL
  print(unclass(x))
}

#' Print method for heritable_ci objects
#'
#' @param x An object of class "heritable_ci"
#' @param digits Number of digits to print
#' @param ... Additional arguments passed to print
#' @noRd
#' @export
print.heritable_ci <- function(x, digits = getOption("digits"), ...) {
  attr(x, "boot_mod") <- NULL
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
#' @importFrom Matrix t diag crossprod tcrossprod
var_diff <- function(V) {
  Matrix::Matrix(outer(diag(V), diag(V), "+") - 2 * V)
}


#' @noRd
#' @export
var_comp.lmerMod <- function(model, target,
                             calc_C22 = TRUE, calc_V = TRUE, calc_C11 = TRUE,
                             marginal = TRUE, stratification = NULL,
                             solver = c("direct", "LMM"), ...) {
  solver <- match.arg(solver)
  X <- lme4::getME(model, "X")
  Z <- lme4::getME(model, "Z")

  sigma2 <- stats::sigma(model)^2
  Lambda <- lme4::getME(model, "Lambda")
  G <- tcrossprod(Lambda) * sigma2
  dimnames(G) <- list(colnames(Z), colnames(Z))

  if(marginal || !is.null(stratification)) marginal <- TRUE
  mapper <- map_target_terms(model, target, marginal)
  W <- m <- mapper$m
  m <- m != 0
  g <- mapper$idx

  # Get BLUP weight
  if (is.null(stratification)) {
    main <- mapper$main
    if (sum(!main) == 0) {
      marginal <- FALSE
    } else if (sum(main) != 0 && !marginal) {
      g <- g[main]
      W <- W[main, , drop = FALSE]
      m <- m[main, , drop = FALSE]
      cli::cli_inform(
        c(
          "Heritability will be evaluated using non-interactive components of the target term.",
          "Interactive effects are excluded from this calculation."
        )
      )
      marginal <- FALSE
    } else {
      cli::cli_inform(
        c(
          "Heritability will be computed as a weighted average across interaction components",
          "involving the target term."
        )
      )
      marginal <- TRUE
    }

    active_g <- NULL # Define genotypes on which heritability is measured.
  } else {
    cli::cli_inform(
      c(
        "Heritability will be computed within the user-specified stratification.",
        "Estimation is restricted to the defined strata."
      )
    )
    marginal <- FALSE

    W <- build_new_Z(model, target, stratification)
    active_g <- attr(W, "active")
  }

  gnames <- colnames(W)
  G_g <- crossprod(W, G[g, g, drop = FALSE]) %*% W
  dimnames(G_g) <- list(gnames, gnames)
  n_g <- length(gnames)

  if (calc_V || calc_C22 || calc_C11) {
    R <- diag(nrow(X)) * sigma2
    V <- R + Z %*% G %*% t(Z)

    if(calc_C22) {

      if(solver == "direct"){
        Vinv <- ginv_sym_sparse(V)
        P <- Vinv - Vinv %*% X %*% ginv_sym_sparse(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv
        C22 <- G - G %*% t(Z) %*% P %*% Z %*% G
      }else{
        C22 <- solve_LMM(X, Z, G, R)$C22
      }

      dimnames(C22) <- list(colnames(Z), colnames(Z))
      C22_g <- crossprod(W, C22[g, g, drop = FALSE]) %*% W
      dimnames(C22_g) <- list(gnames, gnames)
    } else {
      C22_g <- NULL
    }

    if(calc_C11){

      if(solver == "direct"){
        Z_g <-  Z[, g, drop=FALSE]
        V_tilde <- V - Z_g %*% G[g,g,drop=FALSE] %*% t(Z_g)
        X_tilde <- cbind(X, Z_g)

        C11 <- ginv_sym_sparse(t(X_tilde) %*% ginv_sym_sparse(V_tilde) %*% X_tilde)
      } else {
        X_tilde <- cbind(X, Z[, g, drop=FALSE])
        if(length(g) != ncol(Z)){
          Z_tilde <-  Z[, -g, drop=FALSE]
          G_tilde <- G[-g,-g,drop=FALSE]
          C11 <- solve_LMM(X_tilde, Z_tilde, G_tilde, R)$C11
        } else {
          C11 <- ginv_sym_sparse(t(X_tilde) %*% ginv_sym_sparse(R) %*% X_tilde)
        }
      }

      XTX <- crossprod(X_tilde)
      P <- ginv_sym_sparse(XTX) %*% XTX
      W_tilde <- t(crossprod(W, P[-seq_len(ncol(X)),-seq_len(ncol(X)),drop=FALSE]))

      # Genetic variance matrix in the conterpart model
      G_g_tilde <- crossprod(W_tilde, G[g, g, drop = FALSE]) %*% W_tilde
      dimnames(G_g_tilde) <- list(gnames, gnames)

      # Genetic variance matrix in the conterpart model without considering genetic covariance
      ij <- apply(m, 2, function(z) {
        active <- which(z)
        expand.grid(active, active)
      }) |>
        do.call(rbind, args = _)
      G_g_no_cov <- G[g, g, drop=FALSE] * Matrix::sparseMatrix(
        i = ij[, 1],
        j = ij[, 2],
        x = 1,
        dims = c(nrow(m), nrow(m))
      )
      G_g_tilde_no_cov <- crossprod(W_tilde, G_g_no_cov) %*% W_tilde
      dimnames(G_g_tilde) <- list(gnames, gnames)

      C11_g <- crossprod(W, C11[-seq_len(ncol(X)),-seq_len(ncol(X)),drop=FALSE]) %*% W
      dimnames(C11_g) <- list(gnames, gnames)

    } else {
      C11_g <- NULL
      G_g_tilde <- NULL
      G_g_tilde_no_cov <- NULL
    }

    if(!calc_V) {
      V <- NULL
      G <- NULL
      Z <- NULL
      X <- NULL
      idx <- NULL
      W <- NULL
    } else {
      idx <- g
    }

  } else {
    V <- NULL
    G <- NULL
    Z <- NULL
    X <- NULL
    idx <- NULL
    W <- NULL
    C22_g <- NULL
    C11_g <- NULL
    G_g_tilde <- NULL
    G_g_tilde_no_cov <- NULL
  }

  if(!is.null(active_g) && any(!active_g)){
    cli::cli_inform(
      c(
        "Following genotypes are not presented in the specified strat: {.code {gnames[!active_g]}}.",
        "They are excluded from the estimation."
      )
    )

    n_g <- sum(active_g)
    gnames <- gnames[active_g]
    G_g <- G_g[active_g,active_g, drop = FALSE]
    if(calc_C22){
      C22_g <- C22_g[active_g,active_g, drop = FALSE]
    }
    if(calc_C11){
      G_g_tilde <- G_g_tilde[active_g,active_g, drop = FALSE]
      G_g_tilde_no_cov <- G_g_tilde_no_cov[active_g,active_g, drop = FALSE]
      C11_g <- C11_g[active_g,active_g, drop = FALSE]
    }
  }

  list(n_g = n_g, gnames = gnames,
       G_g = G_g, C22_g = C22_g,
       G_g_tilde = G_g_tilde, G_g_tilde_no_cov = G_g_tilde_no_cov, C11_g = C11_g,
       V = V, G = G, Z = Z, X = X, idx = idx, W = W,
       marginal = marginal, stratification = stratification)
}

#' @noRd
#' @export
var_comp.asreml <- function(model, target,
                            calc_C22 = TRUE, calc_V = TRUE, calc_C11 = TRUE,
                            marginal = TRUE, stratification = NULL,
                            solver = c("direct", "LMM"),
                            source = list(), ...) {
  solver <- match.arg(solver)
  model <- check_deisgn_exsits(model)
  design <- model$design

  if(marginal || !is.null(stratification)) marginal <- TRUE
  mapper <- map_target_terms(model, target, marginal)
  W <- m <- mapper$m
  m <- m != 0
  g <- mapper$idx

  # Get random terms
  ran_terms <- lapply(model$G.param, function(x) {
    facnam <- lapply(x[-1], function(y) {
      y <- paste0(y[["facnam"]], "_", y[["levels"]])
      factor(y, levels = y)
    })
    facnam <- do.call(
      expand.grid,
      c(rev(facnam), list(stringsAsFactors = FALSE))
    )
    facnam <- facnam[, rev(seq_len(ncol(facnam))), drop = FALSE]
    facnam <- apply(facnam, 1, function(x) paste0(x, collapse = ":"))
  })
  ran_terms <- do.call(c, ran_terms)

  # Build design matrix
  N <- nrow(design)
  Q <- length(ran_terms)
  Z <- Matrix::Matrix(0, N, Q)
  colnames(Z) <- ran_terms
  common_col <- intersect(colnames(design), ran_terms)
  Z[, common_col] <- design[, common_col]
  X <- design[, !colnames(design) %in% ran_terms, drop = FALSE]

  # Get BLUP weight
  if (is.null(stratification)) {
    main <- mapper$main
    if (sum(!main) == 0) {
      marginal <- FALSE
    } else if (sum(main) != 0 && !marginal) {
      g <- g[main]
      W <- W[main, , drop = FALSE]
      m <- m[main, , drop = FALSE]
      cli::cli_inform(
        c(
          "Heritability will be evaluated using non-interactive components of the target term.",
          "Interactive effects are excluded from this calculation."
        )
      )
      marginal <- FALSE
    } else {
      cli::cli_inform(
        c(
          "Heritability will be computed as a weighted average across interaction components",
          "involving the target term."
        )
      )
      marginal <- TRUE
    }
    active_g <- NULL
  } else {
    cli::cli_inform(
      c(
        "Heritability will be computed within the user-specified stratification.",
        "Estimation is restricted to the defined strata."
      )
    )
    marginal <- FALSE
    W <- build_new_Z(model, target, stratification)
    active_g <- attr(W, "active")
  }

  gnames <- colnames(W)
  n_g <- length(gnames)

  if (calc_V || calc_C22 || calc_C11) {
    n_comp <- length(model$G.param)

    # Build variance
    G_list <- lapply(seq_len(n_comp), function(x) {
      phrase_G(model$G.param[[x]], source = source)
    })
    G <- Matrix::bdiag(G_list)
    G <- G * model$sigma2
    dimnames(G) <- list(ran_terms, ran_terms)

    # Identify missing terms.
    col_design <- colnames(design)
    missing_terms <- which(
      Matrix::colSums(design[, intersect(ran_terms, col_design)] != 0) == 0
    ) |> names()
    missing_terms <- unique(c(missing_terms, setdiff(ran_terms, col_design)))
    if (length(missing_terms) > 0) {
      G[missing_terms, ] <- 0
      G[, missing_terms] <- 0
    }

    G_g <- crossprod(W, G[g, g, drop = FALSE]) %*% W
    dimnames(G_g) <- list(gnames, gnames)

    # Build R matrix
    tmp_data_call <- model$call$data
    model$call$data <- Matrix::Matrix(0, nrow = N, ncol = N)
    R <- asremlPlus::estimateV.asreml(model, which.matrix = "R")|>
      Matrix::Matrix()
    model$call$data <- tmp_data_call

    V <- R + Z %*% G %*% t(Z)

    if(calc_C22) {

      if(solver == "direct"){
        Vinv <- ginv_sym_sparse(V)
        P <- Vinv - Vinv %*% X %*% ginv_sym_sparse(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv
        C22 <- G - G %*% t(Z) %*% P %*% Z %*% G
      }else{
        C22 <- solve_LMM(X, Z, G, R)$C22
      }

      dimnames(C22) <- list(colnames(Z), colnames(Z))
      C22_g <- crossprod(W, C22[g, g, drop = FALSE]) %*% W
      dimnames(C22_g) <- list(gnames, gnames)
    } else {
      C22_g <- NULL
    }

    if(calc_C11){

      if(solver == "direct"){
        Z_g <-  Z[, g, drop=FALSE]
        V_tilde <- V - Z_g %*% G[g,g,drop=FALSE] %*% t(Z_g)
        X_tilde <- cbind(X, Z_g)

        C11 <- ginv_sym_sparse(t(X_tilde) %*% ginv_sym_sparse(V_tilde) %*% X_tilde)
      } else {
        X_tilde <- cbind(X, Z[, g, drop=FALSE])
        if(length(g) != ncol(Z)){
          Z_tilde <-  Z[, -g, drop=FALSE]
          G_tilde <- G[-g,-g,drop=FALSE]
          C11 <- solve_LMM(X_tilde, Z_tilde, G_tilde, R)$C11
        } else {
          C11 <- ginv_sym_sparse(t(X_tilde) %*% ginv_sym_sparse(R) %*% X_tilde)
        }
      }

      XTX <- crossprod(X_tilde)
      P <- ginv_sym_sparse(XTX) %*% XTX
      W_tilde <- t(crossprod(W, P[-seq_len(ncol(X)),-seq_len(ncol(X)),drop=FALSE]))

      # Genetic variance matrix in the conterpart model
      G_g_tilde <- crossprod(W_tilde, G[g, g, drop = FALSE]) %*% W_tilde
      dimnames(G_g_tilde) <- list(gnames, gnames)

      # Genetic variance matrix in the conterpart model without considering genetic covariance
      ij <- apply(m, 2, function(z) {
        active <- which(z)
        expand.grid(active, active)
      }) |>
        do.call(rbind, args = _)
      G_g_no_cov <- G[g, g, drop=FALSE] * Matrix::sparseMatrix(
        i = ij[, 1],
        j = ij[, 2],
        x = 1,
        dims = c(nrow(m), nrow(m))
      )
      G_g_tilde_no_cov <- crossprod(W_tilde, G_g_no_cov) %*% W_tilde
      dimnames(G_g_tilde) <- list(gnames, gnames)

      C11_g <- crossprod(W, C11[-seq_len(ncol(X)),-seq_len(ncol(X)),drop=FALSE]) %*% W
      dimnames(C11_g) <- list(gnames, gnames)
    } else {
      C11_g <- NULL
      G_g_tilde <- NULL
      G_g_tilde_no_cov <- NULL
    }

    if(!calc_V) {
      V <- NULL
      G <- NULL
      Z <- NULL
      X <- NULL
      idx <- NULL
      W <- NULL
    } else {
      idx <- g
    }
    # The Z here is different from design, as design may omit terms.
  } else {
    # Build variance
    matched_grp <- sapply(model$G.param, function(x) {
      facnam <- lapply(x[-1], function(y) names(y[["model"]]) == target)
      any(do.call(c, facnam))
    }) |> which()

    G_list <- lapply(matched_grp, function(x) {
      phrase_G(model$G.param[[x]], source = source)
    })
    G <- Matrix::bdiag(G_list)
    G <- G * model$sigma2
    terms <- names(mapper$idx)
    dimnames(G) <- list(terms, terms)

    # Identify missing terms.
    col_design <- colnames(design)
    missing_terms <- which(
      Matrix::colSums(design[, intersect(terms, col_design)] != 0) == 0
    ) |> names()
    missing_terms <- unique(c(missing_terms, setdiff(terms, col_design)))
    if (length(missing_terms) > 0) {
      G[missing_terms, ] <- 0
      G[, missing_terms] <- 0
    }

    G_g <- crossprod(W, G[rownames(W), rownames(W), drop = FALSE]) %*% W
    dimnames(G_g) <- list(gnames, gnames)

    V <- NULL
    G <- NULL
    Z <- NULL
    X <- NULL
    idx <- NULL
    W <- NULL
    C22_g <- NULL
    C11_g <- NULL
    G_g_tilde <- NULL
    G_g_tilde_no_cov <- NULL
  }

  if(!is.null(active_g) && any(!active_g)){
    cli::cli_inform(
      c(
        "Following genotypes are not presented in the specified strat: {.code {gnames[!active_g]}}.",
        "They are excluded from the estimation."
      )
    )

    n_g <- sum(active_g)
    gnames <- gnames[active_g]
    G_g <- G_g[active_g,active_g, drop = FALSE]
    if(calc_C22){
      C22_g <- C22_g[active_g,active_g, drop = FALSE]
    }
    if(calc_C11){
      G_g_tilde <- G_g_tilde[active_g,active_g, drop = FALSE]
      G_g_tilde_no_cov <- G_g_tilde_no_cov[active_g,active_g, drop = FALSE]
      C11_g <- C11_g[active_g,active_g, drop = FALSE]
    }
  }

  list(n_g = n_g, gnames = gnames,
       G_g = G_g, C22_g = C22_g,
       G_g_tilde = G_g_tilde, G_g_tilde_no_cov = G_g_tilde_no_cov, C11_g = C11_g,
       V = V, G = G, Z = Z, X = X, idx = idx, W = W,
       marginal = marginal, stratification = stratification)
}

#' Extract variance components
#'
#' This function provides the variance matrix quantities needed to evaluate heritability
#' with respect to a target random-effect term in a fitted linear mixed model.
#'
#' @param model Model object of class `lmerMod/merMod` or `asreml`
#' @param target The name of the random effect for which heritability is to be calculated.
#' @param calc_C22 Logical; whether to compute the prediction error
#'   variance (PEV) matrix for the transformed target effect.
#' @param calc_V Logical; whether to retain the marginal variance
#'   matrix of the response and supporting quantities (`V`, `G`, `Z`, `X`, `idx`,
#'   and `m`) in the output.
#' @param calc_C11 Logical; whether to compute the variance
#'   matrix of the fixed-effect counterpart estimator for the transformed target
#'   effect.
#' @param marginal Logical; if `TRUE`, construct marginal (strata-averaged)
#'   mappings so that each genotype receives a single averaged effect per term.
#'   If `FALSE`, mappings will only consider the main genotype effect and ignore the
#'   iteracting terms.
#' @param stratification A one-row data frame defining the stratum in which
#'   genotype effects should be evaluated. The columns must correspond
#'   to model terms that interact with `target`.
#' @param solver A string specifying the solver for the PEV matrix. Can be
#' either `"direct"` (directly invert `V`) or `"LMM"` (Solve the LMM equation).
#' @param ... Additional arguments passed to downstream helper functions.
#' @return A named list with the following elements:
#'
#' \describe{
#'   \item{`n_g`}{Number of transformed genetic coefficients after applying the
#'   mapping matrix `m`.}
#'   \item{`gnames`}{Names of the transformed genetic coefficients.}
#'   \item{`G_g`}{Variance matrix of the transformed target genetic effect.}
#'   \item{`C22_g`}{Prediction error variance matrix of the transformed target
#'   effect, if `calc_C22 = TRUE`; otherwise `NULL`.}
#'   \item{`G_g_tilde`}{Variance matrix of the fixed-effect counterpart estimator for
#'   the transformed target effect, if `calc_C11 = TRUE`; otherwise `NULL`.}
#'   \item{`G_g_tilde_no_cov`}{Variance matrix of the fixed-effect counterpart estimator for
#'   the transformed target effect, without considering target covariance, if `calc_C11 = TRUE`; otherwise `NULL`.}
#'   \item{`C11_g`}{Estimation error variance matrix of the fixed-effect counterpart estimator for
#'   the transformed target effect, if `calc_C11 = TRUE`; otherwise `NULL`.}
#'   \item{`V`}{Marginal covariance matrix of the response, if `calc_V = TRUE`;
#'   otherwise `NULL`.}
#'   \item{`G`}{Variance matrix of all random effects, if `calc_V = TRUE`;
#'   otherwise `NULL`.}
#'   \item{`Z`}{Random-effect design matrix, if `calc_V = TRUE`; otherwise `NULL`.}
#'   \item{`X`}{Fixed-effect design matrix, if `calc_V = TRUE`; otherwise `NULL`.}
#'   \item{`idx`}{Indices of the random-effect coefficients associated with the
#'   target term, if `calc_V = TRUE`; otherwise `NULL`.}
#'   \item{`W`}{Linear mapping from the original target coefficients to the
#'   transformed target effect, if `calc_V = TRUE`; otherwise `NULL`.}
#'   \item{`marginal`}{Logical scalar indicating whether the returned quantities
#'   correspond to a marginal definition of the target effect.}
#'   \item{`stratification`}{The user-supplied stratification object, if any.}
#' }
#' @export
var_comp <- function(model, target,
                     calc_C22 = TRUE, calc_V = TRUE, calc_C11 = TRUE,
                     marginal = TRUE, stratification = NULL, solver = c("direct", "LMM"), ...) {
  UseMethod("var_comp")
}
.S3method("var_comp", "lmerMod", var_comp.lmerMod)
.S3method("var_comp", "asreml", var_comp.asreml)

#' @noRd
#' @keywords internal
map_target_terms.lmerMod <- function(model, target, marginal = TRUE) {
  mf <- stats::model.frame(model)
  mmlist <- lme4::getME(model, "mmList")
  grp_list <- lme4::getME(model, "flist")
  cnms <- lme4::getME(model, "cnms")
  grp_names <- names(cnms)
  Z <- lme4::getME(model, "Z")
  Gp <- lme4::getME(model, "Gp")

  # Match mm and Gp
  mmid <- paste(sub(".*\\|\\s*", "", names(mmlist)),
        sapply(mmlist, function(x) paste(colnames(x), collapse = ":")),
        sep = "|"
        )
  Gpid <- paste(
    names(cnms),
    sapply(cnms, function(x) paste(x, collapse = ":")),
    sep = "|"
  )
  mmlist <- mmlist[match(Gpid, mmid)]

  pattern <- paste0("(^|:)", target, "(:|$)")
  matched_grp <- which(grepl(pattern, grp_names))
  idx <- lapply(matched_grp, function(x) (Gp[x] + 1):Gp[x + 1])
  idx_all <- do.call(c, idx)
  terms <- colnames(Z[, idx_all])
  target_grp <- levels(mf[[target]])
  n_tg <- length(target_grp)

  w_list <- list()
  m_list <- list()
  main_idx <- c()

  for (itr in seq_along(matched_grp)) {
    g_idx <- matched_grp[itr]
    mm <- Matrix::Matrix(mmlist[[g_idx]])
    grp <- levels(grp_list[[grp_names[g_idx]]])
    grp <- factor(grp, levels = grp)

    n <- nrow(mm)
    p <- ncol(mm)
    q <- length(grp)

    # Get weighting matrix
    if (grp_names[g_idx] != target) {
      grp_names_split <- stringr::str_split(grp_names[g_idx], ":", simplify = TRUE)
      target_order <- which(stringr::str_split(grp_names_split, ":", simplify = TRUE) == target)
      grp_split <- stringr::str_split(grp, ":", simplify = TRUE)
      target_key <- rep(grp_split[, target_order], each = p)

      if (marginal) {
        grp_no_target <- grp_split[, -target_order, drop = FALSE]
        grp_key <- apply(grp_no_target, 1, paste, collapse = ":")
        mm_key <- colnames(mm)
        stra_key <- apply(expand.grid(mm_key, grp_key), 1, function(x) paste0(x, collapse = ":"))
        stra_id <- match(stra_key, unique(stra_key))

        z <- Z[, idx[[itr]]]
        w <- numeric(p * q)

        for (id in unique(stra_id)) {
          s <- which(stra_id == id)

          if (length(s) > 0) {
            w[s] <- sum(z[, s]) / n
          }
        }
      } else {
        w <- rep(1, p * q)
      }

      # Get main terms
      main_idx <- c(main_idx, rep(0, p * q))
    } else {
      # BLUP weight
      if (marginal) {
        w <- rep(Matrix::colMeans(mm), q)
      } else {
        w <- rep(1, p * q)
      }
      target_key <- rep(target_grp, each = p)

      # Get main terms
      pi <- which("(Intercept)" %in% colnames(mm))
      z <- rep(0, p * q)
      if (length(pi) == 1) z[pi + p * (seq_len(q) - 1)] <- 1
      main_idx <- c(main_idx, z)
    }

    grp_key <- rep(grp, p)
    names(w) <- grp_key
    m <- build_f_mat(target_key, target_grp)
    dimnames(m) <- list(grp_key, target_grp)

    m <- m * w
    m_list[[itr]] <- m
    w_list[[itr]] <- w
  }

  m <- do.call(rbind, m_list)
  w <- do.call(c, w_list)

  main <- main_idx == 1
  list(
    m = m,
    w = w,
    idx = setNames(idx_all, terms),
    main = setNames(main, terms)
  )
}

#' @keywords internal
#' @noRd
map_target_terms.asreml <- function(model, target, marginal = TRUE) {
  model <- check_deisgn_exsits(model)
  design <- model$design
  mf <- model$mf

  grp_names <- sapply(model$G.param, function(x) {
    facnam <- lapply(x[-1], function(y) names(y[["model"]]))
    paste0(do.call(c, facnam), collapse = ":")
  })

  matched_grp <- sapply(model$G.param, function(x) {
    facnam <- lapply(x[-1], function(y) names(y[["model"]]) == target)
    any(do.call(c, facnam))
  }) |> which()

  grp_terms <- lapply(model$G.param, function(x) {
    facnam <- lapply(x[-1], function(y) {
      y <- paste0(y[["facnam"]], "_", y[["levels"]])
      factor(y, levels = y)
    })
    facnam <- do.call(
      expand.grid,
      c(rev(facnam), list(stringsAsFactors = FALSE))
    )
    facnam <- facnam[, rev(seq_len(ncol(facnam))), drop = FALSE]
    facnam
  })
  Gp <- cumsum(sapply(grp_terms, nrow))
  Gp <- c(0, Gp)
  idx_list <- lapply(matched_grp, function(x) (Gp[x] + 1):Gp[x + 1])
  idx <- do.call(c, idx_list)
  terms_list <- lapply(grp_terms, function(x) {
    apply(x, 1, function(y) paste0(y, collapse = ":"))
  })
  terms <- do.call(c, terms_list)
  terms <- terms[idx]

  g <- mf[[target]]
  if (!is.factor(g)) {
    cli::cli_abort("The target {.code {target}} provided is not a factor.")
  }
  target_grp <- levels(g)
  n_tg <- nlevels(g)
  n <- nrow(mf)

  w_list <- list()
  m_list <- list()
  main_idx <- c()


  for (i in seq_along(matched_grp)) {
    g_idx <- matched_grp[i]
    grp_name <- grp_names[g_idx]
    vars <- sapply(model$G.param[[g_idx]][-1], function(x) names(x[["model"]]))
    target_idx <- which(vars == target)
    stra_idx <- which(vars != target)
    grp <- terms_list[[g_idx]]
    n_g <- length(grp)

    grp_split <- grp_terms[[g_idx]]
    target_key <- grp_split[, target_idx]
    target_lev <- levels(target_key)

    if (grp_name != target) {
      if (marginal) {
        stra_key <- apply(grp_split[, stra_idx, drop = FALSE], 1, function(x) paste0(x, collapse = ":"))
        stra_id <- match(stra_key, unique(stra_key))
        stra_colname <- lapply(
          unique(stra_id),
          function(id) grp[stra_id == id]
        )

        stra_w <- sapply(stra_colname, function(s) {
          s <- intersect(colnames(design), s)
          sum(design[, s, drop = FALSE]) / n
        })
        w <- stra_w[stra_id]
      } else {
        w <- rep(1, n_g)
      }

      # Get main terms
      main_idx <- c(main_idx, rep(0, n_g))
    } else {
      w <- rep(1, n_g)
      main_idx <- c(main_idx, rep(1, n_g))
    }

    m <- build_f_mat(target_key, target_lev)
    m <- m * w
    dimnames(m) <- list(grp, target_grp)
    names(w) <- grp
    m_list[[i]] <- m
    w_list[[i]] <- w
  }

  m <- do.call(rbind, m_list)
  w <- do.call(c, w_list)

  main <- main_idx == 1
  list(
    m = m,
    w = w,
    idx = setNames(idx, terms),
    main = setNames(main, terms)
  )
}


#' Map target-associated random-effect terms to genotype-level effects
#'
#' Identify random-effect terms in a mixed model that involve a target random factor
#' (e.g. a genotype effect) and construct linear mapping matrices that summarise
#' those term-level random effects as genotype-level effects.
#'
#' In models where the target factor appears in multiple random-effect terms
#' (e.g. `gen`, `gen:env`), the fitted random effects live on different
#' coefficient spaces (one per term). This function determines which terms are
#' associated with `target` and builds weighting / aggregation matrices that map
#' each term’s coefficient vector onto a common genotype-level scale (typically an
#' average over the interacting strata).
#'
#' @param model A fitted mixed model object.
#' @param target A character string giving the name of the target random-effect
#'   factor.
#' @param marginal Logical; if `TRUE`, construct marginal (strata-averaged)
#'   mappings so that each genotype receives a single averaged effect per term.
#'   If `FALSE`, mappings will only consider the main genotype effect and ignore the
#'   iteracting terms.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{m}{A numeric/sparse matrix whose rows corresponds to
#'     random-effect coefficients across all matched terms, and columns matche the levels
#'     of `target`.
#'   \item{w}{A numeric vector of giving the per-coefficient weights
#'     applied in `m` (i.e. `m` is formed by multiplying an indicator map by `w`).}
#'   \item{idx}{An integer vector giving the indices of the random effect
#'     variance-covariance matrix that correspond to the matched target-associated coefficients.
#'   \item{main}{A logical vector indicating whether each
#'     matched coefficient corresponds columns for the main `target`
#'     random-effect term. Named consistently with `idx`.}
#' }
#'
#' @keywords internal
#' @noRd
map_target_terms <- function(model, target, marginal = TRUE) {
  UseMethod("map_target_terms")
}
.S3method("map_target_terms", "lmerMod", map_target_terms.lmerMod)
.S3method("map_target_terms", "asreml", map_target_terms.asreml)


#' @keywords internal
#' @noRd
build_new_Z.lmerMod <- function(model, target, new_data) {
  trms <- colnames(new_data)
  mf <- stats::model.frame(model)

  for (trm in trms) {
    if (!trm %in% colnames(mf)) {
      cli::cli_abort("{.code {trm}} in {.code new_data} was not found in the model.")
    }

    if (!is.factor(mf[[trm]]) && !inherits(new_data[, trm], "numeric")) {
      cli::cli_abort("{.code {trm}} should be a numeric.")
    }

    if (is.factor(mf[[trm]]) && !new_data[, trm] %in% levels(mf[[trm]])) {
      cli::cli_abort("Unknow level in {.code {trm}} detected: {.code {new_data[trm]}}.")
    }
  }

  # Target levels
  g <- mf[[target]]
  gnames <- levels(g)
  n_g <- nlevels(g)

  if (length(trms) == 1) {
    new_data <- data.frame(rep(new_data[, 1], n_g))
    colnames(new_data) <- trms
  } else {
    new_data <- new_data[rep(1, n_g), , drop = FALSE]
  }
  new_data[[target]] <- factor(gnames, levels = gnames)

  # Add factor level and contrast
  for (trm in trms) {
    if (is.factor(mf[[trm]])) {
      new_data[[trm]] <- factor(new_data[[trm]], levels = levels(mf[[trm]]))
      stats::contrasts(new_data[[trm]]) <- stats::contrasts(mf[[trm]])
    }
  }

  mmlist <- lme4::getME(model, "mmList")
  mm_names <- names(mmlist)
  frm <- stringr::str_extract(mm_names, "^.+(?=\\ \\|)")
  frm <- paste("~", frm)

  grp_list <- lme4::getME(model, "flist")
  cnms <- lme4::getME(model, "cnms")
  grp_names <- names(cnms)

  # Match mm and Gp
  mmid <- paste(sub(".*\\|\\s*", "", names(mmlist)),
                sapply(mmlist, function(x) paste(colnames(x), collapse = ":")),
                sep = "|"
  )
  Gpid <- paste(
    names(cnms),
    sapply(cnms, function(x) paste(x, collapse = ":")),
    sep = "|"
  )
  mmlist <- mmlist[match(Gpid, mmid)]

  pattern <- paste0("(^|:)", target, "(:|$)")
  matched_grp <- which(grepl(pattern, grp_names))

  # Split grouping-variable names per matched term, e.g. "gen:rep" -> c("gen","rep")
  grp_names_split <- stringr::str_split(grp_names, ":")
  design_names_split <- lapply(mm_names, function(trm) {
    trm <- stringr::str_split(trm, " \\| ")[[1]] |> utils::head(n = 1)
    as.formula(paste("~ ", trm)) |> terms() |> attr("term.labels")
  })
  required_var <- do.call(c, grp_names_split[matched_grp]) |> unique()
  required_var <- c(required_var, do.call(c, design_names_split[matched_grp])) |> unique()
  missing_trms <- required_var[!required_var %in% trms]
  missing_trms <- missing_trms[missing_trms != target]
  if (length(required_var) == 1) {
    cli::cli_abort("No stratification needed.")
  }
  if (length(missing_trms) > 0) {
    cli::cli_abort("Terms {.code {missing_trms}} interact with {.code {target}} but were not provided.")
  }

  Z_list <- list()
  active_g <- rep(TRUE, n_g)

  for (itr in seq_along(matched_grp)) {
    g_idx <- matched_grp[itr]
    term <- grp_names[g_idx]

    # Within-group design for this RE term (n x p)
    mm <- Matrix::sparse.model.matrix(stats::as.formula(frm[g_idx]), new_data)
    grp <- levels(grp_list[[term]])
    grp <- factor(grp, levels = grp)

    n <- nrow(mm)
    p <- ncol(mm)
    q <- length(grp)
    grp_new <- apply(new_data[, grp_names_split[[g_idx]], drop = FALSE], 1, paste, collapse = ":")
    mm_grp <- build_f_mat(grp_new, grp)

    # Check which factor levels are missing
    active_g <- active_g & (grp_new %in% grp)

    z <- Matrix::KhatriRao(t(mm_grp), t(mm)) |> t()
    dimnames(z) <- list(gnames, rep(grp, each = p))
    Z_list[[itr]] <- z
  }
  Z <- do.call(cbind, Z_list) |> t()
  attr(Z, "active") <- active_g
  Z
}

#' @keywords internal
#' @noRd
build_new_Z.asreml <- function(model, target, new_data) {
  model <- check_deisgn_exsits(model)
  design <- model$design
  mf <- model$mf

  grp_names <- sapply(model$G.param, function(x) {
    facnam <- lapply(x[-1], function(y) names(y[["model"]]))
    paste0(do.call(c, facnam), collapse = ":")
  })

  matched_grp <- sapply(model$G.param, function(x) {
    facnam <- lapply(x[-1], function(y) names(y[["model"]]) == target)
    any(do.call(c, facnam))
  }) |> which()

  required_var <- lapply(model$G.param, function(x) {
    facnam <- lapply(x[-1], function(y) names(y[["model"]]))
    unique(do.call(c, facnam))
  })
  required_var <- do.call(c, required_var[matched_grp]) |> unique()

  if (length(required_var) == 1) {
    cli::cli_abort("No stratification needed.")
  }

  # Resolve terms
  grp_terms <- lapply(model$G.param, function(x) {
    facnam <- lapply(x[-1], function(y) {
      y <- paste0(y[["facnam"]], "_", y[["levels"]])
      factor(y, levels = y)
    })
    facnam <- do.call(
      expand.grid,
      c(rev(facnam), list(stringsAsFactors = FALSE))
    )
    facnam <- facnam[, rev(seq_len(ncol(facnam))), drop = FALSE]
    facnam
  })
  Gp <- cumsum(sapply(grp_terms, nrow))
  Gp <- c(0, Gp)
  idx_list <- lapply(matched_grp, function(x) (Gp[x] + 1):Gp[x + 1])
  idx <- do.call(c, idx_list)
  terms_list <- lapply(grp_terms, function(x) {
    apply(x, 1, function(y) paste0(y, collapse = ":"))
  })
  terms <- do.call(c, terms_list)
  terms <- terms[idx]

  trms <- colnames(new_data)

  for (trm in trms) {
    if (!trm %in% colnames(mf)) {
      cli::cli_abort("{.code {trm}} in {.code new_data} was not found in the model.")
    }

    if (!is.factor(mf[[trm]]) && !inherits(new_data[, trm], "numeric")) {
      cli::cli_abort("{.code {trm}} should be a numeric.")
    }

    if (is.factor(mf[[trm]]) && !new_data[, trm] %in% levels(mf[[trm]])) {
      cli::cli_abort("Unknow level in {.code {trm}} detected: {.code {new_data[trm]}}.")
    }
  }

  missing_trms <- required_var[!required_var %in% trms]
  missing_trms <- missing_trms[missing_trms != target]
  if (length(missing_trms) > 0) {
    cli::cli_abort("Terms {.code {missing_trms}} interact with {.code {target}} but were not provided.")
  }

  # Target levels
  g <- mf[[target]]
  if (!is.factor(g)) {
    cli::cli_abort("The target {.code {target}} provided is not a factor.")
  }
  gnames <- levels(g)
  n_g <- nlevels(g)

  if (length(trms) == 1) {
    new_data <- data.frame(rep(new_data[, 1], n_g))
    colnames(new_data) <- trms
  } else {
    new_data <- new_data[rep(1, n_g), , drop = FALSE]
  }
  new_data[[target]] <- factor(gnames, levels = gnames)


  # Add factor level and contrast
  var_type <- c()
  for (trm in trms) {
    var_type[trm] <- class(mf[[trm]])
    if (is.factor(mf[[trm]])) {
      new_data[[trm]] <- factor(new_data[[trm]], levels = levels(mf[[trm]]))
      stats::contrasts(new_data[[trm]]) <- stats::contrasts(mf[[trm]])
    }
  }
  numeric_var <- trms[var_type == "numeric"]

  Z_list <- list()

  active_g <- rep(TRUE, n_g)

  for (i in seq_along(matched_grp)) {
    g_idx <- matched_grp[i]
    term <- grp_names[g_idx]
    vars <- sapply(model$G.param[[g_idx]][-1], function(x) names(x[["model"]]))
    specs <- sapply(model$G.param[[g_idx]][-1], function(x) x[["model"]])
    grp <- terms_list[[g_idx]]

    numeric_var_idx <- which(vars %in% numeric_var)
    factor_idx <- which(!vars %in% numeric_var)
    stra_idx <- which(!vars %in% numeric_var & vars != target)
    target_idx <- which(vars == target)

    grp_split <- grp_terms[[g_idx]]
    target_key <- grp_split[, target_idx]
    target_levs <- levels(target_key)

    # Check which factor levels are missing
    factor_name <- vars[factor_idx]
    available_factor_levs <- apply(mf[,factor_name, drop=FALSE],1, function(x) paste(x, collapse = ":"))
    required_factor_levs <- apply(new_data[,factor_name, drop=FALSE],1, function(x) paste(x, collapse = ":"))
    active_g <- active_g & (required_factor_levs %in% available_factor_levs)

    if (length(numeric_var_idx) > 0 && any(specs[numeric_var_idx] != "id")) {
      cli::cli_inform(
        c(
          "The target interacts with numerical predictors generated by special model syntax:",
          "  {.code {specs[numeric_var_idx]}}",
          "These predictors will be projected into the new space using a nearest-neighbour approximation."
        )
      )

      # Helper function, build basis
      build_basis <- function(Z, stra) {
        sapply(stra, function(s) {
          s <- intersect(colnames(Z), s)
          Matrix::rowSums(Z[, s, drop = FALSE])
        })
      }

      # Group by numeric terms
      numeric_var_key <- apply(grp_split[, numeric_var_idx, drop = FALSE], 1, function(x) paste0(x, collapse = ":"))
      numeric_var_id <- match(numeric_var_key, unique(numeric_var_key))
      numeric_var_colname <- lapply(
        unique(numeric_var_id),
        function(id) grp[numeric_var_id == id]
      )

      # Choose best approximation row in mf, then get weights per group
      best_apprx <- sweep(
        data.frame(mf)[, vars[numeric_var_idx], drop = FALSE],
        2,
        new_data[1, vars[numeric_var_idx]],
        "-"
      )^2 |>
        rowSums() |>
        which.min()

      z <- build_basis(design, numeric_var_colname)[best_apprx, numeric_var_id]

      if (length(stra_idx) > 1) {
        new_strata <- new_data[1, vars[stra_idx]]
        stra_key <- apply(grp_split[, stra_idx, drop = FALSE], 1, function(x) paste0(x, collapse = ":"))
        z <- z * (stra_key == new_strata)
      }
      Z_list[[i]] <- t(z * build_f_mat(target_levs, target_levs))
    } else {
      term <- stringr::str_split(term, ":")[[1]]
      term <- paste0(rev(term), collapse = ":")
      frm <- as.formula(paste0("~0+", term))
      Z_list[[i]] <- Matrix::sparse.model.matrix(frm, new_data)
    }
  }

  Z <- do.call(cbind, Z_list) |> t()
  rownames(Z) <- unname(terms)
  colnames(Z) <- gnames

  attr(Z, "active") <- active_g

  Z
}

#' Construct a strata-specific random-effects design matrix
#'
#' Build a new random-effects design matrix that maps genotype effects
#' within a specified stratum for a fitted mixed model object.
#'
#' For each level of the target random term (e.g. a genotype factor),
#' this function constructs a design matrix that activates the effect
#' corresponding to the supplied `new_data` stratum, while setting all
#' other strata-specific contributions to zero. This enables prediction
#' or extraction of genotype effects within a particular interacting
#' context (e.g. genotype × environment).
#'
#' The returned matrix preserves the ordering and naming of random-effect
#' coefficients as defined by the original mixed-model specification.
#'
#' @param model A fitted mixed model object.
#' @param target A character string giving the name of the target random-effect
#'   factor.
#' @param new_data A one-row data frame defining the stratum in which
#'   genotype effects should be evaluated. The columns must correspond
#'   to model terms that interact with `target`.
#'
#' @return A sparse matrix whose rows correspond to
#'   levels of the target factor and whose columns align with the
#'   original random-effect coefficients involving `target`.
#'
#' @keywords internal
#' @noRd
build_new_Z <- function(model, target, new_data) {
  UseMethod("build_new_Z")
}
.S3method("build_new_Z", "lmerMod", build_new_Z.lmerMod)
.S3method("build_new_Z", "asreml", build_new_Z.asreml)

#' Helper function, build factor matrix
#' @keywords internal
#' @noRd
build_f_mat <- function(x, level) {
  i <- seq_along(x)
  j <- match(x, level)
  keep <- !is.na(j)
  mm_grp <- Matrix::sparseMatrix(
    i    = i[keep],
    j    = j[keep],
    x    = rep(1, sum(keep)),
    dims = c(length(x), length(level))
  )
  mm_grp
}

#' @keywords internal
#' @noRd
ginv_sym_sparse <- function(A, tol = .Machine$double.eps,
                            exact_psd_inv = getOption("exact_psd_inv", TRUE)) {

  if(!exact_psd_inv) {
    A <- Matrix::forceSymmetric(A)
    n <- nrow(A)
    I <- Matrix::Diagonal(n)

    if(n == 1) return(1/A)

    smax <- irlba::irlba(A, nv = 1, nu = 1)$d[1]
    lambda <- tol * smax^2 + tol

    Matrix::solve(
      Matrix::crossprod(A) + lambda * I,
      Matrix::t(A), tol = -Inf
    )
  } else {
    A <- Matrix::forceSymmetric(A)

    e <- eigen(as.matrix(A), symmetric = TRUE)

    d_inv <- ifelse(e$values > sqrt(.Machine$double.eps), 1 / e$values, 0)

    e$vectors %*% (d_inv * t(e$vectors))
  }

}


#' @keywords internal
#' @noRd
solve_LMM <- function(X, Z, G, R){
  zero_idx <- Matrix::diag(G) == 0
  keep_idx <- !zero_idx

  Z_reduce <- Z[, keep_idx, drop = FALSE]
  G_reduce <- G[keep_idx, keep_idx, drop = FALSE]

  Rinv <- ginv_sym_sparse(R)
  Ginv <- ginv_sym_sparse(G_reduce)

  XtRinvX <- crossprod(X, Rinv) %*% X
  XtRinvZ <- crossprod(X, Rinv) %*% Z_reduce
  ZtRinvZ <- crossprod(Z_reduce, Rinv) %*% Z_reduce

  C <- rbind(
    cbind(XtRinvX, XtRinvZ),
    cbind(t(XtRinvZ), ZtRinvZ + Ginv)
  )

  Cinv <- ginv_sym_sparse(Matrix::forceSymmetric(C))

  p <- ncol(X)
  q <- ncol(Z)

  C11 <- Cinv[seq_len(p), seq_len(p), drop = FALSE]

  C22_reduce <- Cinv[p + seq_len(ncol(Z_reduce)),
                     p + seq_len(ncol(Z_reduce)),
                     drop = FALSE]

  C22 <- Matrix::Matrix(0, nrow = q, ncol = q, dimnames = dimnames(G))
  C22[keep_idx, keep_idx] <- C22_reduce
  C22 <- Matrix::forceSymmetric(C22)

  dimnames(C22) <- list(colnames(Z), colnames(Z))

  list(C11 = C11,
       C22 = C22)
}
