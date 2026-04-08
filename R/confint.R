#' Bootstrap confidence interval for heritability
#'
#' @description
#' Computes a confidence interval for a heritability estimate using parametric
#' bootstrap of the underlying mixed model. Parallel computing is supported through
#' [boot::boot()]
#'
#' @param object
#' A heritability object returned by [H2()] (broad-sense) [h2()] (narrow-sense). The object must store the fitted model
#' as an attribute.
#' @param parm a specification of which methods are to be given confidence intervals,
#' either a vector of numbers or a vector of names.
#' If missing, all methods are considered.
#' @param level Confidence level.
#' @param B Integer. Number of bootstrap replicates.
#' @param random_effect Character. Strategy for handling random effects.
#'   \describe{
#'     \item{`"resample"`}{Resample random effects to propagate uncertainty.}
#'     \item{`"conditional"`}{Condition on estimated random effects.}
#'   }
#' @param type Character. Bootstrap interval type; one of `"basic"`,
#' `"norm"`, or `"perc"`.
#' @param return_model Logical. Whether to return to the `boot` object.
#' @param seed Optional random seed.
#' @param ... Additional arguments passed to the bootstrap routine. Check
#' [boot::boot()], as well as the examples below for parallel computation
#' @return
#' A matrix of confidence intervals.
#' @examples
#' \dontrun{
#'   lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
#'   lettuce_asreml <- asreml(
#'     fixed = y ~ rep * pseudo_var1,
#'     random = ~gen,
#'     sparse = ~pseudo_var2,
#'     data = lettuce_subset,
#'     trace = FALSE
#'   )
#'
#'   my_H2 <- H2(lettuce_asreml, "gen", c("Cullis", "Standard"))
#'
#'   my_ci <- confint(my_H2)
#'
#'   # Get bootstrap model
#'   boot_mod <- attr(my_ci, "boot_mod")
#'   boot_mod$t # Check bootstrap statistics
#'
#'   # Parallel computing (On Windows)
#'   # Note that for asreml, `ncpus` can't be larger than the number of asreml
#'   # license available
#'   confint(my_H2, parallel = "snow", ncpus = 3)
#'
#'   # Parallel computing (On non-Windows)
#'   confint(my_H2, parallel = "multicore", ncpus = 3)
#'
#' }
#' @seealso [boot::boot()], [H2()]
#' @export
confint.heritable <- function(object,
                              parm = NULL,
                              level = 0.95,
                              B = 100,
                              random_effect = c("resample", "conditional"),
                              type = c("basic", "norm", "perc"),
                              return_model = TRUE,
                              seed = NULL,
                              ...) {

  parallel <- list(...)[["parallel"]]
  if(!is.null(parallel) && parallel != "no"){
    if(inherits(model, "lmerMod")) require(lme4)
    if(inherits(model, "asreml")) require(asreml)
  }

  # basic: bias corrected percentile interval
  # norm: bias corrected normal interval
  # perc: percentile interval
  type <- match.arg(type)
  random_effect <- match.arg(random_effect)

  if(random_effect == "conditional"){
    resample_u <- FALSE
  } else {
    resample_u <- TRUE
  }

  args <- attr(object, "args")
  args[["options"]] <- list(check = FALSE)
  model <- args[["model"]]
  all_method <- args[["method"]]
  all_method <- all_method[!is.na(object)]

  if(is.null(parm)){
    method <- all_method
  } else {
    method <- intersect(parm, all_method)
    if(length(method) == 0){
      cli::cli_abort("Invalid choice of method, should be one of {.code {all_method}}")
    }
  }

  args[["method"]] <- method

  h2.type <- attr(object, "type")  # Get the heritability type

  if (h2.type == "broad_sense") {
    Fun_use <- function(x) {
      args[["model"]] <- x
      do.call(H2, args)
    }
  } else {
    Fun_use <- function(x) {
      args[["model"]] <- x
      do.call(h2, args)
    }
  }

  if (inherits(model, "lmerMod")) {
    suppressMessages(
      boot_mod <- lme4::bootMer(model,
                                FUN = Fun_use, nsim = B, seed = seed,
                                use.u = !resample_u, ...
      )
    )
    ci <- stats::confint(boot_mod, level = level, type = type)
  } else {
    suppressMessages(
      boot_mod <- bootstrap_asreml(model,
        FUN = Fun_use, nsim = B, seed = seed,
        use.u = !resample_u, source = args$source, ...
      )
    )
    ci <- stats::confint(boot_mod, level = level, type = type)
    ci <- matrix(ci,
                 nrow = length(method),
                 dimnames = list(method, colnames(ci))
          )
  }

  structure(ci,
            class = c("heritable_ci", class(ci)),
            boot_mod = boot_mod
  )

}

#' Parametric bootstrap for an asreml model.
#'
#' @description
#' Simulate \eqn{\hat{y} \sim N(X\hat{\beta}, V)} according to the current asreml fit
#' and then refit to obtain the targeted statistics.
#'
#' @param model An `asreml` fitted model. Must be fitted with `model.frame = TRUE`.
#' @param nsim Integer. Number of bootstrap replicates.
#' @param FUN A function with signature `function(fit)` returning a scalar
#'   (the statistic to bootstrap).
#' @param use.u A logical indicating whether to resample random effects, or only
#' resample residuals.
#' @param source The known genomic relationship matrix (GRM) used in `model` fitted using `asreml::vm()`, provided as a named list.
#' When not provided (an empty list by default), the GRM variable used for `vm` calling will be searched in the global environment.
#' @param seed Optional integer seed for reproducibility.
#' @param ... Additional arguments passed to [boot::boot()].
#'
#' @return A `boot` object.
#'
#' @details
#' Fits parametric bootstrap replicates for an `asreml` model by:
#' - Extracting the fixed-effect fit yhat = X * beta.
#' - Extracting V = Var(y) on the observation scale,
#' - Simulating new responses y* = yhat + L %*% z where L is a Cholesky factor of V,
#' - Refitting the same asreml call on each simulated dataset,
#' - Returning a `boot` object.
#'
#' @examples
#' \dontrun{
#'   lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
#'   lettuce_asreml <- asreml(
#'     fixed = y ~ rep * pseudo_var1,
#'     random = ~gen,
#'     sparse = ~pseudo_var2,
#'     data = lettuce_subset,
#'     trace = FALSE
#'   )
#'
#'   b <- bootstrap_asreml(
#'     lettuce_asreml,
#'     R = 200,
#'     statistic = function(fit) coef(fit)$fixed["(Intercept)", "effect"],
#'     seed = 1
#'   )
#'   boot::boot.ci(b, type = "perc")
#'  }
#' @export
bootstrap_asreml <- function(model,
                             FUN,
                             nsim = 1,
                             use.u = FALSE,
                             source = list(),
                             seed = NULL,
                             ...) {
  if (!inherits(model, "asreml")) {
    stop("`model` must be an `asreml` object.")
  }

  # Get model frame
  design_default <- asreml::asreml.options()$design
  asreml::asreml.options(design = TRUE)
  model <- check_deisgn_exsits(model)
  mf <- model$mf

  mf <- as.data.frame(mf)
  N <- nrow(mf)

  # Get response name if not provided.
  fixed_formula <- model$formulae$fixed
  resp_idx <- attr(terms(fixed_formula), "response")
  response <- all.vars(fixed_formula)[resp_idx]

  # Get the linear predictor of fixed effect.
  if (!use.u) {
    yhat <- get_fixed_fit_asreml(model)
    V <- tryCatch(
      asremlPlus::estimateV.asreml(model, which.matrix = "V"),
      error = function(e) {
        cli::cli_warn(
          c("estimateV.asreml failed: {.code {e$message}}",
                      "Try custom extraction of V due to the failure of estimateV.asreml.")
        )
        NULL
      }
    )

    if (is.null(V)) {
      get_V <- function(model, source){
        grp_names <- sapply(model$G.param, function(x) {
          facnam <- lapply(x[-1], function(y) y[["facnam"]])
          paste0(do.call(c, facnam), collapse = ":")
        })

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

        # Build variance
        G_list <- lapply(seq_along(grp_names), function(x) {
          phrase_G(model$G.param[[x]], source = source)
        })
        G <- Matrix::bdiag(G_list)
        G <- G * model$sigma2
        Q <- ncol(G)
        design <- model$design
        N <- nrow(design)
        Z <- Matrix::Matrix(0, N, Q)
        colnames(Z) <- ran_terms
        common_col <- intersect(colnames(design), ran_terms)
        Z[, common_col] <- design[, common_col]

        R <- asremlPlus::estimateV.asreml(model, which.matrix = "R")|>
          Matrix::Matrix()

        V <- R + Z %*% G %*% t(Z)
      }

      source <- check_GRM_exists(model = model, source = source, return = TRUE)

      V <- tryCatch(
        get_V(model, source),
        error = function(e) {
          cli::cli_warn("Failed to extract V: {.code {e$message}}")
          NULL
        }
      )
    }

    if (is.null(V)) {
      cli::cli_warn("Random effects will not be resampled.")
      use.u <- TRUE
    }
  }

  if (use.u) {
    yhat <- model$linear.predictors
    V <- asremlPlus::estimateV.asreml(model, which.matrix = "R")|>
      Matrix::Matrix()
  }

  # Get V and Cholesky factor L such that L %*% z ~ N(0, V)
  # Use Matrix chol: returns upper triangular U with U' U = V
  # So L = t(U) gives L L' = V
  L <- Matrix::chol(V) |> t()

  # Bootstrap data
  boot_data <- mf
  boot_data[[".yhat"]] <- yhat

  # generator: simulate response into correct column name
  generate_data <- function(data, mle) {
    out <- data
    N <- nrow(out)
    eps <- as.numeric(L %*% stats::rnorm(N))
    out[[response]] <- out[[".yhat"]] + eps
    out
  }

  # Refit wrapper for boot()
  refit_asreml <- function(data, model, FUN, source, design) {
    if(length(source) > 0){
      list2env(source, environment())
    }
    fit <- asreml::update.asreml(model, data = data)
    fit[["design"]] <- design
    FUN(fit)
  }

  if (!is.null(seed)) set.seed(seed)
  boot <- boot::boot(
    data      = boot_data,
    statistic = refit_asreml,
    sim       = "parametric",
    ran.gen   = generate_data,
    R         = nsim,
    model     = model,
    source    = source,
    design    = model$design,
    mle       = L,
    FUN       = FUN,
    ...
  )

  # Set the asreml option back
  asreml::asreml.options(design = design_default)

  boot

}


#' Fixed-effects-only fitted values from an asreml model
#'
#' @description
#' Returns the fitted values based on fixed effects only
#' (\eqn{\hat{y} = X\hat{\beta}}), excluding all random effects.
#'
#' @param model An object of class \code{"asreml"}, fitted with
#' \code{model.frame = TRUE}.
#'
#' @return
#' A numeric vector of length \eqn{N}, giving the fixed-effects-only
#' fitted value for each observation.
#'
#' @details
#' This function reconstructs the fixed-effect design matrix from the
#' stored model frame and multiplies it by the estimated fixed
#' coefficients. Sparse fixed terms (if any) are included.
#'
#' Random effects (BLUPs) are not included.
#'
#' @export
get_fixed_fit_asreml <- function(model) {
  if (!inherits(model, "asreml")) {
    stop("`model` must be an `asreml` object.")
  }

  # Get model frame
  model <- check_deisgn_exsits(model, build_mf = FALSE)
  design <- model$design

  N <- nrow(design)

  # Fixed and sparse term labels
  fixed_terms <- attr(model$formulae$fixed, "term.labels")
  sparse_terms <- attr(model$formulae$sparse, "term.labels")
  has_int <- isTRUE(attr(model$formulae$fixed, "intercept") == 1)

  # Get estimates
  term_names <- c(
    rownames(model$coefficients$fixed),
    rownames(model$coefficients$sparse)
  )
  beta <- c(
    model$coefficients$fixed,
    model$coefficients$sparse
  )

  as.numeric(design[, term_names] %*% beta)
}
