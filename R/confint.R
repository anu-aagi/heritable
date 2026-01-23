#' Bootstrap confidence interval for heritability
#'
#' @description
#' Computes a confidence interval for a heritability estimate using parametric
#' bootstrap of the underlying mixed model.
#'
#' @param object
#' A heritability object returned by [H2()] (broad-sense) or
#' [h2()] (narrow-sense). The object must store the fitted model
#' as an attribute.
#' @param parm a specification of which parameters are to be given confidence intervals,
#' either a vector of numbers or a vector of names.
#' If missing, all parameters are considered.
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
#' @param ... Additional arguments passed to the bootstrap routine.
#'
#' @return
#' A matrix of confidence intervals.
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
#'   my_H2 <- H2(lettuce_asreml, "gen", c("Cullis", "Standard"))
#'
#'   confint(my_H2)
#' }
#'
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

  if(is.null(parm)) parm <- seq_along(object)
  method <- names(object[parm])
  target <- attr(object, "target")
  model <- attr(object, "model") # Get the model
  h2.type <- attr(object, "type")  # Get the heritability type

  if (h2.type == "broad_sense") {
    Fun_use <- function(x) {
     H2(x, target, method, options = list(check = FALSE))
    }
  } else {
    Fun_use <- function(x) {
      h2(x, target, method, options = list(check = FALSE))
    }
  }

  if (inherits(model, "lmerMod")) {
    boot_mod <- lme4::bootMer(model,
      FUN = Fun_use, nsim = B, seed = seed,
      use.u = !resample_u, ...
    )
    ci <- stats::confint(boot_mod, level = level, type = type)
  } else {
    boot_mod <- bootstrap_asreml(model,
      FUN = Fun_use, nsim = B, seed = seed,
      use.u = !resample_u, ...
    )
    ci <- stats::confint(boot_mod, level = level, type = type)
    ci <- matrix(ci,
                 nrow = length(method),
                 dimnames = list(method, colnames(ci))
          )
  }
  attr(ci, "boot_mod") <- boot_mod
  ci
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
                             seed = NULL,
                             ...) {
  if (!inherits(model, "asreml")) {
    stop("`model` must be an `asreml` object.")
  }

  mf <- model$mf
  if (is.null(mf)) {
    stop("Model frame (`model$mf`) not found. Fit with `model.frame = TRUE` so data are stored.")
  }

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
        warning("estimateV.asreml failed: ", e$message)
        NULL
      }
    )
    if (is.null(V)) {
      warning("Random effects will not be resampled.")
      use.u <- FALSE
    }
  }

  if (use.u) {
    yhat <- model$linear.predictors
    V <- asremlPlus::estimateV.asreml(model, which.matrix = "R")
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
  refit_asreml <- function(data, model, FUN) {
    fit <- update(model, data = data)
    FUN(fit)
  }

  if (!is.null(seed)) set.seed(seed)

  boot::boot(
    data      = boot_data,
    statistic = refit_asreml,
    sim       = "parametric",
    ran.gen   = generate_data,
    R         = nsim,
    model     = model,
    mle       = L,
    FUN       = FUN,
    ...
  )
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

  design <- model$design

  # Get the design matrix
  if (is.null(design)) {
    design_default <- asreml::asreml.options()$design
    asreml::asreml.options(design = TRUE)
    model <- asreml::update.asreml(model)
    design <- model$design
    asreml::asreml.options(design = design_default)
  }

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
