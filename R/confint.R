#' Bootstrap confidence interval for heritability
#'
#' @description
#' Computes a confidence interval for a heritability estimate using parametric
#' bootstrap of the underlying mixed model.
#'
#' @param heritable
#' A heritability object returned by \code{heritable::H2()} (broad-sense) or
#' \code{heritable::h2()} (narrow-sense). The object must store the fitted model
#' as an attribute.
#' @param B Integer. Number of bootstrap replicates.
#' @param seed Optional random seed.
#' @param random_effect Character. Indicate whether to resample random
#' effect; one of \code{"resample"} or \code{"fix"}.
#' @param level Confidence level.
#' @param type Character. Bootstrap interval type; one of \code{"basic"},
#' \code{"norm"}, or \code{"perc"}.
#' @param return_model Logical. Whether to return to the \code{boot} object.
#' @param ... Additional arguments passed to the bootstrap routine.
#'
#' @return
#' A confidence interval object as returned by \code{confint()}.
#'
#' @examples
#' lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
#' lettuce_asreml <- asreml(
#'   fixed = y ~ rep * pseudo_var1,
#'   random = ~gen,
#'   sparse = ~pseudo_var2,
#'   data = lettuce_subset,
#'   trace = FALSE
#' )
#'
#' my_H2 <- H2(lettuce_asreml, "gen", c("Cullis", "Standard"))
#'
#' confint(my_H2)
#'
#' @export
confint.heritable <- function(heritable,
                              B = 100,
                              seed = NULL,
                              random_effect = c("resample", "fix"),
                              level = 0.95,
                              type = c("basic", "norm", "perc"),
                              return_model = TRUE,
                              ...) {
  # basic: bias corrected percentile interval
  # norm: bias corrected normal interval
  # perc: percentile interval
  type <- match.arg(type)
  random_effect <- match.arg(random_effect)

  if(random_effect == "fix"){
    resample_u <- FALSE
  } else {
    resample_u <- TRUE
  }

  method <- names(heritable)
  target <- attr(heritable, "target")
  model <- attr(heritable, "model") # Get the model
  h2.type <- attr(heritable, "type")  # Get the heritability type

  if (h2.type == "broad_sense") {
    Fun_use <- function(x) {
      heritable::H2(x, target, method, options = list(check = FALSE))
    }
  } else {
    Fun_use <- function(x) {
      heritable::h2(x, target, method, options = list(check = FALSE))
    }
  }

  if (inherits(model, "lmerMod")) {
    boot_mod <- lme4::bootMer(model,
      FUN = Fun_use, nsim = B, seed = seed,
      use.u = !resample_u, ...
    )
    ci <- confint(boot_mod, level = level, type = type)
  } else {
    boot_mod <- bootstrap_asreml(model,
      FUN = Fun_use, nsim = B, seed = seed,
      use.u = !resample_u, ...
    )
    ci <- confint(boot_mod, level = level, type = type)
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
#' lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
#' lettuce_asreml <- asreml(
#'   fixed = y ~ rep * pseudo_var1,
#'   random = ~gen,
#'   sparse = ~pseudo_var2,
#'   data = lettuce_subset,
#'   trace = FALSE
#' )
#'
#' b <- bootstrap_asreml(
#'   lettuce_asreml,
#'   R = 200,
#'   statistic = function(fit) coef(fit)$fixed["(Intercept)", "effect"],
#'   seed = 1
#' )
#' boot::boot.ci(b, type = "perc")
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
  L <- Matrix::chol(V) %>% t()

  # Bootstrap data
  boot_data <- mf
  boot_data[[".yhat"]] <- yhat

  # Call modified to use the boot "data" argument
  call_asreml <- deparse1(model$call) %>%
    paste0("asreml::", .) %>%
    str2lang()
  call_asreml <- rlang::call_modify(call_asreml, data = quote(data))

  # generator: simulate response into correct column name
  generate_data <- function(data, mle) {
    out <- data
    N <- nrow(out)
    eps <- as.numeric(L %*% rnorm(N))
    out[[response]] <- out[[".yhat"]] + eps
    out
  }

  # Refit wrapper for boot()
  refit_asreml <- function(data, model, FUN) {
    fit <- asreml::update.asreml(model, data = data)
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

  # Geth the design matrix
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
