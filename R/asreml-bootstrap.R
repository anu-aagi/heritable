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
bootstrap_asreml <- function(
  model,
  FUN,
  nsim = 1,
  use.u = FALSE,
  source = list(),
  seed = NULL,
  ...
) {
  if (!inherits(model, "asreml")) {
    stop("`model` must be an `asreml` object.")
  }

  # Get model frame
  design_default <- asreml::asreml.options()$design
  asreml::asreml.options(design = TRUE)
  model <- check_design_exists(model, source = source)
  mf <- model$mf

  mf <- as.data.frame(mf)
  N <- nrow(mf)

  # Get response name if not provided.
  fixed_formula <- model$formulae$fixed
  resp_idx <- attr(terms(fixed_formula), "response")
  response <- deparse(fixed_formula[[2]])

  # Fix for response transformation
  new_formula <- paste0(".y", "~", deparse(fixed_formula[[3]]))
  new_formula <- as.formula(new_formula)

  # Get the linear predictor of fixed effect.
  if (!use.u) {
    yhat <- get_fixed_fit_asreml(model, source)
    V <- tryCatch(
      asremlPlus::estimateV.asreml(model, which.matrix = "V"),
      error = function(e) {
        cli::cli_warn(
          c(
            "estimateV.asreml failed: {.code {e$message}}",
            "Try custom extraction of V due to the failure of estimateV.asreml."
          )
        )
        NULL
      }
    )

    if (is.null(V)) {
      get_V <- function(model, source) {
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

        tmp_data_call <- model$call$data
        model$call$data <- Matrix::Matrix(0, nrow = N, ncol = N)
        R <- asremlPlus::estimateV.asreml(model, which.matrix = "R") |>
          Matrix::Matrix()
        model$call$data <- tmp_data_call

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

    tmp_data_call <- model$call$data
    model$call$data <- Matrix::Matrix(0, nrow = N, ncol = N)
    V <- asremlPlus::estimateV.asreml(model, which.matrix = "R") |>
      Matrix::Matrix()
    model$call$data <- tmp_data_call
  }

  # Get V and Cholesky factor L such that L %*% z ~ N(0, V)
  # Use Matrix chol: returns upper triangular U with U' U = V
  # So L = t(U) gives L L' = V
  L <- Matrix::chol(V) |> t()

  # Bootstrap data
  boot_data <- mf
  boot_data[[".yhat"]] <- yhat
  boot_data[[".y"]] <- mf[, response]

  # generator: simulate response into correct column name
  generate_data <- function(data, mle) {
    out <- data
    N <- nrow(out)
    eps <- as.numeric(L %*% stats::rnorm(N))
    out[[".y"]] <- out[[".yhat"]] + eps
    out
  }

  # Refit wrapper for boot()
  refit_asreml <- function(data, model, FUN, source, design) {
    if (length(source) > 0) {
      list2env(source, environment())
    }
    utils::capture.output(
      fit <- asreml::update.asreml(model, data = data, fixed = new_formula)
    )
    fit[["design"]] <- design
    FUN(fit)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }
  boot <- boot::boot(
    data = boot_data,
    statistic = refit_asreml,
    sim = "parametric",
    ran.gen = generate_data,
    R = nsim,
    model = model,
    source = source,
    design = model$design,
    mle = L,
    FUN = FUN,
    ...
  )

  # Set the asreml option back
  asreml::asreml.options(design = design_default)

  boot
}
