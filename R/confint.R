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
confint.heritable <- function(
  object,
  parm = NULL,
  level = 0.95,
  B = 100,
  random_effect = c("resample", "conditional"),
  type = c("basic", "norm", "perc"),
  return_model = TRUE,
  seed = NULL,
  ...
) {
  # basic: bias corrected percentile interval
  # norm: bias corrected normal interval
  # perc: percentile interval
  type <- match.arg(type)
  random_effect <- match.arg(random_effect)

  if (random_effect == "conditional") {
    resample_u <- FALSE
  } else {
    resample_u <- TRUE
  }

  args <- attr(object, "args")
  args[["options"]] <- list(check = FALSE)
  model <- args[["model"]]
  all_method <- args[["method"]]
  all_method <- all_method[!is.na(object)]

  # Configure parallel computing
  parallel <- list(...)[["parallel"]]
  if (!is.null(parallel) && parallel != "no") {
    if (
      inherits(model, "lmerMod") &&
        !requireNamespace("lme4", quietly = TRUE)
    ) {
      stop(
        "Package 'lme4' is required for parallel bootstrapping of lme4 models."
      )
    }

    if (
      inherits(model, "asreml") &&
        !requireNamespace("asreml", quietly = TRUE)
    ) {
      stop(
        "Package 'asreml' is required for parallel bootstrapping of ",
        "asreml models and must be installed separately."
      )
    }
  }

  if (is.null(parm)) {
    method <- all_method
  } else {
    method <- intersect(parm, all_method)
    if (length(method) == 0) {
      cli::cli_abort(
        "Invalid choice of method, should be one of {.code {all_method}}"
      )
    }
  }

  args[["method"]] <- method

  h2.type <- attr(object, "type") # Get the heritability type

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
      boot_mod <- lme4::bootMer(
        model,
        FUN = Fun_use,
        nsim = B,
        seed = seed,
        use.u = !resample_u,
        ...
      )
    )
    ci <- stats::confint(boot_mod, level = level, type = type)
  } else if (inherits(model, "asreml")) {
    suppressMessages(
      boot_mod <- bootstrap_asreml(
        model,
        FUN = Fun_use,
        nsim = B,
        seed = seed,
        use.u = !resample_u,
        source = args$source,
        ...
      )
    )
    ci <- stats::confint(boot_mod, level = level, type = type)
    ci <- matrix(
      ci,
      nrow = length(method),
      dimnames = list(method, colnames(ci))
    )
  }

  structure(ci, class = c("heritable_ci", class(ci)), boot_mod = boot_mod)
}
