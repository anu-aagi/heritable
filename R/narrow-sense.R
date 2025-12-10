#' Calculate narrow-sense heritability
#' @name h2
#' @param model A fitted model object. Currently only supports models with class `asreml`
#' @param target The name of the random effect for which heritability is to be calculated.
#' @param method The method to use for calculating heritability. Options are "Cullis", "Oakey", "BLUE", "BLUP", "Piepho", "Reg", and "SumDiv". Default is "Cullis".
#' @param options NULL by default, for internal checking of model object before calculations
#' @export
h2 <- function(model, target, method, options) {
  UseMethod("h2")
}

#' @inheritParams h2
#' @export
h2.default <- function(
    model,
    target = NULL,
    method = c("Oakey", "Delta"),
    ...) {
  method <- match.arg(method, several.ok = TRUE)

  initial_checks(model, target, options = NULL)

  h2_values <- sapply(method, function(m) {
    switch(m,
           # Cullis = h2_Cullis(model, target, options = list(check = FALSE)),
           Oakey = h2_Oakey(model, target, options = list(check = FALSE)),
           # Piepho = h2_Piepho(model, target, options = list(check = FALSE)),
           # Delta = h2_Delta(model, target, options = list(check = FALSE)),
           # Standard = h2_Standard(model, target, options = list(check = FALSE)),
           cli::cli_abort(
             "{.fn h2} is not implemented for method {.value m} of class{?es} {.code {class(model)}}"
           )
    )
  })

  # Set names and class
  h2_values <- stats::setNames(h2_values, method)
  structure(h2_values,
            class = c("heritable", class(h2_values)),
            model = model, target = target
  )
}

#' @inheritParams h2
#' @export
h2_Oakey <- function(model, target, options) {
  UseMethod("h2_Oakey")
}

#' Calculate narrow-sense heritability of difference
#' @name h2_Delta
#' @param model A fitted model object. Currently only supports models with class `asreml`
#' @param target The name of the random effect for which heritability is to be calculated.
#' @param type Whether to use BLUEs or BLUPs for calculating heritability.
#' @param aggregate character, when taking means in the calculation, should harmonic or arithmetic mean be used?
#' @param options NULL by default, for internal checking of model object before calculations
#' @export
h2_Delta <- function(model,
                     target = NULL,
                     type = c("BLUE", "BLUP"),
                     aggregate = c("arithmetic", "harmonic"), options) {
  UseMethod("h2_Delta")
}

#' @inheritParams h2_Delta
#' @export
h2_Delta.default <- function(model,
                             target = NULL,
                             type = c("BLUP", "BLUE"),
                             aggregate = c("arithmetic", "harmonic"),
                             options = NULL) {
  aggregate <- match.arg(aggregate)
  type <- match.arg(type)

  H2D_ij <- h2_Delta_pairwise(model, target, type = type)
  delta_values <- H2D_ij[upper.tri(H2D_ij)]

  switch(aggregate,
         "arithmetic" = mean(delta_values),
         "harmonic" = length(delta_values) / sum(1 / delta_values)
  )
}

#' @inheritParams h2_Delta
#' @export
h2_Delta_pairwise <- function(model, target, type = c("BLUE", "BLUP"), options) {
  UseMethod("h2_Delta_pairwise")
}
