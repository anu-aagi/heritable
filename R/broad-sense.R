#' @noRd
#' @export
H2 <- function(model,
               target,
               method = c("Cullis", "Oakey", "Piepho", "Delta", "Standard"),
               options = NULL,
               marginal = TRUE,
               stratification = NULL,
               source = list(),
               ...
               ) {
  UseMethod("H2")
}

#' @importFrom stats setNames
#' @noRd
#' @export
H2.default <- function(model,
                       target,
                       method = c("Cullis", "Oakey", "Piepho", "Delta", "Standard"),
                       options = NULL,
                       marginal = TRUE,
                       stratification = NULL,
                       source = list(),
                       ...
                       ) {
  method <- match.arg(method, several.ok = TRUE)

  initial_checks(model, target, options = options)

  # Check correct model specification.
  if(options$check %||% TRUE){
    check_model_specification(model, target, "broad_sense", source)
  }

  # Check design exists
  model <- check_deisgn_exsits(model, source = source)

  # Build variance component
  vc <- var_comp(model,
                 target = target,
                 marginal = marginal,
                 stratification = stratification,
                 ...)

  # Calculate H2 for each method
  H2_values <- sapply(method, function(m) {
    switch(m,
      Cullis = H2_Cullis(model, target, options = list(check = FALSE), vc = vc, ...),
      Oakey = H2_Oakey(model, target, options = list(check = FALSE), vc = vc, ...),
      Piepho = H2_Piepho(model, target, options = list(check = FALSE), vc = vc, ...),
      Delta = H2_Delta(model, target, options = list(check = FALSE), vc = vc, ...),
      Standard = H2_Standard(model, target, options = list(check = FALSE), vc  = vc, ...),
      cli::cli_abort("{.fn H2} is not implemented for method {.val {m}} of class{?es} {.code {class(model)}}")
    )
  })

  # Set names and class
  H2_values <- stats::setNames(H2_values, method)
  args <- list(
    model = model,
    target = target,
    method = method,
    options = options,
    marginal = marginal,
    stratification = stratification
  )
  dots <- list(...)
  if(length(dots) > 0){
    args <- c(args, dots)
  }

  structure(H2_values,
    class = c("heritable", class(H2_values)),
    type = "broad_sense",
    args = args
  )
}

#' @noRd
#' @export
H2_Standard <- function(model,
                        target,
                        options = NULL,
                        marginal = TRUE,
                        stratification = NULL,
                        vc = NULL,
                        ...) {
  UseMethod("H2_Standard")
}

#' @noRd
#' @export
H2_Cullis <- function(model,
                      target,
                      options = NULL,
                      marginal = TRUE,
                      stratification = NULL,
                      vc = NULL,
                      ...) {
  UseMethod("H2_Cullis")
}

#' @noRd
#' @export
H2_Oakey <- function(model,
                     target,
                     options = NULL,
                     marginal = TRUE,
                     stratification = NULL,
                     vc = NULL,
                     ...) {
  UseMethod("H2_Oakey")
}

#' @noRd
#' @export
H2_Piepho <- function(model,
                      target,
                      options = NULL,
                      marginal = TRUE,
                      stratification = NULL,
                      vc = NULL,
                      ...) {
  UseMethod("H2_Piepho")
}

#' @noRd
#' @export
H2_Delta <- function(
    model,
    target,
    type = c("BLUP", "BLUE"),
    options = NULL,
    marginal = TRUE,
    stratification = NULL,
    vc = NULL,
    ...) {
  UseMethod("H2_Delta")
}

#' @noRd
#' @export
H2_Delta.default <- function(model,
                             target,
                             type = c("BLUP", "BLUE"),
                             options = NULL,
                             marginal = TRUE,
                             stratification = NULL,
                             vc = NULL,
                             ...) {
  type <- match.arg(type)

  H2D_ij <- H2_Delta_pairwise(model,
                              target,
                              type = type,
                              options = options,
                              marginal =  marginal,
                              stratification = stratification,
                              vc = vc,
                              ...)

  if (is.atomic(H2D_ij) && length(H2D_ij) == 1 && is.na(H2D_ij)) {
    return(NA)
  }
  delta_values <- H2D_ij[upper.tri(H2D_ij)]

  switch(type,
    "BLUP" = mean(delta_values),
    "BLUE" = length(delta_values) / sum(1 / delta_values)
  )
}

#' @noRd
#' @export
H2_Delta_by_genotype <- function(model,
                                 target,
                                 type = c("BLUP", "BLUE"),
                                 options = NULL,
                                 marginal = TRUE,
                                 stratification = NULL,
                                 vc = NULL,
                                 ...) {
  UseMethod("H2_Delta_by_genotype")
}

#' @noRd
#' @export
H2_Delta_by_genotype.default <- function(model,
                                         target,
                                         type = c("BLUP", "BLUE"),
                                         options = NULL,
                                         marginal = TRUE,
                                         stratification = NULL,
                                         vc = NULL,
                                         ...) {
  H2D_ij <- H2_Delta_pairwise(model,
                              target,
                              type = type,
                              options = options,
                              marginal =  marginal,
                              stratification = stratification,
                              vc = vc,
                              ...)

  if (is.atomic(H2D_ij) && length(H2D_ij) == 1 && is.na(H2D_ij)) {
    return(NA)
  }

  H2D_i <- as.matrix(H2D_ij) |>
    rowMeans(na.rm = TRUE) |>
    data.frame()

  H2D_i <- setNames(H2D_i, "H2D_i")

  H2D_i_list <- split(H2D_i, rownames(H2D_i))

  return(H2D_i_list)
}

#' @noRd
#' @export
H2_Delta_pairwise <- function(model,
                              target,
                              type =  c("BLUP", "BLUE"),
                              options = NULL,
                              marginal = TRUE,
                              stratification = NULL,
                              vc = NULL,
                              ...) {
  UseMethod("H2_Delta_pairwise")
}


