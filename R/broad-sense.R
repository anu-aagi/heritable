#' @noRd
#' @export
H2 <- function(model,
               target,
               method = c("Cullis", "Oakey", "Piepho", "Delta", "Standard"),
               options = NULL,
               marginal = TRUE,
               stratification = NULL,
               source = list(),
               vc = NULL,
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
                       vc = NULL,
                       ...
                       ) {
  method <- match.arg(
    method,
    c("Cullis", "Oakey", "Piepho", "Delta", "Standard", "Reliability"),
    several.ok = TRUE
  )

  initial_checks(model, target, options = options)

  # Check correct model specification.
  if(options$check %||% TRUE){
    check_model_specification(model, target, "broad_sense", source)
  }

  # Check design exists
  model <- check_deisgn_exsits(model, source = source)

  # Build variance component
  if(is.null(vc)){
    vc <- var_comp(model,
                   target = target,
                   marginal = marginal,
                   stratification = stratification,
                   ...)
  }

  # Calculate H2 for each method
  H2_values <- sapply(method, function(m) {
    switch(m,
      Cullis = H2_Cullis(model, target, options = list(check = FALSE), vc = vc, ...),
      Oakey = H2_Oakey(model, target, options = list(check = FALSE), vc = vc, ...),
      Piepho = H2_Piepho(model, target, options = list(check = FALSE), vc = vc, ...),
      Delta = H2_Delta(model, target, options = list(check = FALSE), vc = vc, ...),
      Standard = H2_Standard(model, target, options = list(check = FALSE), vc  = vc, ...),
      Reliability = H2_Reliability(model, target, options = list(check = FALSE), vc = vc, ...),
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

#' @rdname H2_Reliability
#' @export
H2_Reliability <- function(model,
                           target,
                           options = NULL,
                           marginal = TRUE,
                           stratification = NULL,
                           vc = NULL,
                           ...) {
  UseMethod("H2_Reliability")
}

#' @noRd
#' @export
H2_Reliability.default <- function(model,
                                   target,
                                   options = NULL,
                                   marginal = TRUE,
                                   stratification = NULL,
                                   vc = NULL,
                                   ...) {
  r2 <- H2_Reliability_by_genotype(model,
                                   target,
                                   options = options,
                                   marginal = marginal,
                                   stratification = stratification,
                                   vc = vc,
                                   ...)

  if (is.atomic(r2) && length(r2) == 1 && is.na(r2)) {
    return(NA)
  }

  mean(r2)
}

#' @rdname H2_Reliability
#' @export
H2_Reliability_by_genotype <- function(model,
                                       target,
                                       options = NULL,
                                       marginal = TRUE,
                                       stratification = NULL,
                                       vc = NULL,
                                       ...) {
  UseMethod("H2_Reliability_by_genotype")
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
  delta_g <- attr(H2D_ij, "delta_g")
  delta_pev <- mean(attr(H2D_ij, "delta_pev")[upper.tri(H2D_ij)])

  switch(type,
    "BLUP" = 1 - delta_pev / delta_g,
    "BLUE" = delta_g / (delta_g + delta_pev)
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

  delta_g <- attr(H2D_ij, "delta_g")
  delta_pev <- attr(H2D_ij, "delta_pev")
  delta_pev <- Matrix::rowSums(delta_pev)/(ncol(delta_pev)-1)

  if(type == "BLUP"){
    return(1 - delta_pev / delta_g)
  } else {
    return(delta_g / (delta_g + delta_pev))
  }

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


