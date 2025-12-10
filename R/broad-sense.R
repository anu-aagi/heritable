#' @export
H2 <- function(model,
               target = NULL,
               method = c("Cullis", "Oakey", "Delta", "Piepho", "Standard"),
               options = NULL) {
  UseMethod("H2")
}

#' @importFrom stats setNames
#' @noRd
#' @export
H2.default <- function(model,
                       target = NULL,
                       method = c("Cullis", "Oakey", "Piepho", "Delta", "Standard"),
                       options
                       ) {
  method <- match.arg(method, several.ok = TRUE)

  initial_checks(model, target, options = NULL)

  # Calculate H2 for each method
  H2_values <- sapply(method, function(m) {
    switch(m,
      Cullis = H2_Cullis(model, target, options = list(check = FALSE)),
      Oakey = H2_Oakey(model, target, options = list(check = FALSE)),
      Piepho = H2_Piepho(model, target, options = list(check = FALSE)),
      Delta = H2_Delta(model, target, options = list(check = FALSE)),
      Standard = H2_Standard(model, target, options = list(check = FALSE)),
      cli::cli_abort("{.fn H2} is not implemented for method {.val {m}} of class{?es} {.code {class(model)}}")
    )
  })

  # Set names and class
  H2_values <- stats::setNames(H2_values, method)
  structure(H2_values,
    class = c("heritable", class(H2_values)),
    model = model, target = target
  )
}

#' Calculate standard broad-sense heritability
#' @export
H2_Standard <- function(model, target, options) {
  UseMethod("H2_Standard")
}

#' Calculate broad-sense heritability using Cullis' method
#' @export
H2_Cullis <- function(model, target = NULL, options) {
  UseMethod("H2_Cullis")
}

#' Calculate broad-sense heritability using Oakey's method
#' @export
H2_Oakey <- function(model, target = NULL, options) {
  UseMethod("H2_Oakey")
}

#' Calculate broad-sense heritability using Piepho's method
#' @export
H2_Piepho <- function(model, target = NULL, options) {
  UseMethod("H2_Piepho")
}

#' Calculate broad-sense heritability of differences between genotypes
#' @param type character, whether heritability is calculated using BLUEs or BLUPs
#' @param aggregate character, when taking means in the calculation, should harmonic or arithmetic mean be used?
#' @param options NULL by default, for internal checking of model object before calculations
#' @export
H2_Delta <- function(
    model,
    target = NULL,
    type = c("BLUP", "BLUE"),
    aggregate = c("arithmetic", "harmonic"),
    options) {
  UseMethod("H2_Delta")
}

#' @noRd
#' @export
H2_Delta.default <- function(model,
                             target = NULL,
                             type = c("BLUP", "BLUE"),
                             aggregate = c("arithmetic", "harmonic"),
                             options = NULL) {
  aggregate <- match.arg(aggregate)
  type <- match.arg(type)

  H2D_ij <- H2_Delta_pairwise(model, target, type = type)
  delta_values <- H2D_ij[upper.tri(H2D_ij)]

  switch(aggregate,
    "arithmetic" = mean(delta_values),
    "harmonic" = length(delta_values) / sum(1 / delta_values)
  )
}

#' Calculate average broad-sense heritability of differences for a given genotype
#' @details
#' Additional details...
#'
#' @export
H2_Delta_by_genotype <- function(model, target, type, options) {
  UseMethod("H2_Delta_by_genotype")
}

#' @export
H2_Delta_by_genotype.default <- function(model,
                                         target = NULL,
                                         type = NULL,
                                         options = NULL) {
  H2D_ij <- H2_Delta_pairwise(model, target, type, options)

  H2D_i <- as.matrix(H2D_ij) |>
    rowMeans(na.rm = TRUE) |>
    data.frame()

  H2D_i <- setNames(H2D_i, "H2D_i")

  H2D_i_list <- split(H2D_i, rownames(H2D_i))

  return(H2D_i_list)
}


#' Calculate broad-sense heritability of pariwise differences for all genotypes
#' @export
H2_Delta_pairwise <- function(model, target, type, options) {
  UseMethod("H2_Delta_pairwise")
}


