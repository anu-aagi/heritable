# glmmTMB backend (Gaussian, identity link) ---------------------------------
#
# The high-level heritability methods are backend-agnostic: they operate on the
# output of var_comp() and on dispatched helpers (check_*, map_target_terms,
# fit_counterpart_model), all of which have glmmTMB methods elsewhere. The lme4
# method bodies therefore work unchanged for a glmmTMB model, so each glmmTMB
# method simply delegates to its lme4 counterpart. Delegation is resolved at call
# time, so file collation order does not matter.

#' @noRd
#' @export
h2_Standard.glmmTMB <- function(model, target, options = NULL, marginal = TRUE,
                                stratification = NULL, vc = NULL, ...) {
  h2_Standard.lmerMod(model, target, options = options, marginal = marginal,
                      stratification = stratification, vc = vc, ...)
}

#' @noRd
#' @export
h2_Oakey.glmmTMB <- function(model, target, options = NULL, marginal = TRUE,
                             stratification = NULL, vc = NULL, ...) {
  h2_Oakey.lmerMod(model, target, options = options, marginal = marginal,
                   stratification = stratification, vc = vc, ...)
}

#' @noRd
#' @export
h2_Cullis.glmmTMB <- function(model, target, options = NULL, marginal = TRUE,
                              stratification = NULL, vc = NULL, ...) {
  h2_Cullis.lmerMod(model, target, options = options, marginal = marginal,
                    stratification = stratification, vc = vc, ...)
}

#' @noRd
#' @export
h2_Piepho.glmmTMB <- function(model, target, options = NULL, marginal = TRUE,
                              stratification = NULL, vc = NULL, ...) {
  h2_Piepho.lmerMod(model, target, options = options, marginal = marginal,
                    stratification = stratification, vc = vc, ...)
}

#' @noRd
#' @export
h2_Delta_pairwise.glmmTMB <- function(model, target, type = c("BLUP", "BLUE"),
                                      options = NULL, marginal = TRUE,
                                      stratification = NULL, vc = NULL, ...) {
  h2_Delta_pairwise.lmerMod(model, target, type = type, options = options,
                            marginal = marginal, stratification = stratification,
                            vc = vc, ...)
}

#' @noRd
#' @export
H2_Standard.glmmTMB <- function(model, target, options = NULL, marginal = TRUE,
                                stratification = NULL, vc = NULL, ...) {
  H2_Standard.lmerMod(model, target, options = options, marginal = marginal,
                      stratification = stratification, vc = vc, ...)
}

#' @noRd
#' @export
H2_Cullis.glmmTMB <- function(model, target, options = NULL, marginal = TRUE,
                              stratification = NULL, vc = NULL, ...) {
  H2_Cullis.lmerMod(model, target, options = options, marginal = marginal,
                    stratification = stratification, vc = vc, ...)
}

#' @noRd
#' @export
H2_Oakey.glmmTMB <- function(model, target, options = NULL, marginal = TRUE,
                             stratification = NULL, vc = NULL, ...) {
  H2_Oakey.lmerMod(model, target, options = options, marginal = marginal,
                   stratification = stratification, vc = vc, ...)
}

#' @noRd
#' @export
H2_Piepho.glmmTMB <- function(model, target, options = NULL, marginal = TRUE,
                              stratification = NULL, vc = NULL, ...) {
  H2_Piepho.lmerMod(model, target, options = options, marginal = marginal,
                    stratification = stratification, vc = vc, ...)
}

#' @noRd
#' @export
H2_Delta_pairwise.glmmTMB <- function(model, target, type = c("BLUP", "BLUE"),
                                      options = NULL, marginal = TRUE,
                                      stratification = NULL, vc = NULL, ...) {
  H2_Delta_pairwise.lmerMod(model, target, type = type, options = options,
                            marginal = marginal, stratification = stratification,
                            vc = vc, ...)
}
