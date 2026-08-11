# Changelog

## heritable 0.2.0

- [`h2()`](https://anu-aagi.github.io/heritable/reference/H2.md) and
  [`H2()`](https://anu-aagi.github.io/heritable/reference/H2.md) now
  fail with an informative error for non-Gaussian `lme4` models
  (e.g. [`glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html)) instead of
  a cryptic method-dispatch error
  [\#48](https://github.com/anu-aagi/heritable/issues/48).
- Added narrow-sense heritability estimation through
  [`h2()`](https://anu-aagi.github.io/heritable/reference/H2.md) and the
  `h2_*()` estimator functions for models with a known additive genetic
  covariance structure.
- Added a consistent interface across broad-sense
  ([`H2()`](https://anu-aagi.github.io/heritable/reference/H2.md)),
  narrow-sense
  ([`h2()`](https://anu-aagi.github.io/heritable/reference/H2.md)), and
  their Standard, Cullis, Oakey, Piepho, and Delta estimators.
- Added support for more complex `asreml` and `lme4` random-effect
  structures, including known genomic relationship matrices and
  structured residuals.
- Added richer `heritable` objects that retain the fitted model and
  calculation metadata for downstream methods.

## heritable 0.1.1

- Update Maintainer to Dr. Emi Tanaka
- Refactored `asreml` integration
  [\#37](https://github.com/anu-aagi/heritable/issues/37),
  [\#32](https://github.com/anu-aagi/heritable/issues/32)
  - Added `var_comp.asreml()` for extracting and assembling variance
    components for `asreml` model objects
- Refactored
  [`H2()`](https://anu-aagi.github.io/heritable/reference/H2.md) for
  `asreml`models so that it now takes arguments like `marginal`,
  `stratification`
  - Improved handling of `vm()` specification
- Refactored `H2_standard()`
- Compute confidence intervals for `heritable` objects
  [\#30](https://github.com/anu-aagi/heritable/issues/30)
