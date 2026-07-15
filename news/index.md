# Changelog

## heritable 0.2.0

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
