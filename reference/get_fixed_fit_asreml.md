# Fixed-effects-only fitted values from an asreml model

Returns the fitted values based on fixed effects only (\\\hat{y} =
X\hat{\beta}\\), excluding all random effects.

## Usage

``` r
get_fixed_fit_asreml(model, source = list())
```

## Arguments

- model:

  An object of class `"asreml"`, fitted with `model.frame = TRUE`.

- source:

  The known genomic relationship matrix (GRM) used in `model` fitted
  using `asreml::vm()`, provided as a named list. When not provided (an
  empty list by default), the GRM variable used for `vm` calling will be
  searched in the global environment.

## Value

A numeric vector of length \\N\\, giving the fixed-effects-only fitted
value for each observation.

## Details

This function reconstructs the fixed-effect design matrix from the
stored model frame and multiplies it by the estimated fixed
coefficients. Sparse fixed terms (if any) are included.

Random effects (BLUPs) are not included.
