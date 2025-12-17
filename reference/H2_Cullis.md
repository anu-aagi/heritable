# Calculate Cullis' heritability from model object

Compute "generalised heritability" for unbalanced experimental designs.
See Cullis, Smith and Coombes (2006) for derivation.

## Usage

``` r
H2_Cullis(model, target, options)
```

## Arguments

- model:

  Model object of class `lmerMod/merMod` or `asreml`

- target:

  The name of the random effect for which heritability is to be
  calculated.

- options:

  NULL by default, for internal checking of model object before
  calculations

## Value

Numeric value

## Details

The equation for Cullis heritability is as follow

\$\$H^2\_{Cullis} = 1 - \frac{PEV^{BLUP}\_{\overline\Delta
ij}}{2\sigma^2_g}\$\$

where:

- \\PEV\\ is the prediction error variance matrix of the pairwise
  differences among BLUPS

- \\\sigma^2\\ is the variance attributed to differences between
  genotype

## References

Cullis, B. R., Smith, A. B., & Coombes, N. E. (2006). On the design of
early generation variety trials with correlated data. Journal of
Agricultural, Biological, and Environmental Statistics, 11(4), 381–393.
https://doi.org/10.1198/108571106X154443

## Examples
