# Calculate Piepho's heritability using variance parameters

Compute Piepho's heritability using the variance of differences between
two BLUES.

## Usage

``` r
H2_Piepho_parameters(vc_g, vd_BLUE_avg)
```

## Arguments

- vc_g:

  Numeric. Genotype variance component

- vd_BLUE_avg:

  Numeric. Mean variance of pairwise differences among BLUES

## Value

Single numeric value

## Details

The equation for Piepho's heritability is as follows:

\$\$H^2\_{Piepho} = \frac{\sigma^2_g}{\sigma^2_g +
\overline{PEV\_{BLUE_g}} / 2}\$\$

where:

- \\\overline{PEV\_{BLUE_g}}\\ is the prediction error variance matrix
  for genotype BLUEs

- \\\sigma^2_g\\ is the variance attributed to differences between
  genotype

## References

Piepho, H.-P., & Möhring, J. (2007). Computing Heritability and
Selection Response From Unbalanced Plant Breeding Trials. Genetics,
177(3), 1881–1888. https://doi.org/10.1534/genetics.107.074229

## Examples

``` r
H2_Piepho_parameters(vc_g = 0.25, vd_BLUE_avg = 0.68)
#> [1] 0.4237288
```
