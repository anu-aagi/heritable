# Calculate heritability of pairwise differences using variance parameters

Compute narrow-sense or broad-sense heritability of differences using
the variance of differences between two BLUPs/BLUEs

## Usage

``` r
h2_Delta_BLUP_parameters(G_g, vd_matrix)
h2_Delta_BLUE_parameters(G_g, vd_matrix)

H2_Delta_BLUP_parameters(vc_g, vd_matrix)
H2_Delta_BLUE_parameters(vc_g, vd_matrix)
```

## Arguments

- G_g:

  Numeric. Genotypic variance-covariance matrix.

- vd_matrix:

  Matrix. Variance of pairwise differences among BLUES or BLUPs

- vc_g:

  Numeric. Genotype variance component

## Value

Matrix of pairwise heritability of differences among BLUES or BLUPs

## Details

See
[`h2_Delta()`](https://anu-aagi.github.io/heritable/reference/h2_Delta.md),
[`H2_Delta()`](https://anu-aagi.github.io/heritable/reference/h2_Delta.md)
and reference for full derivation and equation for heritability Delta

## References

Schmidt, P., Hartung, J., Rath, J., & Piepho, H.-P. (2019). Estimating
Broad-Sense Heritability with Unbalanced Data from Agricultural Cultivar
Trials. Crop Science, 59(2), 525–536.
https://doi.org/10.2135/cropsci2018.06.0376

## Examples

``` r
h2_Delta_BLUP_parameters(G_g = diag(0.15, 2, 2), vd_matrix = matrix(c(NA,0.2,0.2,NA),2,2))
#>           [,1]      [,2]
#> [1,]        NA 0.3333333
#> [2,] 0.3333333        NA

H2_Delta_BLUP_parameters(vc_g = 0.01, vd_matrix = matrix(c(NA,0.2,0.2,NA),2,2))
#>      [,1] [,2]
#> [1,]   NA   -9
#> [2,]   -9   NA
```
