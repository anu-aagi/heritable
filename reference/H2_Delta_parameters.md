# Calculate heritability of pairwise differences using variance parameters

Compute broad-sense heritability of differences using the variance of
differences between two BLUPs/BLUEs

## Usage

``` r
h2_Delta_parameters(G_g, vd_matrix, type)

H2_Delta_parameters(vc_g, vd_matrix, type)
```

## Arguments

- vc_g:

  Numeric. Genotype variance component

- vd_matrix:

  Matrix. Variance of pairwise differences among BLUES or BLUPs

- type:

  Character. Either BLUES or BLUPS used to compute the variance of
  pairwise differences.

- G_g:

  Numeric. Genotypic variance-covariance matrix.

## Value

Matrix of pairwise heritability of differences among BLUES or BLUPs

## Details

See
[`H2_Delta()`](https://anu-aagi.github.io/heritable/reference/H2_Delta.md)
and reference for full derivation and equation for heritability Delta

## References

Schmidt, P., Hartung, J., Rath, J., & Piepho, H.-P. (2019). Estimating
Broad-Sense Heritability with Unbalanced Data from Agricultural Cultivar
Trials. Crop Science, 59(2), 525–536.
https://doi.org/10.2135/cropsci2018.06.0376

## Examples

``` r
h2_Delta_parameters(G_g = diag(0.15, 2, 2), vd_matrix = matrix(c(NA,0.2,0.2,NA),2,2), type = "BLUP")
#>           [,1]      [,2]
#> [1,]        NA 0.3333333
#> [2,] 0.3333333        NA

H2_Delta_parameters(vc_g = 0.01, vd_matrix = matrix(c(NA,0.2,0.2,NA),2,2), "BLUE")
#>            [,1]       [,2]
#> [1,]         NA 0.09090909
#> [2,] 0.09090909         NA
```
