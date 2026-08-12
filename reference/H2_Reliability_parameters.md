# Calculate reliability using variance parameters

Compute the reliability (\\r^2\\, also known as the coefficient of
determination) for each genotype from the genotypic variance-covariance
matrix and the prediction error variance (PEV) matrix of the genotype
BLUPs.

## Usage

``` r
H2_Reliability_parameters(G_g, C22_g)
```

## Arguments

- G_g:

  Genotypic variance-covariance matrix.

- C22_g:

  Prediction error variance matrix associated with the genotype effects.

## Value

A named numeric vector of per-genotype reliabilities. Entries for
genotypes with a non-positive genotypic variance are `NA`.

## Details

The reliability of the \\i\\th genotype is

\$\$r^2_i = 1 - \frac{var(\hat{g}^{BLUP}\_i)}{var(g_i)}\$\$

where:

- \\var(\hat{g}^{BLUP}\_i)\\ is the \\i\\th diagonal element of the PEV
  matrix \\C\_{22(g)}\\

- \\var(g_i)\\ is the \\i\\th diagonal element of the genotypic
  variance-covariance matrix \\G\_{(g)}\\

As only the diagonal elements are used, \\r^2_i\\ does not account for
the off-diagonal elements (i.e. covariances) of either matrix. The
overall reliability is obtained as the mean across genotypes,
\\\bar{r}^2 = \frac{1}{n_g}\sum\_{i=1}^{n_g} r^2_i\\.

Genotypes whose genotypic variance (the corresponding diagonal element
of \\G\_{(g)}\\) is not strictly positive yield an undefined reliability
and are returned as `NA` (with a warning) rather than `Inf`/`NaN`.

## References

- Schmidt, P., Hartung, J., Bennewitz, J., & Piepho, H.-P. (2019).
  Heritability in Plant Breeding on a Genotype-Difference Basis.
  Genetics, 212(4), 991–1008.
  https://doi.org/10.1534/genetics.119.302134

- Mrode, R. A. (2014). Linear Models for the Prediction of Animal
  Breeding Values (3rd ed.). CABI.

## Examples

``` r
G_g <- diag(0.15, 3, 3)
C22_g <- matrix(
  c(
    0.08, 0.01, 0.00,
    0.01, 0.07, 0.01,
    0.00, 0.01, 0.09
  ),
  nrow = 3, byrow = TRUE
)
H2_Reliability_parameters(G_g, C22_g)
#> [1] 0.4666667 0.5333333 0.4000000
```
