# Calculate Oakey's heritability from model object

Compute heritability for genotype means using the variance–covariance
matrix of the genotype BLUPs as described by Oakey et al. (2006).

## Usage

``` r
H2_Oakey(model, target, options)
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

Numeric

## Details

\$\$H^2\_{Oakey} = \frac{\sum\_{i = n_z+1}^{n_g}
\lambda_i}{\sum\_{n_g}^{\lambda_i\neq 0}}\$\$ where:

- \\n_g\\ is the number of genotypes

- \\n_z\\ is the number of zero eigenvalues

- \\\lambda_i\\ is the ith eigenvalue of the matrix \\I\_{m} -
  G^{-1}C^{gg}\\

- \\\sigma^2\\ is the variance attributed to differences between
  genotype

See pages 813 and 818 of the reference for full derivation and
explanation for Oakey's heritability

## References

Oakey, H., Verbyla, A., Pitchford, W., Cullis, B., & Kuchel, H. (2006).
Joint modeling of additive and non-additive genetic line effects in
single field trials. Theoretical and Applied Genetics, 113(5), 809–819.
https://doi.org/10.1007/s00122-006-0333-z

## Examples

``` r
# lme4 model
lettuce_phenotypes |> subset(loc == "L2")
#> # A tibble: 267 × 4
#>    loc   gen   rep       y
#>    <fct> <fct> <fct> <dbl>
#>  1 L2    G1    R1      2  
#>  2 L2    G1    R2      1.5
#>  3 L2    G1    R3      2  
#>  4 L2    G2    R1      2  
#>  5 L2    G2    R2      2  
#>  6 L2    G2    R3      1.5
#>  7 L2    G3    R1      1.5
#>  8 L2    G3    R2      2  
#>  9 L2    G3    R3      2.5
#> 10 L2    G4    R1      3.5
#> # ℹ 257 more rows
lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
#> Error: bad 'data': object 'lettuce_subset' not found
H2_Oakey(lettuce_lme4, target = "gen")
#> Error: object 'lettuce_lme4' not found

# asreml model (Requires license)
if (FALSE) { # \dontrun{
lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
                                 random = ~ gen,
                                 data = lettuce_subset,
                                 trace = FALSE
                                 )

H2_Oakey(lettuce_asreml, target = "gen")
} # }
```
