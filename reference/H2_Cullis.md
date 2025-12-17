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
H2_Cullis(lettuce_lme4, target = "gen")
#> Error: object 'lettuce_lme4' not found

# asreml model (Requires license)
if (FALSE) { # \dontrun{
lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
                                 random = ~ gen,
                                 data = lettuce_subset,
                                 trace = FALSE
                                 )

H2_Cullis(lettuce_asreml, target = "gen")
} # }
```
