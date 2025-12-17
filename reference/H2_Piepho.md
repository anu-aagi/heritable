# Calculate Piepho's heritability from model object Compute Piepho's heritability using variance differences between genotype BLUEs

Calculate Piepho's heritability from model object Compute Piepho's
heritability using variance differences between genotype BLUEs

## Usage

``` r
H2_Piepho(model, target, options)
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

## Details

The equation for Piepho's heritability is as follows:

\$\$H^2\_{Piepho} = \frac{\sigma^2_g}{\sigma^2_g +
\overline{PEV\_{BLUE_g}} / 2}\$\$

where:

- \\\overline{PEV\_{BLUE_g}}\\ is the prediction error variance matrix
  for genotype BLUEs

- \\\sigma^2_g\\ is the variance attributed to differences between
  genotype

See reference for full derivation and details.

## References

Piepho, H.-P., & Möhring, J. (2007). Computing Heritability and
Selection Response From Unbalanced Plant Breeding Trials. Genetics,
177(3), 1881–1888. https://doi.org/10.1534/genetics.107.074229

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
H2_Piepho(lettuce_lme4, target = "gen")
#> Error: object 'lettuce_lme4' not found

# asreml model (Requires license)
if (FALSE) { # \dontrun{
lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
                                 random = ~ gen,
                                 data = lettuce_subset,
                                 trace = FALSE
                                 )

H2_Piepho(lettuce_asreml, target = "gen")
} # }
```
