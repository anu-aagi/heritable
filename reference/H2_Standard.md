# Calculate standard heritability from model object

Compute standard heritability using the classic ratio method of
genotypic and phenotypic variance. See Falconer & Mackay (1996)

## Usage

``` r
H2_Standard(model, target, options)
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

The equation used to calculate standard heritability is:
\$\$H^2\_{Standard} = \frac{\sigma^2_g}{\sigma^2_g +
\frac{1}{n_g}\sum\_{n_g}^{i=1} \sigma^2_p / n\_{gi}}\$\$ where:

- \\n_g\\ is the number of genotypes

- \\n\_{gi}\\ is the number of replicate for a given genotype i

- \\\sigma_g\\ is the variance attributed to genotype differences

- \\\sigma_p\\ is the variance attributed to phenotypic differences

## References

Falconer, D. S., & Mackay, T. F. C. (1996). Introduction to quantitative
genetics (4th ed.). Longman.

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
H2_Standard(lettuce_lme4, target = "gen")
#> Error: object 'lettuce_lme4' not found

# asreml model (Requires license)
if (FALSE) { # \dontrun{
lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
                                 random = ~ gen,
                                 data = lettuce_subset,
                                 trace = FALSE
                                 )

H2_Standard(lettuce_asreml, target = "gen")
} # }
```
