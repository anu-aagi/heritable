# Calculate pairwise heritability of differences between genotypes from model object

Instead of computing heritability on a "entry-mean" basis, this method
calculates heritability using "entry-differences". Entry here is
referring to the genotype, line or variety of interest. See reference
for origin and interpretation of `h2/H2_Delta_pairwise` and it's
variants

## Usage

``` r
H2_Delta_pairwise(model, target, type = c("BLUE", "BLUP"), options)
```

## Arguments

- model:

  Model object of class `lmerMod/merMod` or `asreml`

- target:

  The name of the random effect for which heritability is to be
  calculated.

- type:

  character, whether heritability is calculated using BLUEs or BLUPs

- options:

  NULL by default, for internal checking of model object before
  calculations

## Value

A `dspMatrix`

## References

Schmidt, P., Hartung, J., Rath, J., & Piepho, H.-P. (2019). Estimating
Broad-Sense Heritability with Unbalanced Data from Agricultural Cultivar
Trials. Crop Science, 59(2), 525–536.
https://doi.org/10.2135/cropsci2018.06.0376

## See also

[`H2_Delta_by_genotype()`](https://anu-aagi.github.io/heritable/reference/H2_Delta_by_genotype.md),
[`H2_Delta()`](https://anu-aagi.github.io/heritable/reference/H2_Delta.md)

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
H2_Delta_pairwise(lettuce_lme4, target = "gen")
#> Error: object 'lettuce_lme4' not found

# asreml model (Requires license)
if (FALSE) { # \dontrun{
lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
                                 random = ~ gen,
                                 data = lettuce_subset,
                                 trace = FALSE
                                 )

H2_Delta_pairwise(lettuce_asreml, target = "gen")
} # }
```
