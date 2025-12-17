# Calculate average heritability of differences between genotypes from model object

Instead of computing heritability on a "entry-mean" basis, this method
calculates heritability using "entry-differences". Entry here is
referring to the genotype, line or variety of interest. See reference
for origin and interpretation of `H2_Delta` and it's variants

## Usage

``` r
H2_Delta(model,
         target,
         type = c("BLUP", "BLUE"),
         aggregate = c("arithmetic", "harmonic"),
         options
         )
```

## Arguments

- model:

  Model object of class `lmerMod/merMod` or `asreml`

- target:

  The name of the random effect for which heritability is to be
  calculated.

- type:

  character, whether heritability is calculated using BLUEs or BLUPs

- aggregate:

  character, when taking means in the calculation, should harmonic or
  arithmetic mean be used?

- options:

  NULL by default, for internal checking of model object before
  calculations

## Value

Numeric

## Details

The heritability of differences between genotypes is given by:

\$\$H^2\_{\Delta ..} = 1 - \frac{PEV^{BLUP}\_{\overline\Delta
..}}{2\sigma^2_g}\$\$

where:

- \\PEV^{BLUP}\_{\overline\Delta ..}\\ is the mean of the prediction
  error variance matrix for the pairwise differences among BLUPs (BLUEs
  if `method = "BLUE"`) across all genotypes

- \\\sigma^2\\ is the variance attributed to differences between
  genotype

See reference page 995 - 997 for full derivation of this heritability
measure and related variants

## References

Schmidt, P., Hartung, J., Rath, J., & Piepho, H.-P. (2019). Estimating
Broad-Sense Heritability with Unbalanced Data from Agricultural Cultivar
Trials. Crop Science, 59(2), 525–536.
https://doi.org/10.2135/cropsci2018.06.0376

## See also

[`H2_Delta_by_genotype()`](https://anu-aagi.github.io/heritable/reference/H2_Delta_by_genotype.md),
[`H2_Delta_pairwise()`](https://anu-aagi.github.io/heritable/reference/H2_Delta_pairwise.md)

## Examples

``` r
# lme4 model
lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
H2_Delta(lettuce_lme4, target = "gen", type = "BLUP")
#> [1] 0.8294971

# asreml model (Requires license)
if (FALSE) { # \dontrun{
lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
                                 random = ~ gen,
                                 data = lettuce_subset,
                                 trace = FALSE
                                 )

H2_Delta(lettuce_asreml, target = "gen", type = "BLUP")
} # }
```
