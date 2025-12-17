# Calculate heritability of differences for a given genotype from model object

Instead of computing heritability on a "entry-mean" basis, this method
calculates heritability using "entry-differences". Entry here is
referring to the genotype, line or variety of interest. See reference
for origin and interpretation of `h2/H2_Delta_by_genotype` and it's
variants

## Usage

``` r
H2_Delta_by_genotype(model, target, type = c("BLUE", "BLUP"), options)
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

Numeric

Named list, with each element containing a named numeric vector

## Details

The heritability of differences for a given genotype is given by:

\$\$H^2\_{\Delta i.} = 1 - \frac{PEV^{BLUP}\_{\overline\Delta
i.}}{2\sigma^2_g}\$\$

where:

- \\PEV^{BLUP}\_{\overline\Delta i.}\\ is the arithmetic mean of the
  prediction error variance matrix for pairwise differences among BLUPs
  (or BLUEs if `method = "BLUE"`) for a given genotype

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

[`H2_Delta()`](https://anu-aagi.github.io/heritable/reference/H2_Delta.md),
[`H2_Delta_pairwise()`](https://anu-aagi.github.io/heritable/reference/H2_Delta_pairwise.md)

## Examples

``` r
# lme4 model
lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
#> Error: bad 'data': object 'lettuce_subset' not found
H2_Delta_by_genotype(lettuce_lme4, target = "gen")
#> Error: object 'lettuce_lme4' not found

# asreml model (Requires license)
if (FALSE) { # \dontrun{
lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
                                 random = ~ gen,
                                 data = lettuce_subset,
                                 trace = FALSE
                                 )

H2_Delta_by_genotype(lettuce_asreml, target = "gen")
} # }
```
