# Calculate heritability of differences for a given genotype from model object

Instead of computing heritability on a "entry-mean" basis, this method
calculates heritability using "entry-differences". Entry here is
referring to the genotype, line or variety of interest. See reference
for origin and interpretation of `h2/H2_Delta_by_genotype` and it's
variants

## Usage

``` r
h2_Delta_by_genotype(model,
                     target,
                     type = c("BLUP", "BLUE"),
                     options = NULL,
                     marginal = TRUE,
                     stratification = NULL,
                     vc = NULL,
                     ...)
H2_Delta_by_genotype(model,
                     target,
                     type = c("BLUP", "BLUE"),
                     options = NULL,
                     marginal = TRUE,
                     stratification = NULL,
                     vc = NULL,
                     ...)
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

- marginal:

  Logical; if `TRUE`, construct marginal (strata-averaged) mappings so
  that each genotype receives a single averaged effect per term. If
  `FALSE`, mappings will only consider the main genotype effect and
  ignore the iteracting terms.

- stratification:

  A one-row data frame defining the stratum in which genotype effects
  should be evaluated. The columns must correspond to model terms that
  interact with `target`.

- vc:

  A list of precomputed variance components. Should be in the same
  structure as the output of
  [`var_comp()`](https://anu-aagi.github.io/heritable/reference/var_comp.md)

- ...:

  Additional arguments that specify heritability calculation when
  interactions with genotype effects are modelled

## Value

Numeric

Named list, with each element containing a named numeric vector

## Details

The broad-sense heritability of differences between genotypes is given
by:

\$\$H^2\_{\Delta ..} = 1 - \frac{PEV^{BLUP}\_{\overline\Delta
..}}{2\sigma^2_g}\$\$

where:

- \\PEV^{BLUP}\_{\overline\Delta ..}\\ is the mean of the prediction
  error variance matrix for the pairwise differences among BLUPs (BLUEs
  if `method = "BLUE"`) across all genotypes

- \\\sigma^2\\ is the variance attributed to differences between
  genotype

The narrow-sense heritability of differences between genotypes is given
by:

\$\$h^2\_{\Delta ij} = 1 - \frac{PEV^{BLUP}\_{\overline\Delta
ij}}{\operatorname{Var}(g_i - g_j)}\$\$

where:

- \\g_i\\ is the random effect of the \\i^{th}\\ genotype

See reference page 995 - 997 for full derivation of this heritability
measure and related variants

## References

Schmidt, P., Hartung, J., Rath, J., & Piepho, H.-P. (2019). Estimating
Broad-Sense Heritability with Unbalanced Data from Agricultural Cultivar
Trials. Crop Science, 59(2), 525–536.
https://doi.org/10.2135/cropsci2018.06.0376

## See also

[`h2_Delta()`](https://anu-aagi.github.io/heritable/reference/H2_Delta.md),
[`H2_Delta()`](https://anu-aagi.github.io/heritable/reference/H2_Delta.md),
[`h2_Delta_pairwise()`](https://anu-aagi.github.io/heritable/reference/H2_Delta_pairwise.md),
[`H2_Delta_pairwise()`](https://anu-aagi.github.io/heritable/reference/H2_Delta_pairwise.md)

## Examples

``` r
# lme4 model
lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
H2_Delta_by_genotype(lettuce_lme4, target = "gen", type = "BLUP")
#>        G1        G2        G3        G4        G5        G6        G7        G8 
#> 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 
#>        G9       G10       G11       G12       G13       G14       G15       G16 
#> 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 
#>       G17       G18       G19       G20       G21       G22       G23       G24 
#> 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 
#>       G25       G26       G27       G28       G29       G30       G31       G32 
#> 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 
#>       G33       G34       G35       G36       G37       G39       G40       G41 
#> 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 
#>       G42       G43       G44       G45       G46       G47       G49       G50 
#> 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 
#>       G51       G52       G53       G54       G55       G56       G57       G58 
#> 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 
#>       G59       G60       G61       G62       G63       G64       G65       G66 
#> 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 
#>       G67       G68       G69       G70       G71       G72       G73       G74 
#> 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 
#>       G75       G76       G77       G78       G79       G80       G82       G83 
#> 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 
#>       G84       G85       G86       G87       G88       G89       G38       G48 
#> 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971 
#>       G81 
#> 0.8294971 

# asreml model (Requires license)
if (FALSE) { # \dontrun{
lettuce_asreml <- asreml::asreml(fixed = y ~ rep,
                                 random = ~ gen,
                                 data = lettuce_subset,
                                 trace = FALSE
                                 )

H2_Delta_by_genotype(lettuce_asreml, target = "gen", type = "BLUP")
} # }
```
