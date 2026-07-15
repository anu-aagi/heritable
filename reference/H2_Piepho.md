# Calculate Piepho's heritability from model object Compute Piepho's heritability using variance differences between genotype BLUEs

Calculate Piepho's heritability from model object Compute Piepho's
heritability using variance differences between genotype BLUEs

## Usage

``` r
H2_Piepho(model,
          target,
          options = NULL,
          marginal = TRUE,
          stratification = NULL,
          vc = NULL,
          ...)
h2_Piepho(model,
          target,
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
lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
H2_Piepho(lettuce_lme4, target = "gen")
#> [1] 0.8294971

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
