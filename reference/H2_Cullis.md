# Calculate Cullis' heritability from model object

Compute "generalised heritability" for unbalanced experimental designs.
See Cullis, Smith and Coombes (2006) for derivation.

## Usage

``` r
H2_Cullis(model,
          target,
          options = NULL,
          marginal = TRUE,
          stratification = NULL,
          vc = NULL,
          ...)
h2_Cullis(model,
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
lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
H2_Cullis(lettuce_lme4, target = "gen")
#> [1] 0.8294971

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
