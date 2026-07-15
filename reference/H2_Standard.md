# Calculate standard heritability from model object

Compute standard heritability using the classic ratio method of
genotypic and phenotypic variance. See Falconer & Mackay (1996)

## Usage

``` r
h2_Standard(model,
            target,
            options = NULL,
            marginal = TRUE,
            stratification = NULL,
            vc = NULL,
            ...)
H2_Standard(model,
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

The equation used to calculate standard heritability (broad-sense) is:
\$\$H^2\_{Standard} = \frac{\sigma^2_g}{\sigma^2_g +
\frac{1}{n_g}\sum\_{n_g}^{i=1} \sigma^2_p / n\_{gi}}\$\$ where:

- \\n_g\\ is the number of genotypes

- \\n\_{gi}\\ is the number of replicate for a given genotype i

- \\\sigma_g\\ is the variance attributed to genotype differences

- \\\sigma_p\\ is the variance attributed to phenotypic differences

The equation used to calculate standard heritability (narrow-sense) is:
\$\$h^2\_{Standard} = \frac{\operatorname{Var}(g_i -
g_j)}{\operatorname{Var}(y_i.. - y_j..)}\$\$ where:

- \\g_i\\ is the random effect of the \\i^{th}\\ genotype

- \\y_i..\\ is the sample average of the \\i^{th}\\ genotype

## References

Falconer, D. S., & Mackay, T. F. C. (1996). Introduction to quantitative
genetics (4th ed.). Longman.

## See also

`H2_Standard()`, h2_Standard()

## Examples

``` r
# lme4 model
lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
H2_Standard(lettuce_lme4, target = "gen")
#> [1] 0.8294971

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
