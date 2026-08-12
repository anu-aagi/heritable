# Calculate reliability from model object

Compute the reliability (\\r^2\\, also known as the coefficient of
determination) of the genotype BLUPs, as described by Schmidt et al.
(2019). `h2_Reliability()` / `H2_Reliability()` return the overall
(mean) reliability \\\bar{r}^2\\, while `h2_Reliability_by_genotype()` /
`H2_Reliability_by_genotype()` return the per-genotype values \\r^2_i\\.

## Usage

``` r
H2_Reliability(
  model,
  target,
  options = NULL,
  marginal = TRUE,
  stratification = NULL,
  vc = NULL,
  ...
)

H2_Reliability_by_genotype(
  model,
  target,
  options = NULL,
  marginal = TRUE,
  stratification = NULL,
  vc = NULL,
  ...
)

h2_Reliability(
  model,
  target,
  options = NULL,
  marginal = TRUE,
  stratification = NULL,
  vc = NULL,
  ...
)

h2_Reliability_by_genotype(
  model,
  target,
  options = NULL,
  marginal = TRUE,
  stratification = NULL,
  vc = NULL,
  ...
)
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

`h2_Reliability()` / `H2_Reliability()` return a numeric value;
`h2_Reliability_by_genotype()` / `H2_Reliability_by_genotype()` return a
named numeric vector.

## Details

The reliability of the \\i\\th genotype is \$\$r^2_i = 1 -
\frac{var(\hat{g}^{BLUP}\_i)}{var(g_i)}\$\$ where
\\var(\hat{g}^{BLUP}\_i)\\ is the \\i\\th diagonal element of the
prediction error variance matrix \\C\_{22(g)}\\ and \\var(g_i)\\ is the
\\i\\th diagonal element of the genotypic variance-covariance matrix
\\G\_{(g)}\\. As only the diagonal elements are used, \\r^2_i\\ ignores
the off-diagonal elements (covariances) of both matrices. The overall
reliability is the mean across genotypes, \$\$\bar{r}^2 =
\frac{1}{n_g}\sum\_{i=1}^{n_g} r^2_i.\$\$

## References

- Schmidt, P., Hartung, J., Bennewitz, J., & Piepho, H.-P. (2019).
  Heritability in Plant Breeding on a Genotype-Difference Basis.
  Genetics, 212(4), 991–1008.
  https://doi.org/10.1534/genetics.119.302134

- Mrode, R. A. (2014). Linear Models for the Prediction of Animal
  Breeding Values (3rd ed.). CABI.

## See also

[`H2_Reliability_parameters()`](https://anu-aagi.github.io/heritable/reference/H2_Reliability_parameters.md)

## Examples

``` r
# lme4 model
lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
lettuce_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
H2_Reliability(lettuce_lme4, target = "gen")
#> [1] 0.8201769
H2_Reliability_by_genotype(lettuce_lme4, target = "gen")
#>        G1        G2        G3        G4        G5        G6        G7        G8 
#> 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 
#>        G9       G10       G11       G12       G13       G14       G15       G16 
#> 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 
#>       G17       G18       G19       G20       G21       G22       G23       G24 
#> 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 
#>       G25       G26       G27       G28       G29       G30       G31       G32 
#> 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 
#>       G33       G34       G35       G36       G37       G39       G40       G41 
#> 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 
#>       G42       G43       G44       G45       G46       G47       G49       G50 
#> 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 
#>       G51       G52       G53       G54       G55       G56       G57       G58 
#> 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 
#>       G59       G60       G61       G62       G63       G64       G65       G66 
#> 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 
#>       G67       G68       G69       G70       G71       G72       G73       G74 
#> 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 
#>       G75       G76       G77       G78       G79       G80       G82       G83 
#> 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 
#>       G84       G85       G86       G87       G88       G89       G38       G48 
#> 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 0.8201769 
#>       G81 
#> 0.8201769 
```
