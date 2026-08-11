

<!-- README.md is generated from README.qmd. Please edit that file -->

# heritable

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://img.shields.io/cran/v/heritable.svg)](https://CRAN.R-project.org/package=heritable)
[![Codecov test
coverage](https://codecov.io/gh/anu-aagi/heritable/graph/badge.svg)](https://app.codecov.io/gh/anu-aagi/heritable)
<!-- badges: end -->

`heritable` provides flexible and interpretable heritability estimation
from linear mixed models in R. Given genotypes and their observed
traits, it can help quantify the proportion of phenotypic variation that
is attributable to genetic variation. This information is useful for
understanding the genetic architecture of traits, predicting breeding
values, and designing breeding programs.

The package brings commonly used estimators into a unified interface and
supports:

- broad-sense (`H2()`) and narrow-sense (`h2()`) heritability;
- Cullis, Oakey, Piepho, Standard, and Delta estimators, including BLUP-
  and BLUE-based Delta formulations;
- models fitted with [`lme4`](https://cran.r-project.org/package=lme4)
  and [`asreml`](https://vsni.co.uk/software/asreml-r/);
- marginal (environment-averaged) and stratified (environment-specific)
  heritability for genotype-by-environment models; and
- confidence intervals based on parametric bootstrap.

## Installation

The released version is available from CRAN:

``` r
install.packages("heritable")
```

The package is under active development. You can install the development
version from [GitHub](https://github.com/anu-aagi/heritable) with:

``` r
# install.packages("pak")
pak::pak("anu-aagi/heritable")
```

## A simple demo

This basic example calculates broad-sense heritability for a single
trial for a particular trait. The same `H2()` interface works with
`lme4` and `asreml` model objects. Below we fit a simple linear mixed
model using `lme4` in which genotype effects are treated as random and
the `rep` term is included as a fixed effect.

``` r
library(heritable)

lettuce_subset <- subset(lettuce_phenotypes, loc == "L2")
fit <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
```

By the definition of heritability, the genetic effects must be treated
as random. The `target` argument specifies the term encoding genetic
effects in the fitted model. In this case, the random effect term is
`"gen"`.

Uppercase `H2()` calculates broad-sense heritability. Use `target` to
specify the term encoding genetic effects, `method` to select one or
more estimators:

``` r
H2(fit, target = "gen", method = c("Cullis", "Standard"))
#>    Cullis  Standard 
#> 0.8294971 0.8294971
```

`heritable` also supports multi-environment trial (MET) models by
considering either the main, marginal or stratified effects. In the
example below, we use `asreml` to fit a genotype-by-environment (G×E)
model in which genotype, environment, and G×E interaction effects are
all treated as random. Currently the model must contain the main
genotype effects and do not extend to factor analytic structures.

``` r
fit_gxe <- asreml::asreml(
  fixed = y ~ rep,
  random = ~ gen + loc + gen:loc,
  trace = FALSE,
  data = lettuce_phenotypes
)
```

For G×E models, `H2()` supports multiple definitions of heritability.

The classical definition treats the G×E interaction as a nuisance source
of variation. To calculate this form of heritability, set
`marginal = FALSE`:

``` r
H2(fit_gxe, target = "gen", marginal = FALSE)
#>    Cullis     Oakey    Piepho     Delta  Standard 
#> 0.7818483 0.7818483 0.7806164 0.7818483 0.7750432
```

Alternatively, setting `marginal = TRUE` treats the G×E interaction as
part of the genetic variation expressed across environments. In this
case, `H2()` estimates heritability averaged over the target
environments.

``` r
H2(fit_gxe, target = "gen", marginal = TRUE)
#>    Cullis     Oakey    Piepho     Delta  Standard 
#> 0.8918403 0.8816483 0.8895923 0.8918403 0.8895923
```

You can also estimate stratified heritability for a specific environment
by supplying a one-row data frame to the `stratification` argument:

``` r
H2(fit_gxe, target = "gen", stratification = data.frame(loc = "L1"))
#>    Cullis     Oakey    Piepho     Delta  Standard 
#> 0.7870171 0.7780585 0.7250283 0.7870171 0.7250283
```

In this example, heritability is defined and estimated specifically for
environment `"L1"`.

Lowercase `h2()` estimates **narrow-sense heritability** for models that
include a known additive genetic covariance structure, such as a genomic
relationship matrix (GRM). In `asreml`, this is specified using `vm()`.
Below, we fit a G×E model in which the additive genetic effects are
modelled using the GRM.

``` r
fit_grm <- asreml::asreml(
  fixed = y ~ rep,
  random = ~ vm(gen, source = lettuce_GRM) * idv(loc),
  data = lettuce_phenotypes,
  trace = FALSE
)
```

The interface of `h2()` is nearly identical to that of `H2()`. The only
additional argument is `source`, which specifies the genetic covariance
structure used when fitting the model. This information is required
because the fitted model does not retain the covariance matrix itself.

``` r
h2(fit_grm, target = "gen", source = list(lettuce_GRM = lettuce_GRM))
#>    Cullis     Oakey    Piepho     Delta  Standard 
#> 0.9233324 0.6867032 0.8980665 0.9228154 0.8974431
```

Like `H2()`, `h2()` supports both `marginal` and `stratification`
arguments for G×E models.

The returned estimate retains the fitted model, allowing parametric
bootstrap confidence intervals to be obtained using the standard
`confint()` generic.

``` r
estimate <- H2(fit, target = "gen", method = c("Cullis", "Standard"))
confint(estimate, B = 100, seed = 2026)
#>              2.5 %    97.5 %
#> Cullis   0.7792096 0.9395049
#> Standard 0.7792096 0.9395049
```

See the [Get started
guide](https://anu-aagi.github.io/heritable/articles/heritable.html) for
further examples and the
[presentation](https://anu-aagi.github.io/heritable/slides/heritable-slide-for-aagiverse-meeting.html)
for the motivation and statistical framework.

## Too slow to run?

For large or complex models, heritability estimation can be
computationally intensive. Runtime can be reduced by approximating the
required matrix pseudoinverses and computing error-variance matrices
from the mixed-model equations:

``` r
options(exact_psd_inv = FALSE)
H2(fit, target = "gen", solver = "LMM")
#>    Cullis     Oakey    Piepho     Delta  Standard 
#> 0.8294971 0.8294971 0.8294876 0.8294971 0.8294971
```

This approach can substantially reduce runtime, although it may be less
numerically stable than the default calculation. The speed improvement
is generally more apparent for models with many observations, random
effects, or environments.

## Support our work!

``` r
citation("heritable")
#> To cite package 'heritable' in publications use:
#> 
#>   Tanaka E, Deng Y (2026). _heritable: Heritability Estimation from
#>   Linear Mixed Models_. R package version 0.2.0,
#>   <https://anu-aagi.github.io/heritable/>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {heritable: Heritability Estimation from Linear Mixed Models},
#>     author = {Emi Tanaka and Yidi Deng},
#>     year = {2026},
#>     note = {R package version 0.2.0},
#>     url = {https://anu-aagi.github.io/heritable/},
#>   }
```

## Found a bug? Want to contribute?

We welcome feedback and contributions of all shapes and sizes! Take a
look at our existing
[backlog](https://github.com/anu-aagi/heritable/issues), add to a
relevant discussion, or lodge a new issue.
