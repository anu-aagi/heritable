# heritable

`heritable` provides flexible and interpretable heritability estimation
from linear mixed models in R. Given genotypes and their observed
traits, it helps breeders assess whether a trial contains enough genetic
signal for selection.

The package brings commonly used estimators into a unified interface and
supports:

- broad-sense
  ([`H2()`](https://anu-aagi.github.io/heritable/reference/H2.md)) and
  narrow-sense
  ([`h2()`](https://anu-aagi.github.io/heritable/reference/H2.md))
  heritability;
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

This basic example calculates broad-sense heritability for a
single-environment trial. The same
[`H2()`](https://anu-aagi.github.io/heritable/reference/H2.md) interface
works with `lme4` and `asreml` model objects.

``` r

library(heritable)

lettuce_subset <- subset(lettuce_phenotypes, loc == "L2")
fit_lme4 <- lme4::lmer(y ~ rep + (1 | gen), data = lettuce_subset)
```

Uppercase [`H2()`](https://anu-aagi.github.io/heritable/reference/H2.md)
calculates broad-sense heritability. Use `method` to select one or more
estimators:

``` r

H2(
  fit_lme4,
  target = "gen",
  method = c("Cullis", "Oakey", "Piepho", "Delta", "Standard")
)
#>    Cullis     Oakey    Piepho     Delta  Standard 
#> 0.8294971 0.8294971 0.8294971 0.8294971 0.8294971
```

Lowercase [`h2()`](https://anu-aagi.github.io/heritable/reference/H2.md)
calculates narrow-sense heritability when the fitted model includes a
known additive genetic covariance structure, such as a genomic
relationship matrix. With `asreml`, this can be specified using `vm()`:

``` r

fit_genomic <- asreml::asreml(
  fixed = y ~ rep,
  random = ~ asreml::vm(gen, source = lettuce_GRM) * asreml::idv(loc),
  data = lettuce_phenotypes,
  trace = FALSE
)

h2(
  fit_genomic,
  target = "gen",
  marginal = TRUE,
  source = list(lettuce_GRM = lettuce_GRM)
)
```

For genotype-by-environment models, `marginal = TRUE` estimates
heritability averaged across environments. Alternatively, supply a
one-row data frame to `stratification` to estimate heritability within a
particular environment:

``` r

H2(fit_gxe, target = "gen", marginal = TRUE)
H2(fit_gxe, target = "gen", stratification = data.frame(loc = "L1"))
```

Estimates returned by
[`H2()`](https://anu-aagi.github.io/heritable/reference/H2.md) and
[`h2()`](https://anu-aagi.github.io/heritable/reference/H2.md) retain
the fitted model, allowing parametric-bootstrap confidence intervals
through the standard [`confint()`](https://rdrr.io/r/stats/confint.html)
generic:

``` r

estimate <- H2(fit_lme4, target = "gen", method = c("Cullis", "Standard"))
confint(estimate, B = 500, seed = 2026)
```

See the [Get started
guide](https://anu-aagi.github.io/heritable/articles/heritable.html) for
further examples and the
[presentation](https://anu-aagi.github.io/heritable/slides/heritable-slide-for-aagiverse-meeting.html)
for the motivation and statistical framework.

## Support our work!

``` r

citation("heritable")
#> To cite package 'heritable' in publications use:
#> 
#>   Kar F, Deng Y, Li W, Tanaka E (2025). _heritable: R package for heritability
#>   calculations for plant breeding trials_. R package version 0.2.0,
#>   <https://github.com/anu-aagi/heritable>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {heritable: R package for heritability calculations for plant breeding trials},
#>     author = {Fonti Kar and Yidi Deng and Weihao (Patrick) Li and Emi Tanaka},
#>     year = {2026},
#>     note = {R package version 0.2.0},
#>     url = {https://github.com/anu-aagi/heritable},
#>   }
```

## Found a bug? Want to contribute?

We welcome feedback and contributions of all shapes and sizes! Take a
look at our existing
[backlog](https://github.com/anu-aagi/heritable/issues), add to a
relevant discussion, or lodge a new issue.
