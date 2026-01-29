

<!-- README.md is generated from README.qmd. Please edit that file -->

# heritable

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test
coverage](https://codecov.io/gh/anu-aagi/heritable/graph/badge.svg)](https://app.codecov.io/gh/anu-aagi/heritable)
[![CRAN status](https://www.r-pkg.org/badges/version/heritable)](https://CRAN.R-project.org/package=heritable)
<!-- badges: end -->

`heritable` is to be the one-stop shop for heritability calculations in
R. Our goal is to implement existing methods for heritability to aid
reproducibility and reporting of it’s calculations.`heritable` works
with model outputs from `asreml` and `lme4` and support the calculations
of broad and narrow sense heritability of a variety of method for single
environment breeding trials.

## Installation

Note that this package is under active development. You can install the
development version of heritable from
[GitHub](https://github.com/anu-aagi/heritable) with:

``` r
# install.packages("pak")
pak::pak("anu-aagi/heritable")
```

## A simple demo

This is a basic example which shows you how to calculate broad-sense
heritability for a single environment trial using `asreml` and `lme4`.

``` r
library(heritable)  

fit_asreml <- asreml::asreml(yield ~ rep,
    random = ~ gen + rep:block,
    data = agridat::john.alpha,
    trace = FALSE
  )

fit_lme4 <- lme4::lmer(yield ~ rep + (1|gen) + (1|rep:block),
                       data = agridat::john.alpha)
```

The `H2()` function refers to broad-sense heritability and by default,
it will compute all the available heritability methods for you.

``` r
H2(fit_asreml, target = "gen")
#>    Cullis     Oakey    Piepho     Delta  Standard 
#> 0.8090841 0.8090841 0.8029760 0.8090841 0.8400648
H2(fit_lme4, target = "gen")
#>    Cullis     Oakey    Piepho     Delta  Standard 
#> 0.8091339 0.8091339 0.7966376 0.8091339 0.8400679
```

Alternatively, with the help of tidyverse functions, you can return all
the measures for different models as a tibble.

``` r
library(tidyverse)
tibble(model = list(fit_lme4, fit_asreml)) |> 
  mutate(H2 = map(model, ~H2(.x, target = "gen"))) |> 
  unnest_wider(H2)
#> # A tibble: 2 × 6
#>   model     Cullis Oakey Piepho Delta Standard
#>   <list>     <dbl> <dbl>  <dbl> <dbl>    <dbl>
#> 1 <lmerMod>  0.809 0.809  0.797 0.809    0.840
#> 2 <asreml>   0.809 0.809  0.803 0.809    0.840
```

## Support our work!

``` r
citation("heritable")
#> To cite package 'heritable' in publications use:
#> 
#>   Kar F, Tanaka E (2025). _heritable: R package for heritability
#>   calculations for plant breeding trials_. R package version 0.0.9000,
#>   <https://github.com/anu-aagi/heritable>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {heritable: R package for heritability calculations for plant breeding trials},
#>     author = {Fonti Kar and Emi Tanaka},
#>     year = {2025},
#>     note = {R package version 0.0.9000},
#>     url = {https://github.com/anu-aagi/heritable},
#>   }
```

## Found a bug? Want to contribute?

We welcome feedback and contributions of all shapes and sizes! Take a
look at our existing
[backlog](https://github.com/anu-aagi/heritable/issues) and add to a
relevant discussion or lodge a new issue.
