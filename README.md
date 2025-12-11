

<!-- README.md is generated from README.Rmd. Please edit that file -->

# heritable

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/heritable.png)](https://CRAN.R-project.org/package=heritable)
[![Codecov test coverage](https://codecov.io/gh/anu-aagi/heritable/graph/badge.svg)](https://app.codecov.io/gh/anu-aagi/heritable)
[![R-CMD-check](https://github.com/anu-aagi/heritable/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/anu-aagi/heritable/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of heritable is to be the one stop shop for heritability
calculations in R. We currently support model outputs from `asreml` and
support the calculations of broad-sense heritability using the methods
of Cullis et al (2006), Oakey et al (2006) and Schmidt et al (2019).

## Installation

You can install the development version of heritable from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("anu-aagi/heritable")
```

## Example

This is a basic example which shows you how to solve a common problem:

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

``` r
H2(fit_asreml, target = "gen")
#>    Cullis     Oakey    Piepho     Delta  Standard 
#> 0.8090841 0.8090841 0.8029760 0.8090841 0.8400648
H2(fit_lme4, target = "gen")
#>    Cullis     Oakey    Piepho     Delta  Standard 
#> 0.8091338 0.8091338 0.7966375 0.8091338 0.8400678
```

Working as data using tidyverse style of data wrangling.

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
