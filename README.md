<!-- README.md is generated from README.Rmd. Please edit that file -->

# heritable

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN status](https://www.r-pkg.org/badges/version/heritable)](https://CRAN.R-project.org/package=heritable)
<!-- badges: end -->

The goal of heritable is to be the one stop shop for heritability calculations in R. We currently support model outputs from `asreml` and support the calculations of broad-sense heritability using the methods of Cullis et al (2006), Oakey et al (2006) and Schmidt et al (2019).

## Installation

You can install the development version of heritable from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("anu-aagi/heritable")
```

## Example

This is a basic example which shows you how to solve a common problem:


``` r
library(heritable)  

fit <- asreml::asreml(yield ~ rep,
    random = ~ gen + rep:block,
    data = agridat::john.alpha,
    trace = FALSE
  )

 H2(fit, target = "gen") 
```

