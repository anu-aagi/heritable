# Pull fixed and random terms from a model formula

Extract the labels of fixed and random terms from a model object that
exposes a formula with `fixed` and `random` components (for example
objects produced by asreml::asreml). The function returns a named list
containing two character vectors: `fixed` and `random`.

## Usage

``` r
# S3 method for class 'asreml'
pull_terms(model)
```

## Arguments

- model:

  A fitted model object with a `formula` method that returns a list
  containing `fixed` and `random` formula components.

## Value

A named list with components:

- fixed:

  Character vector of labels for fixed-effect terms.

- random:

  Character vector of labels for random-effect terms.
