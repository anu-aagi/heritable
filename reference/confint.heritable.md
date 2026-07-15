# Bootstrap confidence interval for heritability

Computes a confidence interval for a heritability estimate using
parametric bootstrap of the underlying mixed model. Parallel computing
is supported through
[`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html)

## Usage

``` r
# S3 method for class 'heritable'
confint(
  object,
  parm = NULL,
  level = 0.95,
  B = 100,
  random_effect = c("resample", "conditional"),
  type = c("basic", "norm", "perc"),
  return_model = TRUE,
  seed = NULL,
  ...
)
```

## Arguments

- object:

  A heritability object returned by
  [`H2()`](https://anu-aagi.github.io/heritable/reference/H2.md)
  (broad-sense)
  [h2()](https://anu-aagi.github.io/heritable/reference/H2.md)
  (narrow-sense). The object must store the fitted model as an
  attribute.

- parm:

  a specification of which methods are to be given confidence intervals,
  either a vector of numbers or a vector of names. If missing, all
  methods are considered.

- level:

  Confidence level.

- B:

  Integer. Number of bootstrap replicates.

- random_effect:

  Character. Strategy for handling random effects.

  `"resample"`

  :   Resample random effects to propagate uncertainty.

  `"conditional"`

  :   Condition on estimated random effects.

- type:

  Character. Bootstrap interval type; one of `"basic"`, `"norm"`, or
  `"perc"`.

- return_model:

  Logical. Whether to return to the `boot` object.

- seed:

  Optional random seed.

- ...:

  Additional arguments passed to the bootstrap routine. Check
  [`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html), as well as
  the examples below for parallel computation

## Value

A matrix of confidence intervals.

## See also

[`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html),
[`H2()`](https://anu-aagi.github.io/heritable/reference/H2.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  lettuce_subset <- lettuce_phenotypes |> subset(loc == "L2")
  lettuce_asreml <- asreml(
    fixed = y ~ rep * pseudo_var1,
    random = ~gen,
    sparse = ~pseudo_var2,
    data = lettuce_subset,
    trace = FALSE
  )

  my_H2 <- H2(lettuce_asreml, "gen", c("Cullis", "Standard"))

  my_ci <- confint(my_H2)

  # Get bootstrap model
  boot_mod <- attr(my_ci, "boot_mod")
  boot_mod$t # Check bootstrap statistics

  # Parallel computing (On Windows)
  # Note that for asreml, `ncpus` can't be larger than the number of asreml
  # license available
  confint(my_H2, parallel = "snow", ncpus = 3)

  # Parallel computing (On non-Windows)
  confint(my_H2, parallel = "multicore", ncpus = 3)

} # }
```
