# Parametric bootstrap for an asreml model.

Simulate \\\hat{y} \sim N(X\hat{\beta}, V)\\ according to the current
asreml fit and then refit to obtain the targeted statistics.

## Usage

``` r
bootstrap_asreml(
  model,
  FUN,
  nsim = 1,
  use.u = FALSE,
  source = list(),
  seed = NULL,
  ...
)
```

## Arguments

- model:

  An `asreml` fitted model. Must be fitted with `model.frame = TRUE`.

- FUN:

  A function with signature `function(fit)` returning a scalar (the
  statistic to bootstrap).

- nsim:

  Integer. Number of bootstrap replicates.

- use.u:

  A logical indicating whether to resample random effects, or only
  resample residuals.

- source:

  The known genomic relationship matrix (GRM) used in `model` fitted
  using `asreml::vm()`, provided as a named list. When not provided (an
  empty list by default), the GRM variable used for `vm` calling will be
  searched in the global environment.

- seed:

  Optional integer seed for reproducibility.

- ...:

  Additional arguments passed to
  [`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html).

## Value

A `boot` object.

## Details

Fits parametric bootstrap replicates for an `asreml` model by:

- Extracting the fixed-effect fit yhat = X \* beta.

- Extracting V = Var(y) on the observation scale,

- Simulating new responses y\* = yhat + L %\*% z where L is a Cholesky
  factor of V,

- Refitting the same asreml call on each simulated dataset,

- Returning a `boot` object.

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

  b <- bootstrap_asreml(
    lettuce_asreml,
    R = 200,
    statistic = function(fit) coef(fit)$fixed["(Intercept)", "effect"],
    seed = 1
  )
  boot::boot.ci(b, type = "perc")
 } # }
```
