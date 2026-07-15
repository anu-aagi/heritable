# Extract variance components

This function provides the variance matrix quantities needed to evaluate
heritability with respect to a target random-effect term in a fitted
linear mixed model.

## Usage

``` r
var_comp(
  model,
  target,
  calc_C22 = TRUE,
  calc_V = TRUE,
  calc_C11 = TRUE,
  marginal = TRUE,
  stratification = NULL,
  solver = c("direct", "LMM"),
  ...
)
```

## Arguments

- model:

  Model object of class `lmerMod/merMod` or `asreml`

- target:

  The name of the random effect for which heritability is to be
  calculated.

- calc_C22:

  Logical; whether to compute the prediction error variance (PEV) matrix
  for the transformed target effect.

- calc_V:

  Logical; whether to retain the marginal variance matrix of the
  response and supporting quantities (`V`, `G`, `Z`, `X`, `idx`, and
  `m`) in the output.

- calc_C11:

  Logical; whether to compute the variance matrix of the fixed-effect
  counterpart estimator for the transformed target effect.

- marginal:

  Logical; if `TRUE`, construct marginal (strata-averaged) mappings so
  that each genotype receives a single averaged effect per term. If
  `FALSE`, mappings will only consider the main genotype effect and
  ignore the iteracting terms.

- stratification:

  A one-row data frame defining the stratum in which genotype effects
  should be evaluated. The columns must correspond to model terms that
  interact with `target`.

- solver:

  A string specifying the solver for the PEV matrix. Can be either
  `"direct"` (directly invert `V`) or `"LMM"` (Solve the LMM equation).

- ...:

  Additional arguments passed to downstream helper functions.

## Value

A named list with the following elements:

- `n_g`:

  Number of transformed genetic coefficients after applying the mapping
  matrix `m`.

- `gnames`:

  Names of the transformed genetic coefficients.

- `G_g`:

  Variance matrix of the transformed target genetic effect.

- `C22_g`:

  Prediction error variance matrix of the transformed target effect, if
  `calc_C22 = TRUE`; otherwise `NULL`.

- `G_g_tilde`:

  Variance matrix of the fixed-effect counterpart estimator for the
  transformed target effect, if `calc_C11 = TRUE`; otherwise `NULL`.

- `G_g_tilde_no_cov`:

  Variance matrix of the fixed-effect counterpart estimator for the
  transformed target effect, without considering target covariance, if
  `calc_C11 = TRUE`; otherwise `NULL`.

- `C11_g`:

  Estimation error variance matrix of the fixed-effect counterpart
  estimator for the transformed target effect, if `calc_C11 = TRUE`;
  otherwise `NULL`.

- `V`:

  Marginal covariance matrix of the response, if `calc_V = TRUE`;
  otherwise `NULL`.

- `G`:

  Variance matrix of all random effects, if `calc_V = TRUE`; otherwise
  `NULL`.

- `Z`:

  Random-effect design matrix, if `calc_V = TRUE`; otherwise `NULL`.

- `X`:

  Fixed-effect design matrix, if `calc_V = TRUE`; otherwise `NULL`.

- `idx`:

  Indices of the random-effect coefficients associated with the target
  term, if `calc_V = TRUE`; otherwise `NULL`.

- `W`:

  Linear mapping from the original target coefficients to the
  transformed target effect, if `calc_V = TRUE`; otherwise `NULL`.

- `marginal`:

  Logical scalar indicating whether the returned quantities correspond
  to a marginal definition of the target effect.

- `stratification`:

  The user-supplied stratification object, if any.
