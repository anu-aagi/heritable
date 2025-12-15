# Fit the counterpart of an asreml model by swapping a specified term between fixed and random

Fit a "counterpart" model to an existing asreml model by moving a
specified term from the random effects to the fixed effects or vice
versa. This is useful for procedures that require the same term to be
fitted as the opposite effect (for example, calculating Piepho's
heritability).

## Usage

``` r
# S3 method for class 'asreml'
fit_counterpart_model(model, target = NULL)
```

## Arguments

- model:

  An existing fitted asreml model object. The function expects that
  formulas can be retrieved via formula(model)\$fixed and
  formula(model)\$random and that the model frame is available as
  model\$mf.

- target:

  Character(1). Name of the term (e.g. a factor variable) to be switched
  between fixed and random effects. Must match one of the terms present
  in either the model's fixed or random formulas.

## Value

A fitted asreml model object identical to the input model except that
`target` has been moved from the random effects to the fixed effects or
from the fixed effects to the random effects, depending on where it
appeared in the original model. The returned object is produced by
asreml::asreml and can be used for subsequent model comparisons or
variance-component calculations.
