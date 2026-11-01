# Deprecated: get_posterior_predictive_from_NIW_belief

**\[deprecated\]** `get_posterior_predictive_from_NIW_belief()` was
deprecated in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please
use
[`likelihood()`](https://hlplab.github.io/MVBeliefUpdatr/reference/likelihood.md)
instead.

**\[deprecated\]** `get_posterior_predictives_from_NIW_beliefs()` was
deprecated in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please
use
[`likelihood()`](https://hlplab.github.io/MVBeliefUpdatr/reference/likelihood.md)
instead.

## Usage

``` r
get_posterior_predictive_from_NIW_belief(
  x,
  model,
  noise_treatment = "no_noise",
  log = TRUE,
  category = "category",
  category.label = NULL,
  wide = FALSE
)

get_posterior_predictives_from_NIW_beliefs(
  x,
  model,
  noise_treatment = "no_noise",
  log = TRUE,
  category = "category",
  category.label = NULL,
  grouping.var = NULL,
  wide = FALSE
)
```

## Arguments

- x:

  Observations.

- model:

  An NIW belief or adaptor model.

- noise_treatment:

  Noise treatment.

- log:

  Logical; whether log likelihood is returned.

- category:

  Category column name.

- category.label:

  Category labels.

- wide:

  Logical; whether wide format is returned.

- grouping.var:

  Grouping variable name.

## Value

Posterior predictive data frame.

## See also

[`likelihood()`](https://hlplab.github.io/MVBeliefUpdatr/reference/likelihood.md)
