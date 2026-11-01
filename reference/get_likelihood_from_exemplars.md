# Deprecated: get_likelihood_from_exemplars

**\[deprecated\]** `get_likelihood_from_exemplars()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`likelihood()`](https://hlplab.github.io/MVBeliefUpdatr/reference/likelihood.md)
instead.

## Usage

``` r
get_likelihood_from_exemplars(
  x,
  model,
  noise_treatment = "no_noise",
  log = TRUE,
  category = "category",
  category.label = NULL
)
```

## Arguments

- x:

  Observations.

- model:

  Exemplar model.

- noise_treatment:

  Noise treatment.

- log:

  Logical; whether log likelihood is returned.

- category:

  Category column name.

- category.label:

  Category labels.

## Value

Likelihood data frame.

## See also

[`likelihood()`](https://hlplab.github.io/MVBeliefUpdatr/reference/likelihood.md)
