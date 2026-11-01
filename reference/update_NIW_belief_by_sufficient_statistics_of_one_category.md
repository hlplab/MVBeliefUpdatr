# Deprecated: update_NIW_belief_by_sufficient_statistics_of_one_category

**\[deprecated\]**
`update_NIW_belief_by_sufficient_statistics_of_one_category()` was
deprecated in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please
use
[`update_category_representation()`](https://hlplab.github.io/MVBeliefUpdatr/reference/update_category_representation.md)
instead.

## Usage

``` r
update_NIW_belief_by_sufficient_statistics_of_one_category(
  prior_model,
  x_category = NULL,
  x_mean,
  x_SS,
  x_N,
  ...
)
```

## Arguments

- prior_model:

  An NIW S7 category representation.

- x_category:

  Deprecated and ignored for S7 representations.

- x_mean:

  Observation mean.

- x_SS:

  Centered observation sum-of-squares matrix.

- x_N:

  Number of observations.

- ...:

  Deprecated compatibility arguments.

## Value

An updated NIW category representation.

## See also

[`update_category_representation()`](https://hlplab.github.io/MVBeliefUpdatr/reference/update_category_representation.md)
