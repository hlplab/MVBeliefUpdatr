# Deprecated: Get sum of squares from a data frame or matrix

**\[deprecated\]** `get_sum_of_squares_from_df()` and its aliases are
deprecated as of MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0.
Please use
[`get_sufficient_category_statistics()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_sufficient_category_statistics.md)
instead.

## Usage

``` r
get_sum_of_squares_from_df(
  data,
  variables = NULL,
  center = TRUE,
  verbose = FALSE
)

ss(x, center = TRUE)

get_sum_of_uncentered_squares_from_df(data, variables = NULL, verbose = FALSE)

get_sum_of_centered_squares_from_df(data, variables = NULL, verbose = FALSE)
```

## Arguments

- data:

  A `tibble`, `data.frame`, or `matrix`.

- variables:

  Only required if data is not already a `matrix`.

- center:

  Logical; if `TRUE`, compute centered sum of squares.

- verbose:

  Logical; produce verbose output.

## Value

A matrix.
