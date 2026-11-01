# Deprecated: get_limits

**\[deprecated\]** `get_limits()` was deprecated in MVBeliefUpdatr 0.1.0
and will be removed in 0.2.0.

## Usage

``` r
get_limits(data, measure, by = NULL, hdi.prob = 0.99, min = NULL, max = NULL)
```

## Arguments

- data:

  A `tibble` or `data.frame`.

- measure:

  Variable name.

- by:

  Optional grouping variable.

- hdi.prob:

  HDI probability mass.

- min:

  Optional minimum.

- max:

  Optional maximum.

## Value

Numeric vector of length 2.
