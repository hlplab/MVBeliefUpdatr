# Get sufficient category statistics from a data frame

Calculates sufficient statistics for specified cues for all categories
and optional groups in a data frame.

## Usage

``` r
get_sufficient_category_statistics(
  data,
  cues,
  category = "category",
  group = NULL,
  stats = c("n", "mean", "cov", "css", "uss"),
  simplify = TRUE,
  ...
)
```

## Arguments

- data:

  A `data.frame` or `tibble`.

- cues:

  Character vector of column names for cues. Required.

- category:

  Character string specifying category column name. Default:
  `"category"`.

- group:

  Optional character vector of group column name(s). Default: `NULL`.

- stats:

  Character vector of statistics to return. Allowed values are `"n"`,
  `"mean"`, `"cov"`, `"css"`, and `"uss"`. Default:
  `c("n", "mean", "cov", "css", "uss")`.

- simplify:

  Logical. If `TRUE` (default) and the data matches a single category
  and group with a single statistic, returns an atomic scalar, vector,
  or matrix. If `FALSE`, always returns a `tibble`.

- ...:

  Additional arguments (e.g. `model_family`).

## Value

If `simplify = TRUE` and a single statistic is selected for a single
category and group, the atomic numeric value, vector, or matrix.
Otherwise, a `tibble` with the grouping columns and requested

## See also

[`get_exposure_category_statistics`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_exposure_category_statistics.md)

## Examples

``` r
# Example with synthetic phonetic data
df <- data.frame(
  category = rep(c("b", "p"), each = 20),
  vot = c(rnorm(20, 0, 10), rnorm(20, 50, 10)),
  f0 = c(rnorm(20, 200, 20), rnorm(20, 240, 20))
)
get_sufficient_category_statistics(df, cues = c("vot", "f0"))
#> # A tibble: 2 × 6
#>   category     n mean      cov           css           uss          
#>   <chr>    <int> <list>    <list>        <list>        <list>       
#> 1 b           20 <dbl [2]> <dbl [2 × 2]> <dbl [2 × 2]> <dbl [2 × 2]>
#> 2 p           20 <dbl [2]> <dbl [2 × 2]> <dbl [2 × 2]> <dbl [2 × 2]>
```
