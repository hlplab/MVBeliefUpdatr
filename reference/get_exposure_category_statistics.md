# Get exposure category statistics from a Stanfit object or Stanfit input

Extracts or computes exposure category statistics (n, mean, cov, css,
uss) from an `IdealAdaptorStanfitInput`, `MVBU_Stanfit`, or
`MVBU_StanfitPosterior` object.

## Usage

``` r
get_exposure_category_statistics(
  x,
  categories = NULL,
  groups = NULL,
  stats = c("n", "mean", "cov", "css", "uss"),
  simplify = TRUE,
  ...
)

## S7 method for class <IdealAdaptorStanfitInput>/<MVBU_Stanfit>/<MVBU_StanfitPosterior>
get_exposure_category_statistics(
  x,
  categories = NULL,
  groups = NULL,
  stats = c("n", "mean", "cov", "css", "uss"),
  simplify = TRUE,
  ...
)
```

## Arguments

- x:

  An `IdealAdaptorStanfitInput`, `MVBU_Stanfit`, or
  `MVBU_StanfitPosterior` object.

- categories:

  Optional character vector of category names. Default: `NULL`.

- groups:

  Optional character vector of group names. Default: `NULL`.

- stats:

  Character vector of statistics to return. Allowed values are `"n"`,
  `"mean"`, `"cov"`, `"css"`, and `"uss"`. Default:
  `c("n", "mean", "cov", "css", "uss")`.

- simplify:

  Logical. If `TRUE` (default) and the query matches a single category,
  group, and statistic, returns an atomic vector or matrix. If `FALSE`,
  always returns a `tibble`.

- ...:

  Additional arguments passed to methods.

## Value

If `simplify = TRUE` and a single statistic is selected for a single
category and group, the atomic numeric value, vector, or matrix.
Otherwise, a `tibble` with the grouping columns and requested
statistics.

## See also

[`get_sufficient_category_statistics`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_sufficient_category_statistics.md)

## Examples

``` r
sfi <- example_ideal_adaptor_stanfit_input()
get_exposure_category_statistics(sfi, stats = "mean")
#> # A tibble: 4 × 3
#>   group         category mean     
#>   <fct>         <fct>    <list>   
#> 1 left_shifted  /d/      <dbl [1]>
#> 2 right_shifted /d/      <dbl [1]>
#> 3 left_shifted  /t/      <dbl [1]>
#> 4 right_shifted /t/      <dbl [1]>
```
