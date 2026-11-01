# Categorize one or more observations under the model's decision rule

Evaluates the model's response category decisions and response
probabilities.

## Usage

``` r
categorize(x, new_data, decision_rule = "criterion", simplify = NULL, ...)

## S7 method for classes <MVBU_Stanfit>/<MVBU_StanfitPosterior>, <any>, <any>
categorize(x, new_data, decision_rule, simplify = NULL, ...)
```

## Arguments

- x:

  A cognitive model object,
  [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md),
  or
  [`MVBU_StanfitPosterior`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_StanfitPosterior.md).

- new_data:

  A numeric matrix of observations, data frame, or list of observation
  matrices.

- decision_rule:

  An optional override for the model's decision rule (`"criterion"`,
  `"sample"` / `"sampling"`, or `"proportional"`). Defaults to the
  model's configured decision rule.

- simplify:

  Logical scalar; if `TRUE`, returns a character vector of chosen
  `response_category` values. If `FALSE`, returns a data frame with
  columns `response_category` and `response_probability`. Defaults to
  `FALSE` when `decision_rule == "proportional"` and `TRUE` otherwise.

- ...:

  Additional arguments passed to methods.

## Value

If `simplify = TRUE`, a character vector of response categories. If
`simplify = FALSE`, a data frame containing columns `response_category`
and `response_probability` (which equals 1 for deterministic or sampled
decisions, and equals the category posterior probability under
proportional decision rules).

## Examples

``` r
mod <- example_mvg_ideal_observer(n_cues = 2)
categorize(mod, matrix(c(0, 0), nrow = 1))
#> [1] "/b/"
categorize(mod, matrix(c(0, 0), nrow = 1), simplify = FALSE)
#>   response_category response_probability
#> 1               /b/                    1
```
