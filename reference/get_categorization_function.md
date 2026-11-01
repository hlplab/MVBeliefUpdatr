# Get categorization function from a model object

Returns a categorization function that can be used to categorize a set
of cues into categories.

## Usage

``` r
get_categorization_function(
  x,
  lapse_treatment = c("no_lapses", "sample", "marginalize")[3],
  groups = get_group_labels(x, include_prior = TRUE),
  ...
)
```

## Arguments

- x:

  A model object.

- lapse_treatment:

  Should the consequences of attentional lapses be included in the
  categorization function ("marginalize") or not ("no_lapses")?
  (default: "marginalize")

- groups:

  Character vector of groups to include.

- ...:

  Optionally, additional arguments handed to
  [`get_draws`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_draws.md).

## Value

A tibble with categorization functions per group and draw.

## Examples

``` r
# \donttest{
fit <- example_ideal_adaptor_stanfit()
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
cat_fun <- get_categorization_function(fit)
#> Warning: `get_NIW_categorization_function()` was deprecated in MVBeliefUpdatr 0.1.0.
#> ℹ Please use `get_category_posterior_function()` instead.
#> ℹ The deprecated feature was likely used in the MVBeliefUpdatr package.
#>   Please report the issue at <https://github.com/hlplab/MVBeliefUpdatr/issues>.
head(cat_fun)
#> # A tibble: 6 × 3
#> # Groups:   group, .draw [6]
#>   group .draw f     
#>   <fct> <int> <list>
#> 1 prior     1 <fn>  
#> 2 prior     2 <fn>  
#> 3 prior     3 <fn>  
#> 4 prior     4 <fn>  
#> 5 prior     5 <fn>  
#> 6 prior     6 <fn>  
# }
```
