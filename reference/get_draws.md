# Get MCMC prior or posterior draws from a stanfit object

Get MCMC draws of all parameters from incremental Bayesian
belief-updating (IBBU) as a tibble in long format. By default all
post-warmup draws are returned, but if `summarize = TRUE` then just the
mean of each parameter is returned instead.

## Usage

``` r
get_draws(
  fit,
  categories = NULL,
  groups = NULL,
  which = "posterior",
  ndraws = NULL,
  summarize = FALSE,
  nest = TRUE,
  seed = NULL,
  ...
)

## S7 method for class <MVBU_StanfitPosterior>
get_draws(
  fit,
  categories = NULL,
  groups = NULL,
  which = "posterior",
  ndraws = NULL,
  summarize = FALSE,
  nest = TRUE,
  seed = NULL,
  ...
)

## S7 method for class <MVBU_Stanfit>
get_draws(
  fit,
  categories = NULL,
  groups = NULL,
  which = "posterior",
  ndraws = NULL,
  summarize = FALSE,
  nest = TRUE,
  seed = NULL,
  ...
)
```

## Arguments

- fit:

  An
  [`IdealAdaptorStanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
  (or
  [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md))
  object.

- categories:

  Character vector of category names for which draws should be returned.
  (default: all categories in model)

- groups:

  Character vector of group names for which draws should be returned.
  (default: all groups in model including prior)

- which:

  DEPRECATED. Use `groups` instead. Should parameters for the prior,
  posterior, or both be returned? (default: `"posterior"`)

- ndraws:

  Number of random draws or `NULL` if all draws are to be returned.
  (default: `NULL`)

- summarize:

  Should the mean of the draws be returned instead of all of the draws?
  (default: `FALSE`)

- nest:

  Should the category mean vectors and scatter matrices be nested into
  one cell each, or should each element be stored in a separate column?
  (default: `TRUE`)

- seed:

  Optional seed for reproducibility when sampling draws. (default:
  `NULL`)

- ...:

  Additional arguments passed to methods.

## Value

A tibble with MCMC draws of the model parameters.

## Details

By default, the category means and scatter matrices are nested, rather
than each of their elements being stored separately (`nest = TRUE`).

## See also

[`get_parameters`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_parameters.md),
[`get_stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_stanfit.md)

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
draws <- get_draws(fit, ndraws = 5)
head(draws)
#> # A tibble: 6 × 12
#>   .chain .iteration .draw group category kappa    nu m      S        lapse_rate
#>    <int>      <int> <int> <fct> <fct>    <dbl> <dbl> <list> <list>        <dbl>
#> 1      1          2     2 prior /d/       13.2  26.3 <dbl>  <dbl[…]>     0.0277
#> 2      1          2     2 prior /t/       13.2  26.3 <dbl>  <dbl[…]>     0.0277
#> 3      1          3     3 prior /d/       25.4  56.2 <dbl>  <dbl[…]>     0.0223
#> 4      1          3     3 prior /t/       25.4  56.2 <dbl>  <dbl[…]>     0.0223
#> 5      1         16    16 prior /d/       10.2  19.8 <dbl>  <dbl[…]>     0.0762
#> 6      1         16    16 prior /t/       10.2  19.8 <dbl>  <dbl[…]>     0.0762
#> # ℹ 2 more variables: Sigma_exp <dbl>, Sigma_marg <list>
# }
```
