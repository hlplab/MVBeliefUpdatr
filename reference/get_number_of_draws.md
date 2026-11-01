# Get number of post-warmup MCMC samples from stanfit

Get number of post-warmup MCMC samples from stanfit

Get indices for random MCMC draws from stanfit

## Usage

``` r
get_number_of_draws(fit, ...)

get_random_draw_indices(fit, ndraws = NULL, ...)

## S7 method for class <any>
get_number_of_draws(fit, ...)

## S7 method for class <any>
get_random_draw_indices(fit, ndraws = NULL, ...)
```

## Arguments

- fit:

  A [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html)
  object or
  [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
  object.

- ...:

  Additional arguments passed to methods.

- ndraws:

  Number of indices to be returned.

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
get_number_of_draws(fit)
#> [1] 25
get_random_draw_indices(fit, ndraws = 5)
#> [1] 11 10 24  1 19
# }
```
