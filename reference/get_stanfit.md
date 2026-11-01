# Get the stored Stan fit from an S7 fit object

Get the stored Stan fit from an S7 fit object

## Usage

``` r
get_stanfit(x, ...)

set_stanfit(x, stanfit, ...)

## S7 method for class <any>
get_stanfit(x, ...)

## S7 method for class <MVBU_Stanfit>
get_stanfit(x, ...)

## S7 method for classes <any>, <any>
set_stanfit(x, stanfit, ...)

## S7 method for classes <MVBU_Stanfit>, <any>
set_stanfit(x, stanfit, ...)
```

## Arguments

- x:

  An S7 fit object.

- ...:

  Additional arguments passed to methods.

- stanfit:

  An
  [rstan::stanfit](https://mc-stan.org/rstan/reference/stanfit-class.html)
  object.

## Value

The stored
[rstan::stanfit](https://mc-stan.org/rstan/reference/stanfit-class.html)
object when present, or the updated S7 fit object.

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
sf <- get_stanfit(fit)
fit <- set_stanfit(fit, sf)
# }
```
