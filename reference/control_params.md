# Extract control parameters of the NUTS sampler

Extract control parameters of the NUTS sampler

## Usage

``` r
control_params(x, pars = NULL, ...)

## S7 method for class <MVBU_Stanfit>
control_params(x, pars = NULL, ...)
```

## Arguments

- x:

  An S7 fit object or stanfit object.

- pars:

  Optional parameter names.

- ...:

  Additional arguments.

## Value

A named list of control parameters.

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
control_params(fit)
#> $adapt_delta
#> [1] 0.8
#> 
#> $adapt_engaged
#> [1] TRUE
#> 
#> $adapt_gamma
#> [1] 0.05
#> 
#> $adapt_init_buffer
#> [1] 75
#> 
#> $adapt_kappa
#> [1] 0.75
#> 
#> $adapt_t0
#> [1] 10
#> 
#> $adapt_term_buffer
#> [1] 50
#> 
#> $adapt_window
#> [1] 25
#> 
#> $max_treedepth
#> [1] 10
#> 
#> $metric
#> [1] "diag_e"
#> 
#> $stepsize
#> [1] 1
#> 
#> $stepsize_jitter
#> [1] 0
#> 
# }
```
