# Get the stored Stan input from an S7 fit or fit-input object

Get the stored Stan input from an S7 fit or fit-input object

## Usage

``` r
get_staninput(x, ...)

set_staninput(x, staninput, ...)

## S7 method for class <any>
get_staninput(x, ...)

## S7 method for class <MVBU_Stanfit>
get_staninput(x, ...)
```

## Arguments

- x:

  An S7 fit or fit-input object.

- ...:

  Additional arguments passed to methods.

- staninput:

  An S7 Stan input object.

## Value

The stored S7 Stan input object, or the updated S7 object.

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
si <- get_staninput(fit)
fit <- set_staninput(fit, si)
# }
```
