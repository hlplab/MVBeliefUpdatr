# Read and write cached Stanfit objects to disk

Reads or saves a previously fitted
[`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
object (or subclass such as
[`IdealAdaptorStanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md))
from or to an RDS file.

## Usage

``` r
read_stanfit(file)

write_stanfit(x, file, compress = TRUE)
```

## Arguments

- file:

  File path to the cached fit. If the file extension is omitted, `.rds`
  is appended automatically.

- x:

  An object inheriting from
  [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md).

- compress:

  Logical or character string specifying compression passed to
  [`saveRDS`](https://rdrr.io/r/base/readRDS.html) (default: `TRUE`).

## Value

For `read_stanfit`, the loaded
[`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
object with its `@file` slot updated, or `NULL` if the file does not
exist or cannot be read. For `write_stanfit`, the saved object `x`
invisibly with its `@file` slot updated.

## See also

[`fit_ideal_adaptor`](https://hlplab.github.io/MVBeliefUpdatr/reference/fit_ideal_adaptor.md)

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
tmp <- tempfile(fileext = ".rds")
write_stanfit(fit, tmp)
fit_loaded <- read_stanfit(tmp)
is(fit_loaded, "MVBU_Stanfit")
#> [1] TRUE
unlink(tmp)
# }
```
