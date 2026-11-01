# Deprecated: inv_symlog

**\[deprecated\]** `inv_symlog()` was deprecated in MVBeliefUpdatr 0.1.0
and will be removed in 0.2.0. Please use
[`scales::pseudo_log_trans()`](https://scales.r-lib.org/reference/transform_log.html)
instead.

## Usage

``` r
inv_symlog(x, C = 0)
```

## Arguments

- x:

  Numeric vector to invert.

- C:

  Scaling constant. Defaults to 0.

## Value

Inverse transformed numeric vector.

## See also

[`scales::pseudo_log_trans()`](https://scales.r-lib.org/reference/transform_log.html)
