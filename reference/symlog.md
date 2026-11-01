# Deprecated: symlog

**\[deprecated\]** `symlog()` was deprecated in MVBeliefUpdatr 0.1.0 and
will be removed in 0.2.0. Please use
[`scales::pseudo_log_trans()`](https://scales.r-lib.org/reference/transform_log.html)
instead.

## Usage

``` r
symlog(x, C = 0)
```

## Arguments

- x:

  Numeric vector to transform.

- C:

  Scaling constant determining resolution around zero. Defaults to 0.

## Value

Transformed numeric vector.

## See also

[`scales::pseudo_log_trans()`](https://scales.r-lib.org/reference/transform_log.html)
