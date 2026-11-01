# Deprecated: is.NIW_ideal_adaptor

**\[deprecated\]** `is.NIW_ideal_adaptor()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use S7
validators and
[NIW_IdealAdaptor](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md)
instead.

## Usage

``` r
is.NIW_ideal_adaptor(
  x,
  group = NULL,
  category = "category",
  is.long = T,
  with.prior = T,
  with.lapse = if (with.lapse_bias) T else F,
  with.lapse_bias = F,
  verbose = F,
  tolerance = MVBU_PROB_TOL,
  ...
)
```

## Arguments

- x:

  Object to check.

- group:

  Name of one or more group variables. (default: NULL)

- category:

  Name of the category variable. (default: "category")

- is.long:

  Is this check assessing whether the ideal adaptor is in long format?
  (default: `TRUE`)

- with.prior:

  Does this ideal adaptor have a prior? (default: `TRUE`)

- with.lapse:

  Does this ideal adaptor have a lapse rate? (default: `FALSE`)

- with.lapse_bias:

  Does this ideal adaptor have a lapse bias? (default: `FALSE`)

- verbose:

  Should verbose output be provided? (default: `FALSE`)

- tolerance:

  Tolerance for sum-to-one probability checks. (default: MVBU_PROB_TOL)

- ...:

  Additional arguments.

## Value

Logical indicating whether `x` is a valid NIW ideal adaptor.

## See also

[NIW_IdealAdaptor](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md)
