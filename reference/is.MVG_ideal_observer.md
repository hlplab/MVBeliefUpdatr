# Deprecated: is.MVG_ideal_observer

**\[deprecated\]** `is.MVG_ideal_observer()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use S7
validators and
[MVG_IdealObserver](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md)
instead.

## Usage

``` r
is.MVG_ideal_observer(
  x,
  group = NULL,
  category = "category",
  is.long = T,
  with.lapse = if (with.lapse_bias) T else F,
  with.lapse_bias = F,
  verbose = F,
  tolerance = MVBU_PROB_TOL
)
```

## Arguments

- x:

  Object to be checked.

- group:

  Name of one or more group variables, each unique combination of which
  describes an MVG_ideal_observer. (default: NULL)

- category:

  Name of the category variable. (default: "category")

- is.long:

  Is this check assessing whether the ideal observer is in long format
  (`TRUE`) or wide format (`FALSE`)? (default: `TRUE`)

- with.lapse:

  Does this ideal observer have a lapse rate? (default: `FALSE`)

- with.lapse_bias:

  Does this ideal observer have a lapse bias? (default: `FALSE`)

- verbose:

  Should verbose output be provided? (default: `TRUE`)

- tolerance:

  Probability tolerance.

## Value

A logical.

## See also

[MVG_IdealObserver](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md)
