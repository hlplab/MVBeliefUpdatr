# Deprecated: is.NIW_belief

**\[deprecated\]** `is.NIW_belief()` was deprecated in MVBeliefUpdatr
0.1.0 and will be removed in 0.2.0. Please use S7 validators and
[NIW_IdealAdaptor](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md)
instead.

## Usage

``` r
is.NIW_belief(x, group = NULL, category = "category", is.long = T, verbose = F)
```

## Arguments

- x:

  Object to be checked.

- group:

  Name of one or more group variables, each unique combination of which
  describes an NIW_belief. (default: NULL)

- category:

  Name of the category variable. (default: "category")

- is.long:

  Is this check assessing whether the belief is in long format (`TRUE`)
  or wide format (`FALSE`)? (default: `TRUE`)

- verbose:

  Logical. If `TRUE`, emits diagnostics.

## Value

A logical.

## See also

[NIW_IdealAdaptor](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md)
