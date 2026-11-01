# Deprecated: is.MVG

**\[deprecated\]** `is.MVG()` was deprecated in MVBeliefUpdatr 0.1.0 and
will be removed in 0.2.0. Please use S7 validators and
[MVG_IdealObserver](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md)
instead.

## Usage

``` r
is.MVG(x, group = NULL, category = "category", is.long = T, verbose = F)
```

## Arguments

- x:

  Object to be checked.

- group:

  Name of one or more group variables, each unique combination of which
  describes an MVG. (default: NULL)

- category:

  DEPRECATED Name of the category variable. (default: "category")

- is.long:

  Is this check assessing whether the ideal observer is in long format
  (`TRUE`) or wide format (`FALSE`)? (default: `TRUE`)

- verbose:

  Logical. If `TRUE`, emits diagnostics.

## Value

A logical.

## See also

[MVG_IdealObserver](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md)
