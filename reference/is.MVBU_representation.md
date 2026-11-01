# Deprecated: is.MVBU_representation

**\[deprecated\]** `is.MVBU_representation()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use S7
predicates and classes instead.

## Usage

``` r
is.MVBU_representation(x, group = NULL, verbose = F, tolerance = MVBU_PROB_TOL)
```

## Arguments

- x:

  Object to be checked.

- group:

  Name of one or more group variables, each unique combination of which
  describes a model. (default: NULL)

- verbose:

  Should verbose output be provided? (default: `TRUE`)

- tolerance:

  Probability tolerance.

## Value

A logical.

## See also

[MVBU_CategoryRepresentationTemplate](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md)
