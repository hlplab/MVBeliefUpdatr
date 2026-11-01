# Deprecated: is.exemplar_model

**\[deprecated\]** `is.exemplar_model()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use S7
validators and
[Exemplar_Model](https://hlplab.github.io/MVBeliefUpdatr/reference/family-exemplar.md)
instead.

## Usage

``` r
is.exemplar_model(x, group = NULL, verbose = F, tolerance = MVBU_PROB_TOL)
```

## Arguments

- x:

  Object to be checked.

- group:

  Name of one or more group variables, each unique combination of which
  describes an exemplar model. (default: NULL)

- verbose:

  Should verbose output be provided? (default: `TRUE`)

- tolerance:

  Probability tolerance.

## Value

A logical.

## See also

[Exemplar_Model](https://hlplab.github.io/MVBeliefUpdatr/reference/family-exemplar.md)
