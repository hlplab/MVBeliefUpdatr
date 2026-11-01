# Deprecated: get_constructor

**\[deprecated\]** `get_constructor()`, `get_category_constructor()`,
`get_group_constructor()`, `get_cue_constructor()`, and
`get_cue2_constructor()` are deprecated in favor of
[`get_category_labels`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_labels.md),
[`get_group_labels`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_labels.md),
and
[`get_cue_labels`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_labels.md).

## Usage

``` r
get_constructor(x, variable = NULL)

get_category_constructor(x)

get_group_constructor(x)

get_cue_constructor(x)

get_cue2_constructor(x)
```

## Arguments

- x:

  An
  [`IdealAdaptorStanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
  or
  [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
  object.

- variable:

  Either "category", "group", "cue", or "cue2". If set to `NULL` then a
  list of all constructors is returned. (default: `NULL`)

## Value

A constructor function, a list of constructor functions, or `NULL`.
