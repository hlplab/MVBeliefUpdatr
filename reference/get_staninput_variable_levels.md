# Deprecated: get_staninput_variable_levels

**\[deprecated\]** `get_staninput_variable_levels()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`get_category_labels()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_labels.md),
[`get_group_labels()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_labels.md),
or
[`get_cue_labels()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_labels.md)
instead.

## Usage

``` r
get_staninput_variable_levels(
  x,
  variable = c("category", "group", "cue"),
  indices = NULL
)
```

## Arguments

- x:

  A model object.

- variable:

  Either "category", "group", or "cue".

- indices:

  Optional indices.

## See also

[`get_category_labels()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_labels.md),
[`get_group_labels()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_labels.md),
[`get_cue_labels()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_labels.md)
