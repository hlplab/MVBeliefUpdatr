# Deprecated: make_MVG_from_data

**\[deprecated\]** `make_MVG_from_data()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`new_mvg_category_representation_template_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md)
instead.

## Usage

``` r
make_MVG_from_data(
  data,
  group = NULL,
  category = "category",
  cues,
  verbose = F
)
```

## Arguments

- data:

  The tibble or data.frame from which to construct the exemplar
  representations.

- group:

  No longer supported; must be `NULL`.

- category:

  Name of variable in `data` that contains the category information.
  (default: "category")

- cues:

  Name(s) of variables in `data` that contain the cue information.

- verbose:

  If true provides more information. (default: `FALSE`)

## Value

An `MVBU_CategoryRepresentationTemplate` of MVG category
representations.

## See also

[`new_mvg_category_representation_template_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md)
