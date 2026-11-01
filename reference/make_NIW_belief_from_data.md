# Deprecated: make_NIW_belief_from_data

**\[deprecated\]** `make_NIW_belief_from_data()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`new_niw_category_representation_template_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md)
instead.

## Usage

``` r
make_NIW_belief_from_data(
  data,
  group = NULL,
  category = "category",
  cues,
  kappa = nu,
  nu = length(cues) + 2,
  verbose = F
)

make_NIW_prior_from_data(
  data,
  group = NULL,
  category = "category",
  cues,
  kappa = nu,
  nu = length(cues) + 2,
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

- kappa:

  Strength of belief (pseudocount) about the category mean. (default:
  same as `nu`)

- nu:

  Strength of belief (pseudocount) about the category covariance matrix.
  (default: number of cues + 2)

- verbose:

  If true provides more information. (default: `FALSE`)

## Value

An `MVBU_CategoryRepresentationTemplate` of NIW category
representations.

## See also

[`new_niw_category_representation_template_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md)
