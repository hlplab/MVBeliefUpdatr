# Deprecated: update_NIW_belief_by_one_observation

**\[deprecated\]** `update_NIW_belief_by_one_observation()` was
deprecated in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please
use
[`update_category_representation()`](https://hlplab.github.io/MVBeliefUpdatr/reference/update_category_representation.md)
or
[`update_template()`](https://hlplab.github.io/MVBeliefUpdatr/reference/update_template.md)
instead.

## Usage

``` r
update_NIW_belief_by_one_observation(
  prior_model,
  x_category,
  x,
  noise_treatment = "no_noise",
  lapse_treatment = "no_lapses",
  method = "label-certain",
  verbose = FALSE
)
```

## Arguments

- prior_model:

  An NIW category representation or ideal-adaptor model.

- x_category:

  Category label for the observation.

- x:

  Numeric observation vector.

- noise_treatment:

  Noise treatment.

- lapse_treatment:

  Lapse treatment.

- method:

  Updating method.

- verbose:

  Whether to print additional output.

## Value

An updated S7 representation or model.

## See also

[`update_category_representation()`](https://hlplab.github.io/MVBeliefUpdatr/reference/update_category_representation.md),
[`update_template()`](https://hlplab.github.io/MVBeliefUpdatr/reference/update_template.md)
