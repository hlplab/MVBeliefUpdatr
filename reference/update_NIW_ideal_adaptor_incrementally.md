# Deprecated: update_NIW_ideal_adaptor_incrementally

**\[deprecated\]** `update_NIW_ideal_adaptor_incrementally()` was
deprecated in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please
use
[`update_template()`](https://hlplab.github.io/MVBeliefUpdatr/reference/update_template.md)
instead.

**\[deprecated\]** `update_NIW_ideal_adaptor_batch()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`update_template()`](https://hlplab.github.io/MVBeliefUpdatr/reference/update_template.md)
instead.

**\[deprecated\]** `update_NIW_beliefs_incrementally()` was deprecated
in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`update_template()`](https://hlplab.github.io/MVBeliefUpdatr/reference/update_template.md)
instead.

## Usage

``` r
update_NIW_ideal_adaptor_incrementally(
  prior_model,
  exposure,
  exposure.category = "category",
  exposure.cues = get_cue_labels(prior_model),
  exposure.order = NULL,
  noise_treatment = "no_noise",
  lapse_treatment = "no_lapses",
  method = "label-certain",
  keep.update_history = TRUE,
  keep.exposure_data = FALSE,
  verbose = FALSE
)

update_NIW_ideal_adaptor_batch(
  prior_model,
  exposure,
  exposure.category = "category",
  exposure.cues = get_cue_labels(prior_model),
  noise_treatment = "no_noise",
  verbose = FALSE
)

update_NIW_beliefs_incrementally(prior_model, exposure, ...)
```

## Arguments

- prior_model:

  An NIW S7 ideal-adaptor model.

- exposure:

  Observation data.

- exposure.category:

  Name of the category column.

- exposure.cues:

  Cue columns.

- exposure.order:

  Deprecated ordering column.

- noise_treatment:

  Noise treatment.

- lapse_treatment:

  Lapse treatment.

- method:

  Updating method.

- keep.update_history:

  Whether to return intermediate models.

- keep.exposure_data:

  Deprecated and ignored.

- verbose:

  Whether to print additional output.

## Value

An updated S7 model or model history.

## See also

[`update_template()`](https://hlplab.github.io/MVBeliefUpdatr/reference/update_template.md)
