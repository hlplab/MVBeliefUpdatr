# Update a model's category-representation template from observations

Currently implemented for `NIW_IdealAdaptor`; methods for other model
types will be added as their update workflows are migrated to S7.

`update_template()` updates the category template representations of an
ideal adaptor model (`NIW_IdealAdaptor`, `NIX_IdealAdaptor`, or
`MNIX_IdealAdaptor`) given new observations.

## Usage

``` r
update_template(
  x,
  observations,
  updating = c("batch", "incremental"),
  keep_history = FALSE,
  lapse_treatment = "no_lapses",
  noise_treatment = "no_noise",
  update_method = "label-certain",
  ...
)

## S7 method for classes <MVBeliefUpdatr::NIW_IdealAdaptor>, <any>
update_template(
  x,
  observations,
  updating = c("batch", "incremental"),
  keep_history = FALSE,
  lapse_treatment = "no_lapses",
  noise_treatment = "no_noise",
  update_method = "label-certain",
  ...
)

## S7 method for classes <MVBeliefUpdatr::NIX_IdealAdaptor>, <any>
update_template(
  x,
  observations,
  updating = c("batch", "incremental"),
  keep_history = FALSE,
  lapse_treatment = "no_lapses",
  noise_treatment = "no_noise",
  update_method = "label-certain",
  ...
)

## S7 method for classes <MVBeliefUpdatr::MNIX_IdealAdaptor>, <any>
update_template(
  x,
  observations,
  updating = c("batch", "incremental"),
  keep_history = FALSE,
  lapse_treatment = "no_lapses",
  noise_treatment = "no_noise",
  update_method = "label-certain",
  ...
)
```

## Arguments

- x:

  An `NIW_IdealAdaptor`, `NIX_IdealAdaptor`, or `MNIX_IdealAdaptor`
  object.

- observations:

  A data.frame or tibble containing cue columns and, for
  `update_method = "label-certain"`, a `category` column.

- updating:

  Either `"batch"` or `"incremental"`. Batch updates are currently
  supported only for `update_method = "label-certain"`. Default is
  `"batch"`.

- keep_history:

  If `TRUE` with incremental updating, return a list of all intermediate
  models followed by the final model. Default is `FALSE`.

- lapse_treatment:

  Treatment of attentional lapses during updating for label-uncertain
  models: `"no_lapses"`, `"sample"`, or `"marginalize"`.

  - `"no_lapses"`: Lapses are ignored; updating reflects pure normative
    evidence evaluation (even if `lapse_rate > 0`).

  - `"sample"`: Simulates stochastic attentional lapses on individual
    observations. With probability `lapse_rate`, a lapse occurs and the
    observation is treated as uninformative during updating.

  - `"marginalize"`: Evaluates updating in expectation across infinitely
    many encounters; evidence is downweighted by \\1 - \lambda\\.

  Default is `"no_lapses"`.

- noise_treatment:

  Treatment of perceptual noise during categorization for
  label-uncertain updating: `"no_noise"`, `"sample"`, or
  `"marginalize"`.

  - `"no_noise"`: No perceptual noise is applied to input cues, and no
    noise is assumed during categorization.

  - `"sample"`: Simulates sensory noise on a specific observation. Noise
    drawn from \\\mathcal{N}(0, \Sigma\_{\text{noise}})\\ is added to
    input cues (and projected into composite cue space in early
    integration), and the observer accounts for their own noise
    uncertainty during categorization by adding
    \\\Sigma\_{\text{noise}}\\ to their likelihood scale.

  - `"marginalize"`: Evaluates categorization in the limit across
    infinitely many noisy encounters (expected consequences). No noise
    is sampled onto the input stimulus, but \\\Sigma\_{\text{noise}}\\
    is added to the category covariance scale.

  Default is `"no_noise"`.

- update_method:

  One of `"no-updating"`, `"label-certain"`, `"nolabel-criterion"`,
  `"nolabel-sampling"`, `"nolabel-proportional"`, or
  `"nolabel-uniform"`. Default is `"label-certain"`.

- ...:

  Additional arguments passed to methods.

## Value

An updated model, or a history of updated models when requested by a
method.

The updated model, or a list of intermediate models when
`updating = "incremental"` and `keep_history = TRUE`.

## See also

[`update_category_representation`](https://hlplab.github.io/MVBeliefUpdatr/reference/update_category_representation.md)

## Examples

``` r
mod <- example_niw_ideal_adaptor(n_cues = 2)
obs <- tibble::tibble(
  category = factor(c("/b/", "/d/")),
  VOT = c(10, 15),
  f0_semitones = c(-1, -0.5)
)
updated <- update_template(mod, obs)
```
