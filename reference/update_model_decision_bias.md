# Update model's decision biases based on labelled input

Returns model with updated decision biases. This updating process is
order-sensitive. For each input observation x, a posterior expectation
for each category is obtained based on the specified noise and lapse
treatment. This expectation is used to calculate the surprisal (in bits)
of the category label x_category. The surprisal is weighted by the
learning rate beta to increase the (log-odds of the) decision bias for
that category.

## Usage

``` r
update_model_decision_bias_by_one_observation(
  model,
  beta,
  x,
  x_category,
  noise_treatment = .infer_noise_treatment(dplyr::first(model$Sigma_noise)),
  lapse_treatment = "no_lapses",
  update_prior = T
)

update_model_decision_bias_incrementally(
  model,
  beta,
  exposure,
  exposure.category = "category",
  exposure.cues = get_cue_labels_from_model(model),
  exposure.order = NULL,
  noise_treatment = .infer_noise_treatment(dplyr::first(model$Sigma_noise)),
  lapse_treatment = "no_lapses",
  keep.update_history = TRUE,
  keep.exposure_data = FALSE,
  verbose = FALSE
)
```

## Arguments

- model:

  A model object with decision biases.

- beta:

  The learning rate with which decision biases change. Might be
  integrated into model objects in the future.

- x:

  The cues of single observation.

- x_category:

  The category label(s) of one or more observations.

- noise_treatment:

  Determines whether perceptual noise is considered during
  categorization, and how. Can be `"no_noise"`, `"sample"`, or
  `"marginalize"`.

- lapse_treatment:

  Determines whether attentional lapses can occur during which no
  updating occurs. Can be `"no_lapses"`, `"sample"`, or `"marginalize"`.

- update_prior:

  Should the prior probability of each category be updated along with
  the decision bias? (default: `TRUE`)

- exposure:

  `data.frame` or `tibble` with exposure data. Each row is assumed to
  contain one observation.

- exposure.category:

  Name of variable in `data` that contains the category information.
  (default: `"category"`)

- exposure.cues:

  Name(s) of variables in `data` that contain the cue information. By
  default these cue names are extracted from the prior object.

- exposure.order:

  Name of variable in `data` that contains the order of the exposure
  data. If `NULL` the exposure data is assumed to be in the order in
  which it should be presented.

- keep.update_history:

  Should the history of the updating be stored and returned? (default:
  `TRUE`)

- keep.exposure_data:

  Should the input data be included in the output? (default: `FALSE`)

- verbose:

  Should more informative output be provided? (default: `FALSE`)

## Value

A model object with updated decision biases.

## Examples

``` r
if (FALSE) { # \dontrun{
# Updating decision biases for a cognitive model
updated_model <- update_model_decision_bias_by_one_observation(
  model = model,
  beta = 0.1,
  x = c(0, 0),
  x_category = "b"
)
} # }
```
