# Deprecated: get_NIW_categorization_function

**\[deprecated\]** `get_NIW_categorization_function()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`get_category_posterior_function()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_category_posterior_function.md)
instead.

**\[deprecated\]**
`get_categorization_function_from_NIW_ideal_adaptor()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`get_category_posterior_function()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_category_posterior_function.md)
instead.

## Usage

``` r
get_NIW_categorization_function(
  ms,
  Ss,
  kappas,
  nus,
  priors = rep(1/length(ms), length(ms)),
  lapse_rate = 0,
  lapse_biases = rep(1/length(ms), length(ms)),
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = if (lapse_rate > 0) "marginalize" else "no_lapses"
)

get_categorization_function_from_NIW_ideal_adaptor(model, ...)
```

## Arguments

- ms:

  Means of the multivariate normal distributions over category means.

- Ss:

  Scatter matrices of the inverse Wishart distribution over category
  covariance matrices.

- kappas:

  Strength of the beliefs into the distribution over category means.

- nus:

  Strength of the beliefs into the distribution over category covariance
  matrices.

- priors:

  Vector of categories' prior probabilities.

- lapse_rate:

  A lapse rate for the categorization responses.

- lapse_biases:

  A lapse bias for the categorization responses.

- Sigma_noise:

  A noise matrix.

- noise_treatment:

  Noise treatment.

- lapse_treatment:

  Lapse treatment.

- model:

  A model object.

- ...:

  Additional arguments.

## Value

A categorization function.

A categorization function.

## See also

[`get_category_posterior_function()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_category_posterior_function.md)
