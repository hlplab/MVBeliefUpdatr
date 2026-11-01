# Deprecated: get_MVG_likelihood

**\[deprecated\]** `get_MVG_likelihood()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`likelihood()`](https://hlplab.github.io/MVBeliefUpdatr/reference/likelihood.md)
instead.

**\[deprecated\]** `get_likelihood_from_MVG()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`likelihood()`](https://hlplab.github.io/MVBeliefUpdatr/reference/likelihood.md)
instead.

## Usage

``` r
get_MVG_likelihood(
  x,
  mu,
  Sigma,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  log = TRUE
)

get_likelihood_from_MVG(
  x,
  model,
  noise_treatment = "no_noise",
  log = TRUE,
  category = "category",
  category.label = NULL,
  wide = FALSE
)
```

## Arguments

- x:

  Observations.

- mu:

  Mean vector.

- Sigma:

  Covariance matrix.

- Sigma_noise:

  Optional noise covariance matrix.

- noise_treatment:

  Noise treatment.

- log:

  Logical; whether log likelihood is returned.

- model:

  An MVG model.

- category:

  Category column name.

- category.label:

  Category labels.

- wide:

  Logical; whether wide format is returned.

## Value

Numeric vector of likelihoods.

Likelihood data frame.

## See also

[`likelihood()`](https://hlplab.github.io/MVBeliefUpdatr/reference/likelihood.md)
