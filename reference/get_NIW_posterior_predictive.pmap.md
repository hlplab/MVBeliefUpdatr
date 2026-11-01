# Deprecated: get_NIW_posterior_predictive.pmap

**\[deprecated\]** `get_NIW_posterior_predictive.pmap()` was deprecated
in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`get_NIW_posterior_predictive()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_NIW_posterior_predictive.md)
instead.

## Usage

``` r
get_NIW_posterior_predictive.pmap(x, m, S, kappa, nu, ...)
```

## Arguments

- x:

  Observation(s). Can be a vector with k elements, a matrix with k
  columns, a data frame, or a list of numeric vectors.

- m:

  Mean vector of length k.

- S:

  Scatter matrix of dimension k x k.

- kappa:

  Strength of belief on the mean (pseudocounts).

- nu:

  Strength of belief on the covariance matrix (pseudocounts).

- ...:

  Additional arguments.

## See also

[`get_NIW_posterior_predictive()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_NIW_posterior_predictive.md),
[`likelihood()`](https://hlplab.github.io/MVBeliefUpdatr/reference/likelihood.md)
