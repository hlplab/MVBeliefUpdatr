# Deprecated: get_NIW_posterior_predictive

**\[deprecated\]** Get posterior predictive density of observations x
given the Normal-Inverse-Wishart (NIW) parameters m, S, kappa, and nu.
This is the density of a multivariate Student-T distribution (see Murphy
2012 p. 134) .

## Usage

``` r
get_NIW_posterior_predictive(
  x,
  m,
  S,
  kappa,
  nu,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  log = TRUE
)
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

- Sigma_noise:

  Optional perceptual noise covariance matrix (default: `NULL`).

- noise_treatment:

  Noise treatment: `"no_noise"`, `"sample"`, or `"marginalize"`.

- log:

  Logical; if `TRUE`, return log-transformed density (default: `TRUE`).

## References

Murphy KP (2012). *Machine learning: a probabilistic perspective*. MIT
Press.
