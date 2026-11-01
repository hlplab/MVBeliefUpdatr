# Deprecated: get_NIX_posterior_predictive

**\[deprecated\]** Get posterior predictive density of 1D observations x
given the Normal-Inverse-Chi-Squared (NIX) parameters m, sigma2 (or S),
kappa, and nu. In 1D, this evaluates a Student-T density.

## Usage

``` r
get_NIX_posterior_predictive(
  x,
  m,
  sigma2,
  kappa,
  nu,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  log = TRUE
)
```

## Arguments

- x:

  Observation(s) (numeric vector, matrix with 1 column, or 1D data
  frame).

- m:

  Prior mean location.

- sigma2:

  Prior expected variance parameter (or variance scale \\\sigma_0^2\\).

- kappa:

  Strength of belief on the mean (pseudocounts).

- nu:

  Degrees of freedom / strength of belief on variance (pseudocounts).

- Sigma_noise:

  Optional 1x1 noise covariance / variance (default: `NULL`).

- noise_treatment:

  Noise treatment: `"no_noise"`, `"sample"`, or `"marginalize"`.

- log:

  Logical; if `TRUE`, return log-transformed density (default: `TRUE`).
