# Convert Between NIW/NIX Hyperparameters and Expected/Marginal Moments

Functions to convert between conjugate Normal-Inverse-Wishart (NIW) or
Normal-Inverse-Chi-Squared (NIX) hyperparameters (\\m, S\\) and the
expected or marginal distribution moments (\\\mu, \Sigma\\) (see Murphy
2012 p. 134) .

## Usage

``` r
get_expected_Sigma_from_S(S, nu)

get_S_from_expected_Sigma(Sigma, nu)

get_marginal_Sigma_from_S(S, nu, kappa)

get_expected_mu_from_m(m)

get_m_from_expected_mu(mu)
```

## Arguments

- S:

  Scatter matrix \\S\\ (numeric matrix, scalar, or list of matrices).

- nu:

  Strength of belief (pseudocount) about \\\Sigma\\ (\\\nu \> D + 1\\).

- Sigma:

  Expected category covariance matrix \\\Sigma\\ (numeric matrix,
  scalar, or list of matrices).

- kappa:

  Strength of belief (pseudocount) about \\\mu\\.

- m:

  Location hyperparameter \\m\\ (numeric vector, matrix, or list).

- mu:

  Expected category mean vector \\\mu\\ (numeric vector, matrix, or
  list).

## Value

For `get_expected_Sigma_from_S`, the expected covariance matrix \\\Sigma
= S / (\nu - D - 1)\\. For `get_S_from_expected_Sigma`, the scatter
matrix \\S = \Sigma \cdot (\nu - D - 1)\\. For
`get_marginal_Sigma_from_S`, the marginal category covariance matrix
\\\Sigma = (\kappa + 1) / \kappa \cdot S / (\nu - D - 1)\\. For
`get_expected_mu_from_m` and `get_m_from_expected_mu`, the mean vector
or location hyperparameter.

## Details

Under the NIW model with \\D\\ dimensions, scatter matrix \\S\\, and
degrees of freedom / pseudocounts \\\nu \> D + 1\\, the expected
covariance matrix is: \$\$\Sigma = \frac{S}{\nu - D - 1}\$\$ Conversely,
the scatter matrix is obtained from \\\Sigma\\ by: \$\$S = \Sigma \cdot
(\nu - D - 1)\$\$. For conjugate NIW and NIX models, the expected mean
\\\mu\\ is identically equal to the location parameter \\m\\: \$\$\mu =
m\$\$

## References

Murphy KP (2012). *Machine learning: a probabilistic perspective*. MIT
Press.

## See also

[`get_expected_sigma`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_parameters.md),
[`get_expected_mu`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_parameters.md)

## Examples

``` r
# 2D example
S <- matrix(c(4, 1, 1, 3), nrow = 2)
nu <- 10
kappa <- 5

# Convert scatter matrix to expected covariance matrix (D = 2)
Sigma <- get_expected_Sigma_from_S(S, nu = nu)
Sigma
#>           [,1]      [,2]
#> [1,] 0.5714286 0.1428571
#> [2,] 0.1428571 0.4285714

# Convert expected covariance matrix back to scatter matrix
get_S_from_expected_Sigma(Sigma, nu = nu)
#>      [,1] [,2]
#> [1,]    4    1
#> [2,]    1    3

# Marginal covariance matrix
get_marginal_Sigma_from_S(S, nu = nu, kappa = kappa)
#> [1] 0.6857143

# Mean and location conversions
m <- c(0, 1)
mu <- get_expected_mu_from_m(m)
get_m_from_expected_mu(mu)
#> [1] 0 1
```
