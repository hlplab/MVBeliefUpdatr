# Get covariance matrix from correlation matrix and standard deviations

Convert between a correlation matrix with standard deviations and a
covariance matrix.

## Usage

``` r
cor2cov(omega, tau)

cov2tau(v)
```

## Arguments

- omega:

  A correlation matrix.

- tau:

  A vector of standard deviations.

- v:

  A covariance matrix.

## Value

For `cor2cov`, a covariance matrix. For `cov2tau`, a numeric vector of
standard deviations.

## References

<https://stats.stackexchange.com/questions/62850/obtaining-covariance-matrix-from-correlation-matrix>

## See also

[stats::cov2cor](https://rdrr.io/r/stats/cor.html), `cor2cov`, `cov2tau`

## Examples

``` r
omega <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
tau <- c(2, 3)
# Correlation matrix + standard deviations to covariance matrix
Sigma <- cor2cov(omega, tau)
Sigma
#>      [,1] [,2]
#> [1,]    4    3
#> [2,]    3    9

# Extract standard deviations from covariance matrix
cov2tau(Sigma)
#> [1] 2 3
```
