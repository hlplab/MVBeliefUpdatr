# Convert sum-of-square and covariance matrices

Convert centered sum-of-square matrices (css), uncentered sum-of-square
matrices (uss), and covariance matrices (cov) into each other.

## Usage

``` r
uss2css(uss, n, mean)

uss2cov(uss, n, mean)

css2uss(css, n, mean)

css2cov(css, n)

cov2css(cov, n)
```

## Arguments

- uss, css, cov:

  Uncentered sum-of-square, centered sum-of-square, and covariance
  matrix.

- n:

  Number of observations that have gone into the sum-of-square matrix.

- mean:

  Mean of the observations that have gone into the sum-of-square matrix.

## Value

A square matrix.

## See also

[`css`](https://hlplab.github.io/MVBeliefUpdatr/reference/css.md),
[`uss`](https://hlplab.github.io/MVBeliefUpdatr/reference/css.md),
[`cov2cor`](https://rdrr.io/r/stats/cor.html),
[`cor2cov`](https://hlplab.github.io/MVBeliefUpdatr/reference/cor2cov.md),
[`cov2tau`](https://hlplab.github.io/MVBeliefUpdatr/reference/cor2cov.md)

## Examples

``` r
x <- matrix(rnorm(20), nrow = 10, ncol = 2)
n <- nrow(x)
m <- colMeans(x)
u <- uss(x)
c_mat <- css(x)

# Convert uncentered sum-of-squares to centered sum-of-squares
uss2css(u, n = n, mean = m)
#>          [,1]     [,2]
#> [1,] 6.400554 1.185887
#> [2,] 1.185887 9.154892

# Convert centered sum-of-squares to covariance matrix
css2cov(c_mat, n = n)
#>           [,1]      [,2]
#> [1,] 0.7111726 0.1317652
#> [2,] 0.1317652 1.0172102

# Convert covariance matrix back to centered sum-of-squares
cov2css(cov(x), n = n)
#>          [,1]     [,2]
#> [1,] 6.400554 1.185887
#> [2,] 1.185887 9.154892

# Convert centered sum-of-squares back to uncentered
css2uss(c_mat, n = n, mean = m)
#>          [,1]     [,2]
#> [1,] 6.434553 1.102713
#> [2,] 1.102713 9.358366

# Convert uncentered sum-of-squares directly to covariance
uss2cov(u, n = n, mean = m)
#>           [,1]      [,2]
#> [1,] 0.7111726 0.1317652
#> [2,] 0.1317652 1.0172102
```
