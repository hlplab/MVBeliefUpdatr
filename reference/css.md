# Get sum-of-squares matrix

Calculate centered (`css()`) or uncentered (`uss()`) sum-of-squares
matrix for a matrix of observations.

## Usage

``` r
css(x)

uss(x)
```

## Arguments

- x:

  A matrix of observations, with each row being a vector observation.

## Value

A square matrix.

## See also

[`css2cov`](https://hlplab.github.io/MVBeliefUpdatr/reference/uss2css.md),
[`cov2css`](https://hlplab.github.io/MVBeliefUpdatr/reference/uss2css.md),
[`uss2css`](https://hlplab.github.io/MVBeliefUpdatr/reference/uss2css.md)

## Examples

``` r
x <- matrix(rnorm(20), nrow = 10, ncol = 2)
# Centered sum-of-squares matrix
css(x)
#>           [,1]      [,2]
#> [1,] 5.7041628 0.5126964
#> [2,] 0.5126964 6.5593209
# Uncentered sum-of-squares matrix
uss(x)
#>           [,1]      [,2]
#> [1,] 5.8087603 0.8942042
#> [2,] 0.8942042 7.9508290
```
