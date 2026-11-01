# Convert atomic values or lists into arrays with inner/outer dimensions

This helper is used to shape prior information and other Stan inputs
into the dimensionality expected by the Stan programs.

## Usage

``` r
to_array(
  x,
  inner_dims = NULL,
  outer_dims = NULL,
  dimnames = NULL,
  simplify = TRUE
)
```

## Arguments

- x:

  Input value. Can be atomic, list, or NULL.

- inner_dims:

  Integer vector describing the inner dimensions of the array.

- outer_dims:

  Integer vector describing the outer dimensions of the array.

- dimnames:

  Optional dimnames for the resulting array.

- simplify:

  Logical controlling whether inner dimensions of length 1 are
  simplified.

## Value

An array.

## Examples

``` r
# Convert a vector into an array with outer dimension 3
to_array(c(1, 2), inner_dims = 2, outer_dims = 3)
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    1    2
#> [3,]    1    2

# Convert a list of matrices into an array
mat_list <- list(matrix(1:4, 2, 2), matrix(5:8, 2, 2))
to_array(mat_list, inner_dims = c(2, 2), outer_dims = 2)
#> , , 1
#> 
#>      [,1] [,2]
#> [1,]    1    3
#> [2,]    2    4
#> 
#> , , 2
#> 
#>      [,1] [,2]
#> [1,]    5    7
#> [2,]    6    8
#> 
```
