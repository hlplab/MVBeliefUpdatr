# Invert an affine transformation

Computes the exact algebraic inverse of an affine transformation \\f(x)
= A(x + b)\\. The inverse transformation is itself an affine
transformation \\f^{-1}(y) = A^{-1}(y - A b) = A\_{inv}(y + b\_{inv})\\,
where \\A\_{inv} = A^{-1}\\ and \\b\_{inv} = -A b\\.

## Usage

``` r
invert_transform(transform, ...)

## S7 method for class <MVBU_TransformInformation>
invert_transform(transform, ...)
```

## Arguments

- transform:

  An
  [`MVBU_TransformInformation`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-TransformInformation-class.md)
  object.

- ...:

  Additional arguments passed to methods.

## Value

An inverted
[`MVBU_TransformInformation`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-TransformInformation-class.md)
object.

## Examples

``` r
trans <- MVBU_TransformInformation(
  cue.labels = c("c1", "c2"),
  type = "center",
  shift = c(1, -2),
  SCALE = diag(c(2, 0.5))
)
inv_trans <- invert_transform(trans)
inv_trans@shift
#> [1] -2  1
inv_trans@SCALE
#>      [,1] [,2]
#> [1,]  0.5    0
#> [2,]  0.0    2
```
