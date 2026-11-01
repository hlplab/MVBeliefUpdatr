# Get affine transformation

Computes an affine transformation for a set of `cues` in `data`. Returns
an
[`MVBU_TransformInformation`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-TransformInformation-class.md)
object representing a linear map \\f(x) = A(x + b)\\, where \\b\\ is the
shift vector and \\A\\ is the scale/rotation matrix.

## Usage

``` r
get_affine_transform(
  data,
  cues,
  type = c("identity", "center", "standardize", "PCA whiten", "ZCA whiten")[1]
)
```

## Arguments

- data:

  A `tibble` or `data.frame`.

- cues:

  A character vector of column names in `data` that should be
  transformed.

- type:

  The type of transformation to apply. Can be one of `"identity"`,
  `"center"`, `"standardize"`, `"PCA whiten"`, or `"ZCA whiten"`. Except
  for the identity transform, all transforms center. Standardize
  additionally divides by the standard deviations. PCA whitening rotates
  the data into the space of principal components with unit variance.
  ZCA whitening achieves unit variance while maintaining original
  orientation as close as possible. (default: `"identity"`)

## Value

An
[`MVBU_TransformInformation`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-TransformInformation-class.md)
object that can be used to transform data or models.

## Examples

``` r
df <- data.frame(
  vot = rnorm(50, 20, 10),
  f0 = rnorm(50, 220, 25)
)
# Standardize transformation
trans <- get_affine_transform(df, cues = c("vot", "f0"),
                              type = "standardize")
trans
#> <MVBU_TransformInformation>
```
