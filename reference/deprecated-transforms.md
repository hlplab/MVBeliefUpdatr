# Deprecated: Transform and untransform cues by applying or undoing PCA, centering, and/or scaling

**\[deprecated\]** `transform_cues()` and `untransform_cues()` are
deprecated as of MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0.
Please use
[`get_affine_transform()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_affine_transform.md)
and
[`invert_transform()`](https://hlplab.github.io/MVBeliefUpdatr/reference/invert_transform.md)
instead.

## Usage

``` r
transform_cues(
  data,
  cues,
  center = if (pca) T else F,
  scale = F,
  pca = F,
  attach = T,
  transform.parameters = NULL,
  return.transformed.data = T,
  return.transform.parameters = F,
  return.transform.function = F,
  return.untransform.function = F
)

untransform_cues(
  data,
  cues,
  uncenter = NULL,
  unscale = NULL,
  unpca = NULL,
  attach = T,
  transform.parameters = NULL,
  return.untransformed.data = T,
  return.untransform.function = F
)
```

## Arguments

- data:

  `tibble` or `data.frame`.

- cues:

  Vector of characters with names of cue variables.

- center, uncenter:

  Should the data be (un)centered? (default: `TRUE` unless `pca = TRUE`)

- scale, unscale:

  Should the data be (un)standardized? (default: `TRUE` unless
  `pca = TRUE`)

- pca, unpca:

  Should the data be transformed into/back from orthogonal principal
  components?

- attach:

  Should the transformed cues be attached to `data`?

- transform.parameters:

  List of transforms (default: `NULL`).

- return.transformed.data, return.untransformed.data:

  Should (un)transformed data be returned?

- return.transform.parameters:

  Should the list of transforms be returned?

- return.transform.function, return.untransform.function:

  Should a transform function be returned?

## Value

By default a `data.frame`.
