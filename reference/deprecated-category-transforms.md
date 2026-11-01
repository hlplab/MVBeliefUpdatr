# Deprecated: Transform and untransform category means and covariances

**\[deprecated\]** `transform_category_mean()`,
`untransform_category_mean()`, `transform_category_cov()`, and
`untransform_category_cov()` are deprecated as of MVBeliefUpdatr 0.1.0
and will be removed in 0.2.0. Please use
[`transform_mean()`](https://hlplab.github.io/MVBeliefUpdatr/reference/transform_model.md),
[`untransform_mean()`](https://hlplab.github.io/MVBeliefUpdatr/reference/transform_model.md),
[`transform_cov()`](https://hlplab.github.io/MVBeliefUpdatr/reference/transform_model.md),
and
[`untransform_cov()`](https://hlplab.github.io/MVBeliefUpdatr/reference/transform_model.md)
instead.

## Usage

``` r
transform_category_mean(m, transform)

untransform_category_mean(m, transform)

transform_category_cov(S, transform)

untransform_category_cov(S, transform)
```

## Arguments

- m:

  A mean vector.

- transform:

  An
  [`MVBU_TransformInformation`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-TransformInformation-class.md)
  object.

- S:

  A covariance matrix.

## Value

A transformed or untransformed mean vector or covariance matrix.
