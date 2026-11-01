# S7 class for transformation information

Represents an invertible affine transformation \\f(x) = A(x + b)\\,
where \\b\\ is the shift vector and \\A\\ is the scale/rotation matrix.

## Usage

``` r
MVBU_TransformInformation(
  cue.labels = character(),
  type = "identity",
  shift = numeric(),
  SCALE = 1,
  transform.function = NULL
)
```

## Arguments

- cue.labels:

  Character vector of cue variable names.

- type:

  Transformation type (e.g., `"identity"`, `"center"`, `"standardize"`,
  `"PCA whiten"`, `"ZCA whiten"`).

- shift:

  Numeric shift vector added before scaling.

- SCALE:

  Numeric square invertible matrix or scalar scale factor.

- transform.function:

  Function to apply transformation to data.

## Examples

``` r
trans <- MVBU_TransformInformation(
  cue.labels = c("c1", "c2"),
  type = "center",
  shift = c(0, 0),
  SCALE = diag(2)
)
trans
#> <MVBU_TransformInformation>
```
