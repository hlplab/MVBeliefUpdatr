# Construct S7 Category Representations, Templates, and Models from Data

Dynamic constructors that inspect the input `type` (or `family`)
argument and dispatch to the appropriate family-specific constructor
from data.

## Usage

``` r
new_category_representation_from_data(
  data,
  type,
  category = "category",
  cues,
  ...
)

new_category_representation_template_from_data(
  data,
  type,
  category = "category",
  cues,
  ...
)

new_model_from_data(data, type, category = "category", cues, ...)
```

## Arguments

- data:

  A data frame containing observation cues and category labels. For
  individual category representations, data must contain observations
  from exactly one category.

- type:

  Family type string: `"UVG"`, `"NIX"`, `"MUVG"`, `"MNIX"`, `"MVG"`,
  `"NIW"`, or `"EXEMPLAR"`.

- category:

  Character string giving the category column name in `data`. Defaults
  to `"category"`.

- cues:

  Character vector of cue column names in `data`.

- ...:

  Additional model and family-specific parameters passed to
  constructors. For cognitive models, this includes `category_prior`
  (prior category probabilities), `decision_rule` (`"sampling"`,
  `"criterion"`, or `"proportional"`), `lapse_rate` (scalar in
  `[0, 1]`), `lapse_bias`, `Sigma_noise` (perceptual noise covariance),
  `noise_treatment` (`"no_noise"`, `"sample"`, or `"marginalize"`), and
  `lapse_treatment` (`"no_lapses"`, `"sample"`, or `"marginalize"`). For
  specific representation families, this includes parameters such as
  `kappa_0`, `nu_0`, `m_0`, `S_0`, `c`, or `bandwidth`.

## Value

An S7 representation, template, or cognitive model object.

## See also

[family-uvg](https://hlplab.github.io/MVBeliefUpdatr/reference/family-uvg.md),
[family-nix](https://hlplab.github.io/MVBeliefUpdatr/reference/family-nix.md),
[family-muvg](https://hlplab.github.io/MVBeliefUpdatr/reference/family-muvg.md),
[family-mnix](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mnix.md),
[family-mvg](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md),
[family-niw](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md),
[family-exemplar](https://hlplab.github.io/MVBeliefUpdatr/reference/family-exemplar.md)

## Examples

``` r
df <- data.frame(
  category = c("A", "A", "B", "B"),
  x = c(1, 2, 5, 6),
  y = c(2, 1, 6, 5)
)

# Create an MVG category representation for category "A"
cat_a <- new_category_representation_from_data(
  df[df$category == "A", ],
  type = "MVG",
  category = "category",
  cues = c("x", "y")
)
cat_a
#> <MVBeliefUpdatr::MVG_CategoryRepresentation>
#>   Category: A
#>   Cues (2): x, y
#>   mu: c(1.5, 1.5)
#>   Sigma:
#>      x    y
#> x  0.5 -0.5
#> y -0.5  0.5

# Create a full category template from data
tpl <- new_category_representation_template_from_data(
  df,
  type = "MVG",
  category = "category",
  cues = c("x", "y")
)
tpl
#> <MVBU_CategoryRepresentationTemplate>
#>   Categories (2): A, B
#>   Cues (2): x, y
#>   Representations (2):
#>     A: MVG(mu = vector(2), Sigma = matrix(2, 2))
#>     B: MVG(mu = vector(2), Sigma = matrix(2, 2))

# Create a full cognitive model from data
mod <- new_model_from_data(
  df,
  type = "MVG",
  category = "category",
  cues = c("x", "y")
)
mod
#> <MVBeliefUpdatr::MVG_IdealObserver>
#>   Categories (2): A, B
#>   Cues (2): x, y
#>   Decision rule: sampling
#>   Category prior: A: 0.5, B: 0.5
#>   Lapse rate: 0 (treatment: no_lapses)
#>   Perceptual noise: none (treatment: no_noise)
#>   Category representations (2):
#>     A: MVG(mu = vector(2), Sigma = matrix(2, 2))
#>     B: MVG(mu = vector(2), Sigma = matrix(2, 2))
```
