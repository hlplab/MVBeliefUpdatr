# Deprecated: lift_likelihood_to_model

**\[deprecated\]** `lift_likelihood_to_model()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`new_mvg_ideal_observer()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md),
[`new_niw_ideal_adaptor()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md),
or
[`new_exemplar_model()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-exemplar.md)
with an S7 `category_template` instead.

**\[deprecated\]** `lift_exemplars_to_exemplar_model()` was deprecated
in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`new_exemplar_model()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-exemplar.md)
instead.

**\[deprecated\]** `lift_MVG_to_MVG_ideal_observer()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`new_mvg_ideal_observer()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md)
instead.

**\[deprecated\]** `lift_NIW_belief_to_NIW_ideal_adaptor()` was
deprecated in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please
use
[`new_niw_ideal_adaptor()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md)
instead.

**\[deprecated\]** `lift_MVG_ideal_observer_to_NIW_ideal_adaptor()` was
deprecated in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please
use
[`as_niw_ideal_adaptor()`](https://hlplab.github.io/MVBeliefUpdatr/reference/as_family_coercion.md)
instead.

## Usage

``` r
lift_likelihood_to_model(
  x,
  group = NULL,
  prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  verbose = F
)

lift_exemplars_to_exemplar_model(
  x,
  group = NULL,
  prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  verbose = F
)

lift_MVG_to_MVG_ideal_observer(
  x,
  group = NULL,
  prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  verbose = F
)

lift_NIW_belief_to_NIW_ideal_adaptor(
  x,
  group = NULL,
  prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  verbose = F
)

lift_MVG_ideal_observer_to_NIW_ideal_adaptor(
  x,
  group = NULL,
  kappa,
  nu,
  verbose = F
)
```

## Arguments

- x:

  An S7
  [`MVBU_CategoryRepresentationTemplate`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md)
  object.

- group:

  No longer supported; must be `NULL`.

- prior:

  Optional category-prior vector. (default: uniform over categories)

- lapse_rate:

  Optional lapse rate. (default: 0)

- lapse_bias:

  Optional lapse-bias vector. (default: same as `prior`)

- Sigma_noise:

  Optional perceptual-noise covariance matrix. (default: `NULL`)

- verbose:

  Ignored.

## Value

An S7 cognitive-model object of the family matching `x`.

## See also

[`new_mvg_ideal_observer()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md),
[`new_niw_ideal_adaptor()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md),
[`new_exemplar_model()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-exemplar.md)

[`new_exemplar_model()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-exemplar.md)

[`new_mvg_ideal_observer()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md)

[`new_niw_ideal_adaptor()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md)

[`as_niw_ideal_adaptor()`](https://hlplab.github.io/MVBeliefUpdatr/reference/as_family_coercion.md)
