# Deprecated: get_posterior_from_MVG_ideal_observer

**\[deprecated\]** `get_posterior_from_MVG_ideal_observer()` was
deprecated in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please
use
[`posterior()`](https://hlplab.github.io/MVBeliefUpdatr/reference/posterior.md)
instead.

**\[deprecated\]** `get_posterior_from_model()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`posterior()`](https://hlplab.github.io/MVBeliefUpdatr/reference/posterior.md)
instead.

## Usage

``` r
get_posterior_from_MVG_ideal_observer(
  x,
  model,
  decision_rule = "sampling",
  noise_treatment = if (decision_rule == "sampling") "sample" else "no_noise",
  lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize"
)

get_posterior_from_model(model, ...)
```

## Arguments

- x:

  Observations.

- model:

  A model object.

- decision_rule:

  Decision rule.

- noise_treatment:

  Noise treatment.

- lapse_treatment:

  Lapse treatment.

- ...:

  Arguments passed to
  [`posterior`](https://hlplab.github.io/MVBeliefUpdatr/reference/posterior.md).

## Value

Posterior data frame.

## See also

[`posterior()`](https://hlplab.github.io/MVBeliefUpdatr/reference/posterior.md)
