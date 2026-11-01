# Deprecated: get_categorization_from_MVG_ideal_observer

**\[deprecated\]** `get_categorization_from_MVG_ideal_observer()` was
deprecated in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please
use
[`categorize()`](https://hlplab.github.io/MVBeliefUpdatr/reference/categorize.md)
instead.

**\[deprecated\]** `get_categorization_from_NIW_ideal_adaptor()` was
deprecated in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please
use
[`categorize()`](https://hlplab.github.io/MVBeliefUpdatr/reference/categorize.md)
instead.

**\[deprecated\]** `get_categorization_from_exemplar_model()` was
deprecated in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please
use
[`categorize()`](https://hlplab.github.io/MVBeliefUpdatr/reference/categorize.md)
instead.

**\[deprecated\]** `get_categorization_from_model()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`categorize()`](https://hlplab.github.io/MVBeliefUpdatr/reference/categorize.md)
instead.

## Usage

``` r
get_categorization_from_MVG_ideal_observer(
  x,
  model,
  decision_rule = "sampling",
  noise_treatment = if (decision_rule == "sampling") "sample" else "no_noise",
  lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  simplify = FALSE
)

get_categorization_from_NIW_ideal_adaptor(
  x,
  model,
  decision_rule = "sampling",
  noise_treatment = if (decision_rule == "sampling") "sample" else "no_noise",
  lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  simplify = FALSE,
  verbose = FALSE
)

get_categorization_from_exemplar_model(
  x,
  model,
  decision_rule,
  noise_treatment = if (decision_rule == "sampling") "sample" else "no_noise",
  lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  simplify = FALSE
)

get_categorization_from_model(model, decision_rule = "sampling", ...)
```

## Arguments

- x:

  Observations.

- model:

  A model object.

- decision_rule:

  Decision rule to use.

- noise_treatment:

  Noise treatment.

- lapse_treatment:

  Lapse treatment.

- simplify:

  Logical; whether to simplify output.

- verbose:

  Logical; verbosity.

- ...:

  Arguments passed to
  [`categorize`](https://hlplab.github.io/MVBeliefUpdatr/reference/categorize.md).

## Value

Categorization data frame or vector.

Categorization data frame or vector.

Categorization data frame or vector.

## See also

[`categorize()`](https://hlplab.github.io/MVBeliefUpdatr/reference/categorize.md)
