# Deprecated: make_MVG_ideal_observer_from_data

**\[deprecated\]** `make_MVG_ideal_observer_from_data()` was deprecated
in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`new_mvg_ideal_observer_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md)
instead.

## Usage

``` r
make_MVG_ideal_observer_from_data(
  data,
  group = NULL,
  category = "category",
  cues,
  prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  verbose = F
)
```

## Arguments

- data:

  The tibble or data.frame from which to construct the exemplar
  representations.

- group:

  No longer supported; must be `NULL`.

- category:

  Name of variable in `data` that contains the category information.
  (default: "category")

- cues:

  Name(s) of variables in `data` that contain the cue information.

- prior:

  Optional category-prior vector, forwarded as `category_prior`.
  (default: uniform over categories)

- lapse_rate:

  Optional lapse rate. (default: 0)

- lapse_bias:

  Optional lapse-bias vector. (default: same as `prior`)

- Sigma_noise:

  Optional perceptual-noise covariance matrix. (default: `NULL`)

- verbose:

  If true provides more information. (default: `FALSE`)

## Value

An `MVG_IdealObserver` object.

## See also

[`new_mvg_ideal_observer_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md)
