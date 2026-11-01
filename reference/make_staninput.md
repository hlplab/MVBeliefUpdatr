# Deprecated: make_staninput

**\[deprecated\]** `make_staninput()` was deprecated in MVBeliefUpdatr
0.1.0 and will be removed in 0.2.0. Please use
[`new_ideal_adaptor_stanfit_input()`](https://hlplab.github.io/MVBeliefUpdatr/reference/new_ideal_adaptor_stanfit_input.md)
instead.

**\[deprecated\]** `make_ideal_adaptor_stanfit_input()` was deprecated
in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`new_ideal_adaptor_stanfit_input()`](https://hlplab.github.io/MVBeliefUpdatr/reference/new_ideal_adaptor_stanfit_input.md)
instead.

## Usage

``` r
make_staninput(
  exposure,
  test,
  cues,
  category = "category",
  response = "response",
  group = "group",
  group.unique = NULL,
  lapse_rate = NULL,
  mu_0 = NULL,
  Sigma_0 = NULL,
  tau_scale = NULL,
  L_omega_eta = NULL,
  control = control_staninput(),
  stanmodel = "NIW_ideal_adaptor",
  verbose = FALSE
)

make_ideal_adaptor_stanfit_input(
  exposure,
  test,
  cues,
  category = "category",
  response = "response",
  group = "group",
  group.unique = NULL,
  lapse_rate = NULL,
  mu_0 = NULL,
  Sigma_0 = NULL,
  tau_scale = NULL,
  L_omega_eta = NULL,
  control = control_staninput(),
  stanmodel = "NIW_ideal_adaptor",
  verbose = FALSE
)
```

## Arguments

- exposure:

  A data frame or tibble with exposure data. Each row is assumed to
  contain one observation.

- test:

  A data frame or tibble with test data. Each row is assumed to contain
  one observation.

- cues:

  Names of columns with cue values. These columns must be present in
  both `exposure` and `test`.

- category:

  Name of the column in `exposure` that stores the category information.
  (default: "category")

- response:

  Name of the column in `test` that stores the response information.
  (default: "response_category")

- group:

  Name of the column in `exposure` and `test` that stores the grouping
  information. (default: "group")

- tau_scale:

  **\[deprecated\]** Deprecated. Use
  `priors = parameter_prior(..., class = "tau_0")` in
  [`fit_ideal_adaptor`](https://hlplab.github.io/MVBeliefUpdatr/reference/fit_ideal_adaptor.md)
  instead.

- L_omega_eta:

  **\[deprecated\]** Deprecated. Use
  `priors = parameter_prior(..., class = "cor_0")` in
  [`fit_ideal_adaptor`](https://hlplab.github.io/MVBeliefUpdatr/reference/fit_ideal_adaptor.md)
  instead.

- control:

  A list of control parameters for the constructor, such as
  `list(transform_type = "standardize")` or produced by
  [`control_staninput()`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_StanfitInput.md).

## Value

An object of class
[`IdealAdaptorStanfitInput`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_StanfitInput.md).

An object of class
[`IdealAdaptorStanfitInput`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_StanfitInput.md).

## See also

[`new_ideal_adaptor_stanfit_input()`](https://hlplab.github.io/MVBeliefUpdatr/reference/new_ideal_adaptor_stanfit_input.md)
