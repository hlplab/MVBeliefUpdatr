# Deprecated: is.ideal_adaptor_stanfit

**\[deprecated\]** `is.ideal_adaptor_stanfit()` is deprecated; use
`S7::S7_inherits(x, IdealAdaptorStanfit)` instead.

**\[deprecated\]** `ideal_adaptor_stanfit()` is deprecated; use the S7
constructor
[`MVBU_Stanfit()`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
or
[`fit_ideal_adaptor()`](https://hlplab.github.io/MVBeliefUpdatr/reference/fit_ideal_adaptor.md)
instead.

## Usage

``` r
is.ideal_adaptor_stanfit(x, verbose = FALSE)

ideal_adaptor_stanfit(
  data = data.frame(),
  staninput = NULL,
  stanvars = NULL,
  backend = "rstan",
  save_pars = NULL,
  stan_args = list(),
  stanfit = NULL,
  transform_information = NULL,
  criteria = list(),
  file = NULL,
  version = NULL,
  metadata = list()
)
```

## Arguments

- x:

  Object to be checked.

- verbose:

  Currently being ignored.

- data:

  A named `list` containing `$original` and `$transformed` data frames.
  Default: [`list()`](https://rdrr.io/r/base/list.html).

- staninput:

  A Stan input object containing the data handed to rstan (an
  [`IdealAdaptorStaninput`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Staninput.md)
  object). Default: `NULL`.

- stanvars:

  A `stanvars` object or `NULL`. Default: `NULL`.

- backend:

  The name of the backend used to fit the model (character). Default:
  `"rstan"`.

- save_pars:

  Optional storage for saved parameter names. Default: `NULL`.

- stan_args:

  Named list of additional control arguments that were passed to the
  Stan backend directly. Default:
  [`list()`](https://rdrr.io/r/base/list.html).

- stanfit:

  An object of class
  [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html)
  containing posterior draws. Default: `NULL`.

- transform_information:

  An object of type
  [`MVBU_TransformInformation`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-TransformInformation-class.md)
  or `NULL`. Default: `NULL`.

- criteria:

  An empty `list` for adding model fit criteria after estimation of the
  model. Default: [`list()`](https://rdrr.io/r/base/list.html).

- file:

  Optional name of a file in which the model object was stored in or
  loaded from. Default: `NULL`.

- version:

  The versions of MVBeliefUpdatr and rstan with which the model was
  fitted. Default: `NULL`.

- metadata:

  List containing auxiliary information including `label_information`.
  Default: [`list()`](https://rdrr.io/r/base/list.html).

## Value

A logical.

An object of class `MVBU_Stanfit`.
