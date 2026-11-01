# MVBeliefUpdatr Stanfit Classes and Constructors

An S7 base class for Stan fit objects.

## Usage

``` r
MVBU_Stanfit(
  data = list(),
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
  metadata = list(),
  cache = list()
)

IdealAdaptorStanfit(
  data = list(),
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
  metadata = list(),
  cache = list()
)

MNIX_IdealAdaptorStanfit(
  data = list(),
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
  metadata = list(),
  cache = list()
)

IdealIntegrator_MNIX_IdealAdaptorStanfit(
  data = list(),
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
  metadata = list(),
  cache = list()
)

WeightedIntegrator_MNIX_IdealAdaptorStanfit(
  data = list(),
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
  metadata = list(),
  cache = list()
)

CategoryWeightedIntegrator_MNIX_IdealAdaptorStanfit(
  data = list(),
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
  metadata = list(),
  cache = list()
)

GroupWeightedIntegrator_MNIX_IdealAdaptorStanfit(
  data = list(),
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
  metadata = list(),
  cache = list()
)

GroupCategoryWeightedIntegrator_MNIX_IdealAdaptorStanfit(
  data = list(),
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
  metadata = list(),
  cache = list()
)

InformativityWeightedIntegrator_MNIX_IdealAdaptorStanfit(
  data = list(),
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
  metadata = list(),
  cache = list()
)

MNIX_IdealAdaptor_IdealIntegratorStanfit(
  data = list(),
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
  metadata = list(),
  cache = list()
)

MNIX_IdealAdaptor_WeightedIntegratorStanfit(
  data = list(),
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
  metadata = list(),
  cache = list()
)

MNIX_IdealAdaptor_CategoryWeightedIntegratorStanfit(
  data = list(),
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
  metadata = list(),
  cache = list()
)

MNIX_IdealAdaptor_GroupWeightedIntegratorStanfit(
  data = list(),
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
  metadata = list(),
  cache = list()
)

MNIX_IdealAdaptor_GroupCategoryWeightedIntegratorStanfit(
  data = list(),
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
  metadata = list(),
  cache = list()
)

MNIX_IdealAdaptor_InformativityWeightedIntegratorStanfit(
  data = list(),
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
  metadata = list(),
  cache = list()
)
```

## Arguments

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

- cache:

  Internal cache list storing pre-extracted `$draws` and `$summary`
  objects populated by
  [`add_posterior_latents`](https://hlplab.github.io/MVBeliefUpdatr/reference/add_posterior_latents.md).
  Default: [`list()`](https://rdrr.io/r/base/list.html).

## See also

[`add_posterior_latents`](https://hlplab.github.io/MVBeliefUpdatr/reference/add_posterior_latents.md),
[`fit_ideal_adaptor`](https://hlplab.github.io/MVBeliefUpdatr/reference/fit_ideal_adaptor.md),
[`vignette("fitting-and-working-with-stanfit-models")`](https://hlplab.github.io/MVBeliefUpdatr/articles/fitting-and-working-with-stanfit-models.md)

## Examples

``` r
# Instantiate an empty Stanfit container object
fit <- MVBU_Stanfit()
is(fit, "MVBU_Stanfit")
#> [1] TRUE

# \donttest{
# Load a pre-compiled example Stanfit object
ex_fit <- example_ideal_adaptor_stanfit()
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
class(ex_fit)
#> [1] "NIW_IdealAdaptorStanfit" "IdealAdaptorStanfit"    
#> [3] "MVBU_Stanfit"            "MVBU_Object"            
#> [5] "S7_object"              
# }
```
