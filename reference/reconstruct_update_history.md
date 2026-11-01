# Reconstruct exposure update history as a model list

Reconstructs sequential cognitive model states or prior-to-posterior
update checkpoints into an
[MVBU_ModelList](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_ModelList.md).

## Usage

``` r
reconstruct_update_history(
  object,
  uncertainty_treatment = c("marginalize", "discard"),
  ndraws = 20L,
  step_size = 10L,
  groups = NULL,
  categories = NULL,
  parallel = FALSE,
  seed = 42L,
  ...
)

## S7 method for class <MVBU_Stanfit>
reconstruct_update_history(
  object,
  uncertainty_treatment = c("marginalize", "discard"),
  ndraws = 20L,
  step_size = 10L,
  groups = NULL,
  categories = NULL,
  parallel = FALSE,
  seed = 42L,
  ...
)
```

## Arguments

- object:

  An
  [MVBU_Stanfit](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
  object or cognitive model.

- uncertainty_treatment:

  Character string specifying treatment of posterior parameter
  uncertainty: `"marginalize"` (default) reconstructs trajectories
  across posterior draws and summarizes them, while `"discard"`
  reconstructs a single trajectory from the prior point estimate.

- ndraws:

  Number of posterior draws to use when
  `uncertainty_treatment = "marginalize"`. Defaults to `20L`.

- step_size:

  Observation step size between checkpoints.

- groups:

  Groups to reconstruct.

- categories:

  Categories to reconstruct.

- parallel:

  Logical; whether to parallelize draw reconstruction.

- seed:

  Random seed.

- ...:

  Additional arguments passed to methods.

## Value

An
[MVBU_ModelList](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_ModelList.md)
object.

## Examples

``` r
# \donttest{
fit <- example_ideal_adaptor_stanfit()
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
hist <- reconstruct_update_history(fit, ndraws = 2)
# }
```
