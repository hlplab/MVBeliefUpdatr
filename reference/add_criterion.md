# Add Information Criteria to a Stanfit Object (brms-compatible)

Computes and stores model fit criteria (such as `"loo"`, `"waic"`,
`"kfold"`, `"loo_subsample"`, `"bayes_R2"`, `"loo_R2"`, or `"marglik"`)
in an
[`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
object's `@criteria` slot, following the interface convention of brms.

## Usage

``` r
add_criterion(x, criterion = "loo", model_name = NULL, overwrite = TRUE, ...)
```

## Arguments

- x:

  An
  [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
  object.

- criterion:

  Character vector of criteria to add. Supported: `"loo"`, `"waic"`,
  `"kfold"`, `"loo_subsample"`, `"bayes_R2"`, `"loo_R2"`, `"marglik"`.
  Default: `"loo"`.

- model_name:

  Optional character string giving a name for the model.

- overwrite:

  Logical; whether to overwrite existing criteria. Default: `TRUE`.

- ...:

  Additional arguments passed to underlying criterion evaluation
  functions (such as
  [`loo`](https://mc-stan.org/loo/reference/loo.html),
  [`waic`](https://mc-stan.org/loo/reference/waic.html),
  [`kfold`](https://mc-stan.org/loo/reference/kfold-generic.html), or
  [`bridgesampling::bridge_sampler`](https://rdrr.io/pkg/bridgesampling/man/bridge_sampler.html)).

## Value

The updated
[`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
object with criteria stored in `x@criteria`.

## See also

[`loo_compare.MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/loo_compare.MVBU_Stanfit.md),
[`loo.MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/loo.MVBU_Stanfit.md)

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
fit <- add_criterion(fit, criterion = "loo")
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
fit@criteria$loo
#> 
#> Computed from 25 by 15 log-likelihood matrix.
#> 
#>          Estimate  SE
#> elpd_loo    -10.4 2.7
#> p_loo         1.7 0.7
#> looic        20.8 5.5
#> ------
#> MCSE of elpd_loo is NA.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.8, 1.3]).
#> 
#> Pareto k diagnostic values:
#>                           Count Pct.    Min. ESS
#> (-Inf, 0.28]   (good)      0      0.0%  <NA>    
#>    (0.28, 1]   (bad)       0      0.0%  <NA>    
#>     (1, Inf)   (very bad) 15    100.0%  <NA>    
#> See help('pareto-k-diagnostic') for details.
# }
```
