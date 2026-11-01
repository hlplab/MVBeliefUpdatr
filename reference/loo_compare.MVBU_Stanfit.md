# Compare Models Fitted with Stan

Compare fitted
[`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
objects using Leave-One-Out Cross-Validation (`loo`) or WAIC via
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html).

## Usage

``` r
# S3 method for class 'MVBU_Stanfit'
loo_compare(x, ..., criterion = "loo")
```

## Arguments

- x:

  An
  [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
  object.

- ...:

  Additional
  [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
  objects or arguments passed to
  [`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html).

- criterion:

  Character string specifying the criterion to use for model comparison
  (`"loo"` or `"waic"`). Default: `"loo"`.

## Value

A matrix of class `compare.loo` (see
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html))
containing model comparison metrics.

## See also

[`add_criterion`](https://hlplab.github.io/MVBeliefUpdatr/reference/add_criterion.md),
[`loo.MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/loo.MVBU_Stanfit.md)

## Examples

``` r
# \donttest{
fit1 <- example_ideal_adaptor_stanfit()
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
fit2 <- example_ideal_adaptor_stanfit(model_family = "NIX")
#> Warning: There were 1 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#> https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 2.1, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
loo::loo_compare(fit1, fit2, criterion = "loo")
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
#>   model elpd_diff se_diff p_worse diag_diff        diag_elpd
#>  model1       0.0     0.0      NA           15 k_psis > 0.28
#>  model2      -1.5     0.9    0.96   N < 100 15 k_psis > 0.28
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
# }
```
