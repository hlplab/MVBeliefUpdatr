# Coerce a Stanfit model object to a lightweight StanfitPosterior representation

Extracts posterior draws of category representation parameters and
creates a lightweight, self-contained
[MVBU_StanfitPosterior](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_StanfitPosterior.md)
object suitable for rapid categorization and visualization without
retaining the full Stan sampling engine.

## Usage

``` r
as_MVBU_stanfit_posterior(x, ...)

## S7 method for class <MVBU_Stanfit>
as_MVBU_stanfit_posterior(x, ...)

## S7 method for class <MVBU_StanfitPosterior>
as_MVBU_stanfit_posterior(x, ...)
```

## Arguments

- x:

  An
  [MVBU_Stanfit](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
  object.

- ...:

  Additional arguments (currently unused).

## Value

An
[MVBU_StanfitPosterior](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_StanfitPosterior.md)
object.

An
[MVBU_StanfitPosterior](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_StanfitPosterior.md)
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
post <- as_MVBU_stanfit_posterior(fit)
# }
```
