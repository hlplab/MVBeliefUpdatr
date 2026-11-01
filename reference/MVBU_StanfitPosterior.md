# S7 class for Stanfit posterior parameter draws

Represents a posterior distribution over a cognitive model, as obtained
for example from an
[MVBU_Stanfit](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
object. Holds MCMC parameter draw distributions without carrying
sampling chains, diagnostics, or raw model source structures, and is
thus more lightweight than the
[MVBU_Stanfit](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
class.

## Usage

``` r
MVBU_StanfitPosterior(
  draws = list(),
  metadata = list(),
  data = list(),
  cache = list(),
  transform_information = NULL
)
```

## Arguments

- draws:

  A list containing `raw` and `summary` posterior draws. Default:
  [`list()`](https://rdrr.io/r/base/list.html).

- metadata:

  A list containing auxiliary information including `label_information`
  (cues, categories, group), `model_family`, `parameters`,
  `original_variable_names`, and `exposure_category_statistics`.
  Default: [`list()`](https://rdrr.io/r/base/list.html).

- data:

  A named `list` containing `$original` and `$transformed` data frames.
  Default: [`list()`](https://rdrr.io/r/base/list.html).

- cache:

  A list storing cached computations and summaries. Default:
  [`list()`](https://rdrr.io/r/base/list.html).

- transform_information:

  An optional
  [MVBU_TransformInformation](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-TransformInformation-class.md)
  object. Default: `NULL`.

## Examples

``` r
# Construct a lightweight StanfitPosterior object
post <- MVBU_StanfitPosterior()
post
#> MVBeliefUpdatr Lightweight Stanfit Posterior
#> Model Family     : MVBU_StanfitPosterior 
#> Group Context    :  
#> Cues             :  
#> Categories       :  

# \donttest{
# Coerce from a fitted Stan model
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
post
#> MVBeliefUpdatr Lightweight Stanfit Posterior
#> Model Family     : NIW 
#> Group Context    : no_exposure, left_shifted, right_shifted 
#> Cues             : VOT 
#> Categories       : /d/, /t/ 
# }
```
