# Plot Stan MCMC diagnostics

Plot Stan MCMC diagnostics

## Usage

``` r
plot_diagnostics(x, ...)

## S7 method for class <MVBU_Stanfit>
plot_diagnostics(x, ...)
```

## Arguments

- x:

  An MVBU object
  ([`MVBU_CategoryRepresentation`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md),
  [`MVBU_CategoryRepresentationTemplate`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md),
  [`MVBU_CognitiveModel`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md),
  or
  [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)).

- ...:

  Additional arguments passed to specific plotting methods:

  - `uncertainty_treatment`: Treatment of posterior parameter
    uncertainty in fitted Stan models (`"marginalize"` or `"sample"`).

  - `ndraws`: Integer specifying the number of posterior MCMC draws to
    use from
    [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
    (defaults to `20L`).

  - `show_exposure_data`, `show_test_data`: Logical; if `TRUE`, overlays
    observed exposure/test data points.

  - `parallel`, `n_cores`: Logical and integer controlling parallel grid
    evaluation.

  - `noise_treatment`: Optional character string specifying perceptual
    noise treatment (`"no_noise"`, `"sample"`, or `"marginalize"`).

  - `interactive`: Logical; if `TRUE`, produces interactive 3D WebGL
    scenes via Plotly.

  - `slices`, `slice_cue`: Configuration for 3D sliced category plots.

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
plot_diagnostics(fit)
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
#> Warning: Removed 3 rows containing non-finite outside the scale range (`stat_bin()`).

# }
```
