# Add Pre-extracted Posterior Latents and Evaluation Closures

Extracts and structures draw-level latent parameters (\\m\\, \\S\\,
\\\kappa\\, \\\nu\\, \\\Sigma\_{\text{exp}}\\,
\\\Sigma\_{\text{marg}}\\, and lapse rates) and compiles unified
likelihood and posterior evaluation closures. Pre-computed outputs are
stored strictly in `model@cache$draws` and `model@cache$summary`:

- `model@cache$draws`: Contains draw-level representations:

  - `parameters`:

    A nested tibble of raw posterior draws across iterations/draws,
    including \\m\\, \\S\\, \\\kappa\\, \\\nu\\,
    \\\Sigma\_{\text{exp}}\\, \\\Sigma\_{\text{marg}}\\, and optional
    `lapse_rate`.

  - `likelihood_function`:

    Closure evaluating draw-wise category likelihoods with model-family
    predictive densities (NIW: multivariate Student-t; NIX: univariate
    Student-t; MNIX: independent Student-t).

  - `posterior_function`:

    Closure evaluating draw-wise posterior category probabilities.

- `model@cache$summary`: Contains point estimates and across-draw
  summaries:

  - `parameters`:

    Point estimates obtained by taking the mean across all posterior
    draws for each latent parameter (\\m\\, \\S\\, \\\kappa\\, \\\nu\\,
    \\\Sigma\_{\text{exp}}\\, \\\Sigma\_{\text{marg}}\\, and
    `lapse_rate`) per group and category.

  - `expected_moments`:

    Distribution summary table (mean, standard deviation, and quantile
    credible intervals specified by `probs`) across posterior draws for
    the category mean \\\mu\\ and expected category covariance
    \\\Sigma\_{\text{exp}} = S / (\nu - D - 1)\\.

  - `marginal_moments`:

    Distribution summary table across posterior draws for the category
    mean \\\mu\\ and marginal predictive covariance
    \\\Sigma\_{\text{marg}} = \frac{\kappa + 1} {\kappa}
    \Sigma\_{\text{exp}}\\.

  - `likelihood_function`:

    Closure evaluating category likelihoods from aggregated point
    estimates.

  - `posterior_function`:

    Closure evaluating category posterior probabilities from aggregated
    point estimates.

## Usage

``` r
add_posterior_latents(model, probs = c(0.025, 0.5, 0.975), ...)
```

## Arguments

- model:

  An
  [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
  object.

- probs:

  Numeric vector of probabilities for credible interval quantiles.
  (default: `c(0.025, 0.5, 0.975)`)

- ...:

  Additional arguments passed to
  [`get_draws`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_draws.md).

## Value

The updated
[`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
object with pre-extracted latents stored in `model@cache$draws` and
`model@cache$summary`.

## See also

[`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md),
[`get_draws`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_draws.md),
[`add_criterion`](https://hlplab.github.io/MVBeliefUpdatr/reference/add_criterion.md),
[`vignette("fitting-and-working-with-stanfit-models")`](https://hlplab.github.io/MVBeliefUpdatr/articles/fitting-and-working-with-stanfit-models.md)

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
fit <- add_posterior_latents(fit)
fit
#> <NIW_IdealAdaptorStanfit>
#>   Model family: NIW
#>   Categories (2): /d/, /t/
#>   Cues (1): VOT
#>   Groups (3): no_exposure, left_shifted, right_shifted
#>   Draws: 25
# }
```
