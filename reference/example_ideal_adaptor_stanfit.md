# Example IdealAdaptorStanfit object

Fits or loads a pre-fitted example
[`IdealAdaptorStanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
model. If `file` is specified, the model is loaded from or saved to disk
using
[`fit_ideal_adaptor`](https://hlplab.github.io/MVBeliefUpdatr/reference/fit_ideal_adaptor.md)
according to `file_refit`.

## Usage

``` r
example_ideal_adaptor_stanfit(
  model_family = c("NIW", "NIX", "MNIX"),
  n_cues = NULL,
  categories = c("/d/", "/t/"),
  seed = 42L,
  staninput_control = MVBeliefUpdatr::control_staninput(),
  file = NULL,
  file_refit = c("on_change", "never", "always"),
  stanmodel = NULL,
  chains = 1L,
  iter = 50L,
  warmup = 25L,
  refresh = 0L,
  ...
)
```

## Arguments

- model_family:

  Character string indicating the model family: `"NIW"` (default),
  `"NIX"`, or `"MNIX"`.

- n_cues:

  Integer number of acoustic cues: 1 (`"VOT"`), 2 (`"VOT"` and
  `"f0_semitones"`), or 3 (`"VOT"`, `"f0_semitones"`, and
  `"vowel_duration"`). Default: 1 (or 2 for `"MNIX"`).

- categories:

  Character vector of length 2 naming the categories. Default:
  `c("/d/", "/t/")`.

- seed:

  Random seed for reproducible data generation and sampling. Default:
  42L.

- staninput_control:

  Named list of Stan input control parameters created by
  [`control_staninput`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_StanfitInput.md).
  Default:
  [`control_staninput()`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_StanfitInput.md).

- file:

  Either `NULL` or a character string naming the file (without or with
  `.rds` extension) to load or save the model.

- file_refit:

  Character string controlling cache behavior when `file` is specified:
  `"on_change"` (default), `"never"`, or `"always"`. See
  [`fit_ideal_adaptor`](https://hlplab.github.io/MVBeliefUpdatr/reference/fit_ideal_adaptor.md).

- stanmodel:

  Optional character string naming the Stan model. Defaults to
  `"MNIX_ideal_adaptor_ideal_integrator"` when `model_family == "MNIX"`
  and `"<model_family>_ideal_adaptor"` otherwise.

- chains:

  Number of Markov chains. Default: 1L.

- iter:

  Total number of iterations per chain. Default: 50L.

- warmup:

  Number of warmup iterations per chain. Default: 25L.

- refresh:

  Frequency of progress reports during sampling (0 disables). Default:
  0L.

- ...:

  Additional arguments passed to
  [`fit_ideal_adaptor`](https://hlplab.github.io/MVBeliefUpdatr/reference/fit_ideal_adaptor.md)
  (e.g., sampler arguments such as
  `control = list(adapt_delta = 0.99)`).

## Value

An object of class
[`IdealAdaptorStanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md).

## See also

[`example_ideal_adaptor_stanfit_input`](https://hlplab.github.io/MVBeliefUpdatr/reference/example_ideal_adaptor_stanfit_input.md),
[`fit_ideal_adaptor`](https://hlplab.github.io/MVBeliefUpdatr/reference/fit_ideal_adaptor.md)

## Examples

``` r
# \donttest{
# Fit a quick 1-cue NIX ideal adaptor (wrapped in \donttest for CRAN)
fit <- example_ideal_adaptor_stanfit(
  model_family = "NIX",
  n_cues = 1
)
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
fit
#> <NIX_IdealAdaptorStanfit>
#>   Model family: NIX
#>   Categories (2): /d/, /t/
#>   Cues (1): VOT
#>   Groups (3): no_exposure, left_shifted, right_shifted
#>   Draws: 25
# }
```
