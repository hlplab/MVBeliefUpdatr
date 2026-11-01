# Fit ideal adaptor

Infers a prior and posterior distribution of ideal adaptors from the
input data using Stan. Currently, three types of ideal adaptor models
are available, each using a conjugate prior over Gaussian or
multivariate Gaussian categories:

## Usage

``` r
fit_ideal_adaptor(
  stanfit_input,
  stanmodel = "NIW_ideal_adaptor",
  fixed_parameters = NULL,
  priors = NULL,
  file = NULL,
  file_refit = "never",
  file_compress = TRUE,
  stanvars = NULL,
  backend = "rstan",
  save_pars = NULL,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  init = "random",
  control = NULL,
  silent = 1,
  threads = NULL,
  stan_model_args = list(),
  rename = TRUE,
  ...
)
```

## Arguments

- stanfit_input:

  An object of class
  [`IdealAdaptorStanfitInput`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_StanfitInput.md)
  containing the prepared data and transform metadata.

- stanmodel:

  Character string naming the Stan model to use. Defaults to
  `"NIW_ideal_adaptor"`.

- fixed_parameters:

  An optional
  [`MVBU_FixedParameter`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_FixedParameter.md)
  object specifying fixed cognitive parameters created via
  [`fix_parameter`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_FixedParameter.md).
  Supports fixing global parameters (e.g., `"lapse_rate"`,
  `"category_prior"`, `"weights"`), prior hyperparameters (`"m_0"`,
  `"S_0"`, `"kappa_0"`, `"nu_0"`, `"expected_Sigma_0"`), and truly known
  category moments (`"mu_0"`, `"Sigma_0"`). When fixing `"mu_0"` and/or
  `"Sigma_0"`, the model automatically selects the appropriate conjugate
  update and predictive likelihood (Normal vs. Student-\\t\\); see
  [`fix_parameter`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_FixedParameter.md)
  for details. Default: `NULL`.

- priors:

  An optional
  [`MVBU_Prior`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Prior.md)
  object specifying prior distributions created via
  [`parameter_prior`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Prior.md).
  Default: `NULL`.

- file:

  Either NULL or a character string. In the latter case, the fitted
  model object is saved via saveRDS in a file named after the string
  supplied in file. The .rds extension is added automatically. If the
  file already exists, the model from that file will be loaded and
  returned instead of refitting the model. As existing files won't be
  overwritten, you have to manually remove the file in order to refit
  and save the model under an existing file name. The file name is
  stored in the `ideal_adaptor_stanfit` object for later use. (default:
  `NULL`)

- file_refit:

  Modifies when the fit stored via the `file` argument is re-used. For
  `"never"` (default) the fit is always loaded if it exists and fitting
  is skipped. For `"always"` the model is always refitted. If set to
  `"on_change"`, model will be refit if data passed to Stan differ from
  what is stored in the file. Refit will not be triggered for changes in
  additional parameters of the fit (e.g., initial values, number of
  iterations, control arguments, ...).

- file_compress:

  Logical or a character string, specifying one of the compression
  algorithms supported by
  [`saveRDS`](https://rdrr.io/r/base/readRDS.html). If the `file`
  argument is provided, this compression will be used when saving the
  fitted model object.

- stanvars:

  Optional `stanvars` object. Default: `NULL`.

- backend:

  Character string naming the package to use as the backend for fitting
  the Stan model. Options are `"rstan"` (the default) or `"cmdstanr"`.
  Details on the rstan and cmdstanr packages are available at
  <https://mc-stan.org/rstan/> and <https://mc-stan.org/cmdstanr/>,
  respectively.

- save_pars:

  Storage for parameter names to save. Default: `NULL`.

- chains:

  Number of Markov chains. Default: 4.

- iter:

  Total number of iterations per chain. Default: 2000.

- warmup:

  Number of warmup (burnin) iterations per chain. Default: 1000.

- init:

  Initial values specification. Default: `"random"`.

- control:

  Named list of control parameters passed to Stan sampler. Default:
  `NULL`.

- silent:

  Verbosity level between `0` and `2`. If `1` (the default), most of the
  informational messages of compiler and sampler are suppressed. If `2`,
  even more messages are suppressed. The actual sampling progress is
  still printed. Set `refresh = 0` to turn this off as well. If using
  `backend = "rstan"` you can also set `open_progress = FALSE` to
  prevent opening additional progress bars.

- threads:

  Optional number of threads per chain (or list containing `threads`)
  for multithreading within Stan. Default: `NULL`.

- stan_model_args:

  A `list` of further arguments passed to
  [rstan::stan_model](https://mc-stan.org/rstan/reference/stan_model.html)
  for `backend = "rstan"` or to
  [`cmdstanr::cmdstan_model`](https://mc-stan.org/cmdstanr/reference/cmdstan_model.html)
  for `backend = "cmdstanr"`, which allows to change how models are
  compiled.

- rename:

  For internal use only.

- ...:

  Additional parameters are passed to
  [rstan::sampling](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html),
  including the following:

## Value

An object of class `ideal_adaptor_stanfit` with the fitted stan model.

## Details

- `NIW_ideal_adaptor`:

  A Normal-Inverse-Wishart (NIW) prior over the ideal adaptor. This is
  the default model with multivariate Gaussian categories. Accepts
  univariate and multivariate input, though the NIX model should be
  faster for univariate input.

- `NIX_ideal_adaptor`:

  A Normal-Inverse-Chisquare (NIX) prior over the ideal adaptor with
  univariate Gaussian categories. Accepts only univariate input.

- `MNIX_ideal_adaptor`:

  Separate NIXs for each of multiple cues that are integrated over
  during categorization assuming ideal cue weights based on the relative
  informativity of each cue.

In interpreting the inferred parameters, it should be kept in mind that
the *inferred* scatter matrix (e.g., S_0 for the NIW model) includes
variability from internal perceptual and/or external environmental
noise, *in addition* to the motor noise that is reflected in production
data. This also implies that, **if `Sigma_0` is provided by the user it
should be arguably convolved with an estimate of perceptual noise**.
This is particularly important if the data you're fitting contains test
phases without exposure (e.g., pre-exposure tests).

Use `methods(class = "ideal_adaptor_stanfit")` for an overview on
available methods.

## See also

[`IdealAdaptorStanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
for information about ideal adaptor stanfit objects,
[`get_draws`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_draws.md)
to draw samples from the stanfit.

## Examples

``` r
# \donttest{
# Prepare minimal input data
exp_df <- data.frame(
  group = "g1",
  category = c("b", "p"),
  vot = c(0, 50)
)
test_df <- data.frame(
  group = "g1",
  response_category = c("b", "p"),
  vot = c(10, 40)
)
sfi <- new_ideal_adaptor_stanfit_input(
  exposure = exp_df,
  test = test_df,
  cues = "vot"
)

# Fit NIW ideal adaptor model
fit <- fit_ideal_adaptor(
  stanfit_input = sfi,
  chains = 1,
  iter = 50,
  warmup = 25
)
#> 
#> SAMPLING FOR MODEL 'NIW_ideal_adaptor' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 6.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.63 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: There aren't enough warmup iterations to fit the
#> Chain 1:          three stages of adaptation as currently configured.
#> Chain 1:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 1:          the given number of warmup iterations:
#> Chain 1:            init_buffer = 3
#> Chain 1:            adapt_window = 20
#> Chain 1:            term_buffer = 2
#> Chain 1: 
#> Chain 1: Iteration:  1 / 50 [  2%]  (Warmup)
#> Chain 1: Iteration:  5 / 50 [ 10%]  (Warmup)
#> Chain 1: Iteration: 10 / 50 [ 20%]  (Warmup)
#> Chain 1: Iteration: 15 / 50 [ 30%]  (Warmup)
#> Chain 1: Iteration: 20 / 50 [ 40%]  (Warmup)
#> Chain 1: Iteration: 25 / 50 [ 50%]  (Warmup)
#> Chain 1: Iteration: 26 / 50 [ 52%]  (Sampling)
#> Chain 1: Iteration: 30 / 50 [ 60%]  (Sampling)
#> Chain 1: Iteration: 35 / 50 [ 70%]  (Sampling)
#> Chain 1: Iteration: 40 / 50 [ 80%]  (Sampling)
#> Chain 1: Iteration: 45 / 50 [ 90%]  (Sampling)
#> Chain 1: Iteration: 50 / 50 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.027 seconds (Warm-up)
#> Chain 1:                0.04 seconds (Sampling)
#> Chain 1:                0.067 seconds (Total)
#> Chain 1: 
#> Warning: There were 2 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
fit
#> <NIW_IdealAdaptorStanfit>
#>   Model family: NIW
#>   Categories (2): b, p
#>   Cues (1): vot
#>   Groups (1): g1
#>   Draws: 25
# }
```
