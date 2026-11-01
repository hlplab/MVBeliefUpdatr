# Extract parameters from models, templates, or representations

Extract parameters, parameter names, or priors from objects

## Usage

``` r
get_parameters(x, space = c("original", "transformed"), ...)

get_stan_priors(x, space = c("original", "transformed"), ...)

get_parameter_names(x, original_pars = FALSE, ...)

get_expected_mu(x, ...)

get_expected_sigma(x, ...)

get_marginal_mu(x, ...)

get_marginal_sigma(x, ...)

## S7 method for class <MVBU_Prior>
get_stan_priors(x, space = c("original", "transformed"), ...)

## S7 method for class <MVBU_Object>
get_parameters(x, space = c("original", "transformed"), ...)

## S7 method for class <MVBU_Object>
get_stan_priors(x, space = c("original", "transformed"), ...)

## S7 method for class <MVBU_Stanfit>
get_parameter_names(x, original_pars = FALSE, ...)

## S7 method for class <any>
get_parameter_names(x, original_pars = FALSE, ...)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>
get_expected_mu(x, ...)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>
get_expected_sigma(x, ...)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>
get_marginal_mu(x, ...)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>
get_marginal_sigma(x, ...)
```

## Arguments

- x:

  A category representation, template, cognitive model, stanfit, or
  prior object.

- space:

  Character string indicating which parameter space to return:
  `"original"` (default) or `"transformed"`.

- ...:

  Additional arguments passed to methods.

- original_pars:

  Logical scalar; if `TRUE`, return original (unmapped / top-level) Stan
  model parameter names from the underlying Stan program (as stored in
  `stanfit@model_pars`) rather than the expanded, indexed, and
  post-processed/renamed parameter draw names (as returned by
  `names(stanfit)`). Default: `FALSE`.

## Value

A list, vector, matrix, data frame, or
[`MVBU_Prior`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Prior.md)
object.

## See also

[`parameter_prior`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Prior.md),
[`get_data`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_data.md)

## Examples

``` r
mod <- example_mvg_ideal_observer(n_cues = 2)
get_parameters(mod)
#> $`/b/`
#> $`/b/`$mu
#> [1]  8.477551 -1.368047
#> 
#> $`/b/`$Sigma
#>                      VOT f0_semitones
#> VOT          23.20282907  -0.08870048
#> f0_semitones -0.08870048  14.68964368
#> 
#> 
#> $`/d/`
#> $`/d/`$mu
#> [1] 13.917921 -0.580694
#> 
#> $`/d/`$Sigma
#>                    VOT f0_semitones
#> VOT          75.914604    -1.445236
#> f0_semitones -1.445236    16.963832
#> 
#> 
#> $`/g/`
#> $`/g/`$mu
#> [1] 16.9919098 -0.7449837
#> 
#> $`/g/`$Sigma
#>                      VOT f0_semitones
#> VOT          107.3474463   -0.1125058
#> f0_semitones  -0.1125058   20.0938994
#> 
#> 
#> $`/k/`
#> $`/k/`$mu
#> [1] 54.544491  0.576098
#> 
#> $`/k/`$Sigma
#>                    VOT f0_semitones
#> VOT          409.69309    -11.39129
#> f0_semitones -11.39129     13.89484
#> 
#> 
#> $`/p/`
#> $`/p/`$mu
#> [1] 50.6026289 -0.4758416
#> 
#> $`/p/`$Sigma
#>                     VOT f0_semitones
#> VOT          408.755971    -6.736239
#> f0_semitones  -6.736239    17.966702
#> 
#> 
#> $`/t/`
#> $`/t/`$mu
#> [1] 61.2998157 -0.2026473
#> 
#> $`/t/`$Sigma
#>                     VOT f0_semitones
#> VOT          436.650466    -9.659716
#> f0_semitones  -9.659716    15.159686
#> 
#> 
get_parameter_names(mod)
#> [1] "mu"             "Sigma"          "category_prior" "lapse_rate"    
#> [5] "lapse_bias"     "Sigma_noise"   
get_expected_mu(mod)
#> $`/b/`
#> [1]  8.477551 -1.368047
#> 
#> $`/d/`
#> [1] 13.917921 -0.580694
#> 
#> $`/g/`
#> [1] 16.9919098 -0.7449837
#> 
#> $`/k/`
#> [1] 54.544491  0.576098
#> 
#> $`/p/`
#> [1] 50.6026289 -0.4758416
#> 
#> $`/t/`
#> [1] 61.2998157 -0.2026473
#> 
get_expected_sigma(mod)
#> $`/b/`
#>                      VOT f0_semitones
#> VOT          23.20282907  -0.08870048
#> f0_semitones -0.08870048  14.68964368
#> 
#> $`/d/`
#>                    VOT f0_semitones
#> VOT          75.914604    -1.445236
#> f0_semitones -1.445236    16.963832
#> 
#> $`/g/`
#>                      VOT f0_semitones
#> VOT          107.3474463   -0.1125058
#> f0_semitones  -0.1125058   20.0938994
#> 
#> $`/k/`
#>                    VOT f0_semitones
#> VOT          409.69309    -11.39129
#> f0_semitones -11.39129     13.89484
#> 
#> $`/p/`
#>                     VOT f0_semitones
#> VOT          408.755971    -6.736239
#> f0_semitones  -6.736239    17.966702
#> 
#> $`/t/`
#>                     VOT f0_semitones
#> VOT          436.650466    -9.659716
#> f0_semitones  -9.659716    15.159686
#> 
```
