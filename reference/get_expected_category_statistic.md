# Extract expected category statistics from an S7 model, representation, or fit object

Computes the expected value of category parameters such as mean vectors
(\\\mu\\) or covariance/scatter matrices (\\\Sigma\\ / \\S\\) under the
model's distribution.

## Usage

``` r
get_expected_category_statistic(
  x,
  categories = NULL,
  groups = NULL,
  statistic = c("mu", "Sigma"),
  ...
)

get_marginal_category_statistic(
  x,
  categories = NULL,
  groups = NULL,
  statistic = c("mu", "Sigma"),
  ...
)

## S7 method for class <MVBU_StanfitPosterior>
get_expected_category_statistic(
  x,
  categories = NULL,
  groups = NULL,
  statistic = c("mu", "Sigma"),
  ...
)

## S7 method for class <MVBU_StanfitPosterior>
get_marginal_category_statistic(
  x,
  categories = NULL,
  groups = NULL,
  statistic = c("mu", "Sigma"),
  ...
)

## S7 method for class <MVBU_Stanfit>
get_expected_category_statistic(
  x,
  categories = NULL,
  groups = NULL,
  statistic = c("mu", "Sigma"),
  ...
)

## S7 method for class <MVBU_Stanfit>
get_marginal_category_statistic(
  x,
  categories = NULL,
  groups = NULL,
  statistic = c("mu", "Sigma"),
  ...
)
```

## Arguments

- x:

  An S7 object (cognitive model, category representation, template, or
  stanfit).

- categories:

  Optional character vector of category names. Default: `NULL` (all
  categories).

- groups:

  Optional character vector of group names. Default: `NULL` (all
  groups).

- statistic:

  Character scalar indicating the statistic to extract: `"mu"` (expected
  mean), `"Sigma"` or `"sigma"` (expected covariance), `"m"`, `"S"`,
  `"kappa"`, or `"nu"`.

- ...:

  Additional arguments passed to methods.

## Value

A numeric vector, matrix, or list of expected values.

## See also

[`get_expected_mu`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_parameters.md),
[`get_expected_sigma`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_parameters.md),
`get_marginal_category_statistic`,
[`get_parameters`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_parameters.md)

## Examples

``` r
mod <- example_mvg_ideal_observer(n_cues = 2)
get_expected_category_statistic(mod, statistic = "mu")
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
get_expected_category_statistic(mod, statistic = "Sigma")
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
get_marginal_category_statistic(mod, statistic = "mu")
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
```
