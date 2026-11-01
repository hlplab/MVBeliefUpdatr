# Compute category likelihoods for one or more observations

For NIX, MNIX, and NIW families the likelihood is the posterior
predictive of the category, so this generic supersedes the
family-specific posterior predictive helpers.

## Usage

``` r
likelihood(x, new_data, categories, ...)

## S7 method for classes <MVBU_Stanfit>/<MVBU_StanfitPosterior>, <any>, <any>
likelihood(x, new_data, categories, ...)
```

## Arguments

- x:

  A category representation, representation template, or cognitive
  model.

- new_data:

  A numeric matrix of observations, a vector, or a list of observations.

- categories:

  An optional subset of categories to restrict the result to.

- ...:

  Additional arguments passed to methods.

## Value

A matrix of category likelihoods with one row per observation and one
column per category.

## Examples

``` r
mod <- example_mvg_ideal_observer(n_cues = 2)
likelihood(mod, matrix(c(0, 0), nrow = 1))
#>              /b/         /d/         /g/         /k/          /p/          /t/
#> [1,] 0.001724177 0.001235615 0.000881319 4.81127e-05 8.097206e-05 2.549558e-05
```
