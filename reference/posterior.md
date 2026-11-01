# Compute posterior category probabilities for one or more observations

Compute posterior category probabilities for one or more observations

## Usage

``` r
posterior(x, new_data, categories, ...)

## S7 method for classes <MVBU_Stanfit>/<MVBU_StanfitPosterior>, <any>, <any>
posterior(x, new_data, categories, ...)
```

## Arguments

- x:

  A cognitive model object.

- new_data:

  A numeric matrix of observations or a list of matrices.

- categories:

  An optional subset of categories to restrict the posterior to.

- ...:

  Additional arguments passed to methods.

## Value

A matrix of posterior category probabilities with one row per
observation.

## Examples

``` r
mod <- example_mvg_ideal_observer(n_cues = 2)
posterior(mod, matrix(c(0, 0), nrow = 1))
#>            /b/       /d/       /g/        /k/        /p/         /t/
#> [1,] 0.4315091 0.3092369 0.2205673 0.01204114 0.02026484 0.006380768
```
