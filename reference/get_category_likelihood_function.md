# Extract the category-likelihood function from a representation, template, or model

For NIX, MNIX, and NIW families the returned function evaluates the
posterior predictive, which is the category likelihood for those
families.

## Usage

``` r
get_category_likelihood_function(x, ...)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>
get_category_likelihood_function(x, ...)
```

## Arguments

- x:

  A category representation, representation template, or cognitive
  model.

- ...:

  Additional arguments passed to methods.

## Value

A function of `(new_data, log, noise_treatment, Sigma_noise)` returning
category likelihoods.

## Examples

``` r
mod <- example_mvg_ideal_observer(n_cues = 2)
lik_fn <- get_category_likelihood_function(mod)
lik_fn(matrix(c(0, 0), nrow = 1))
#>              /b/         /d/         /g/         /k/          /p/          /t/
#> [1,] 0.001724177 0.001235615 0.000881319 4.81127e-05 8.097206e-05 2.549558e-05
```
