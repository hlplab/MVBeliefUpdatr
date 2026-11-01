# Deprecated: get_priors_from_model

**\[deprecated\]** `get_priors_from_model()` was deprecated in
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
[`get_category_prior()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_category_prior.md)
instead.

## Usage

``` r
get_priors_from_model(model, categories = NULL)
```

## Arguments

- model:

  A model object.

- categories:

  Optional vector of category labels to subset.

## Value

Numeric vector of prior probabilities.

## See also

[`get_category_prior()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_category_prior.md)
