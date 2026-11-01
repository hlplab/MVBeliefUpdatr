# Summarize an IdealAdaptorStanfitInput object

Provides a detailed overview of the exposure sufficient statistics
across group and category combinations as well as key properties of the
test dataset.

Specifies reasonable defaults for the parameters to be summarized for
the stanfit object.

## Usage

``` r
# S3 method for class 'IdealAdaptorStanfitInput'
summary(object, categories = NULL, groups = NULL, ...)

# S3 method for class 'MVBU_Stanfit'
summary(
  object,
  pars = NULL,
  sufficient_pars_only = TRUE,
  indices_as_names = TRUE,
  include_transformed_pars = FALSE,
  probs = c(0.025, 0.5, 0.975),
  ...
)
```

## Arguments

- object:

  An
  [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
  object.

- categories:

  Optional character vector of category labels to summarize. Defaults to
  all categories in `object`.

- groups:

  Optional character vector of group labels to summarize. Defaults to
  all groups in `object`.

- ...:

  Additional arguments passed to `rstan::summary`.

- pars:

  A character vector of parameter names to be summarized. If `NULL`, all
  model parameters are summarized. (default: `NULL`)

- sufficient_pars_only:

  Should only the sufficient parameters be summarized? (default: `TRUE`)

- indices_as_names:

  Should the indices of the parameters be translated into level names in
  the summary? This output is substantially more readable, with names
  shown in separate columns (distribution, group, category, cue).
  (default: `TRUE`)

- include_transformed_pars:

  Should transformed parameters be included in the summary? (default:
  `FALSE`)

- probs:

  Numeric vector of probabilities for credible interval quantiles.
  (default: `c(0.025, 0.5, 0.975)`)
