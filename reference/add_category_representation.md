# Add a category representation to a template, model, or model list

Adds a category representation to an existing template, cognitive model,
or each model in a model list. When adding to a cognitive model,
category priors and lapse biases for the new representation default to 0
(rescaling existing categories proportionally), or can be specified for
the new category or across all categories.

## Usage

``` r
add_category_representation(
  x,
  representation,
  name = NULL,
  category_prior = NULL,
  lapse_bias = NULL,
  ...
)

## S7 method for class <MVBU_CategoryRepresentationTemplate>
add_category_representation(
  x,
  representation,
  name = NULL,
  category_prior = NULL,
  lapse_bias = NULL,
  ...
)

## S7 method for class <MVBU_CognitiveModel>
add_category_representation(
  x,
  representation,
  name = NULL,
  category_prior = NULL,
  lapse_bias = NULL,
  ...
)
```

## Arguments

- x:

  An
  [MVBU_CategoryRepresentationTemplate](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md),
  [MVBU_CognitiveModel](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md),
  or
  [MVBU_ModelList](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_ModelList.md)
  object.

- representation:

  An
  [MVBU_CategoryRepresentation](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md)
  object to add.

- name:

  Optional character scalar specifying the name of the category.
  Default: `NULL`.

- category_prior:

  Optional numeric vector of category prior probabilities. Can be a
  single value for the new category (defaulting to 0, which rescales
  existing categories by `previous_value * (1 - category_prior)`), or a
  full vector of prior probabilities across all categories (including
  the new one) summing to 1. Default: `NULL`.

- lapse_bias:

  Optional numeric vector of lapse bias probabilities. Can be a single
  value for the new category (defaulting to 0, which rescales existing
  categories by `previous_value * (1 - lapse_bias)`), or a full vector
  of lapse biases across all categories (including the new one) summing
  to 1. Default: `NULL`.

- ...:

  Additional arguments passed to methods.

## Value

An updated object of the same class with the category representation
added.

## See also

[`get_category_labels()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_labels.md)

## Examples

``` r
mod <- example_mvg_ideal_observer(n_cues = 2)
new_rep <- new_mvg_category_representation(
  category_labels = "/s/",
  cue_labels = c("VOT", "f0_semitones"),
  mu = c(30, 0),
  Sigma = diag(c(25, 10))
)
mod_new <- add_category_representation(mod, new_rep, name = "/s/")
get_category_labels(mod_new)
#> [1] "/b/" "/d/" "/g/" "/k/" "/p/" "/t/" "/s/"
```
