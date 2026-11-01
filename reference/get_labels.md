# Get cue, category, response category, or group labels from an object

Query or modify the label mappings (cues, categories, response
categories, or participant groups) stored in representations, templates,
cognitive models, or Stan fit/input objects.

## Usage

``` r
get_labels(x, ...)

get_cue_labels(x, indices = NULL, ...)

get_category_labels(x, indices = NULL, ...)

get_response_category_labels(x, indices = NULL, ...)

get_group_labels(x, indices = NULL, include_prior = FALSE, ...)

set_labels(
  x,
  cue = NULL,
  category = NULL,
  response_category = NULL,
  group = NULL,
  ...
)

## S7 method for class <MVBU_ModelList>
get_cue_labels(x, indices = NULL, ...)

## S7 method for class <MVBU_ModelList>
get_category_labels(x, indices = NULL, ...)

## S7 method for class <list>
set_labels(
  x,
  cue = NULL,
  category = NULL,
  response_category = NULL,
  group = NULL,
  ...
)

## S7 method for class <MVBU_Object>
set_labels(
  x,
  cue = NULL,
  category = NULL,
  response_category = NULL,
  group = NULL,
  ...
)

## S7 method for class <MVBU_Object>
get_labels(x, ...)

## S7 method for class <MVBU_Object>
get_cue_labels(x, indices = NULL, ...)

## S7 method for class <MVBeliefUpdatr::MNIX_EarlyIdealAdaptor>
get_cue_labels(x, indices = NULL, ...)

## S7 method for class <MVBU_Object>
get_category_labels(x, indices = NULL, ...)

## S7 method for class <MVBU_Object>
get_response_category_labels(x, indices = NULL, ...)

## S7 method for class <MVBU_Object>
get_group_labels(x, indices = NULL, include_prior = FALSE, ...)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>
get_cue_labels(x, indices = NULL, ...)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>
get_category_labels(x, indices = NULL, ...)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>
get_response_category_labels(x, indices = NULL, ...)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>
get_group_labels(x, indices = NULL, include_prior = FALSE, ...)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>
get_labels(x, ...)
```

## Arguments

- x:

  A representation, representation template, cognitive model, stanfit,
  or fit-input object.

- ...:

  Additional arguments passed to methods.

- indices:

  An optional integer vector of indices to subset the returned labels.

- include_prior:

  Logical; whether to include the prior group label. Default: `FALSE`.

- cue:

  Character vector of cue labels. Defaults to `character(0)`.

- category:

  Character vector of category labels. Defaults to `character(0)`.

- response_category:

  Character vector of response category labels. Defaults to `category`.

- group:

  Character vector of group labels. Defaults to `character(0)`.

## Value

A character vector of labels, a named list of labels (`get_labels`), or
the modified object with updated labels (`set_labels`).

## See also

[`get_parameters`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_parameters.md)

## Examples

``` r
mod <- example_mvg_ideal_observer(n_cues = 2)
get_labels(mod)
#> $cue
#> [1] "VOT"          "f0_semitones"
#> 
#> $category
#> [1] "/b/" "/d/" "/g/" "/k/" "/p/" "/t/"
#> 
#> $response_category
#> [1] "/b/" "/d/" "/g/" "/k/" "/p/" "/t/"
#> 
#> $group
#> character(0)
#> 
get_cue_labels(mod)
#> [1] "VOT"          "f0_semitones"
get_category_labels(mod)
#> [1] "/b/" "/d/" "/g/" "/k/" "/p/" "/t/"
mod <- set_labels(mod, cue = c("c1", "c2"))
```
