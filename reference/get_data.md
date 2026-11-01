# Get original variable names from an S7 model or stanfit input object

Retrieves the original variable names that were mapped to standard roles
(`"group"`, `"group_unique"`, `"category"`, `"response_category"`, and
`"cues"`) when the model or input object was created.

Retrieves underlying experimental data (all observations across phases,
exposure observations only, or test observations only) stored in an S7
[`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md),
[`IdealAdaptorStanfitInput`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_StanfitInput.md),
or
[`MVBU_StanfitPosterior`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_StanfitPosterior.md)
object.

## Usage

``` r
get_original_variable_names(
  x,
  variable = c("group", "group_unique", "category", "response_category", "cues"),
  ...
)

get_data(
  x,
  space = c("original", "transformed"),
  groups = NULL,
  categories = NULL,
  response_categories = NULL,
  n_samples = NULL,
  original_names = FALSE,
  ...
)

get_exposure_data(
  x,
  space = c("original", "transformed"),
  groups = NULL,
  categories = NULL,
  n_samples = NULL,
  original_names = FALSE,
  ...
)

get_test_data(
  x,
  space = c("original", "transformed"),
  groups = NULL,
  response_categories = NULL,
  n_samples = NULL,
  original_names = FALSE,
  ...
)

## S7 method for class <IdealAdaptorStanfitInput>
get_data(
  x,
  space = c("original", "transformed"),
  groups = NULL,
  categories = NULL,
  response_categories = NULL,
  n_samples = NULL,
  original_names = FALSE,
  ...
)

## S7 method for class <IdealAdaptorStanfitInput>
get_exposure_data(
  x,
  space = c("original", "transformed"),
  groups = NULL,
  categories = NULL,
  n_samples = NULL,
  original_names = FALSE,
  ...
)

## S7 method for class <IdealAdaptorStanfitInput>
get_test_data(
  x,
  space = c("original", "transformed"),
  groups = NULL,
  response_categories = NULL,
  n_samples = NULL,
  original_names = FALSE,
  ...
)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>
get_original_variable_names(
  x,
  variable = c("group", "group_unique", "category", "response_category", "cues"),
  ...
)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>
get_data(
  x,
  space = c("original", "transformed"),
  groups = NULL,
  categories = NULL,
  response_categories = NULL,
  n_samples = NULL,
  original_names = FALSE,
  ...
)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>
get_exposure_data(
  x,
  space = c("original", "transformed"),
  groups = NULL,
  categories = NULL,
  n_samples = NULL,
  original_names = FALSE,
  ...
)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>
get_test_data(
  x,
  space = c("original", "transformed"),
  groups = NULL,
  response_categories = NULL,
  n_samples = NULL,
  original_names = FALSE,
  ...
)
```

## Arguments

- x:

  An S7 stanfit, stanfit-input, or stanfit-posterior object.

- variable:

  Character vector subsetting the variable types to retrieve. Must be
  one or more of `"group"`, `"group_unique"`, `"category"`,
  `"response_category"`, or `"cues"`. Defaults to all of them.

- ...:

  Additional arguments passed to methods.

- space:

  Character string indicating which cue space to return: `"original"`
  (default) or `"transformed"`.

- groups:

  Optional character vector of group labels to filter observations by.
  Defaults to all available groups.

- categories:

  Optional character vector of category labels to filter exposure
  observations by. Defaults to all available category levels.

- response_categories:

  Optional character vector of response category labels to filter test
  observations by. Defaults to all available response category levels.

- n_samples:

  Optional integer; subsample up to `n_samples` rows. Defaults to `NULL`
  (all observations).

- original_names:

  Logical; if `TRUE`, renames standardized columns back to original
  variable names. Defaults to `FALSE`.

## Value

A named list of character vectors (or a character vector if a single
variable is requested) containing the original variable names.

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html)
containing the requested data.

## See also

`get_data`, `get_exposure_data`, `get_test_data`

`get_original_variable_names`

## Examples

``` r
sfi <- example_ideal_adaptor_stanfit_input()
get_original_variable_names(sfi)
#> $group
#> [1] "group"
#> 
#> $group_unique
#> [1] "group"
#> 
#> $category
#> [1] "category"
#> 
#> $response_category
#> [1] "response_category"
#> 
#> $cues
#> [1] "VOT"
#> 
head(get_data(sfi))
#> # A tibble: 6 × 6
#>   group        group_unique category response_category     VOT Phase   
#>   <fct>        <fct>        <fct>    <chr>               <dbl> <fct>   
#> 1 left_shifted left_shifted /d/      NA                  1.06  exposure
#> 2 left_shifted left_shifted /d/      NA                 -0.317 exposure
#> 3 left_shifted left_shifted /d/      NA                -18.4   exposure
#> 4 left_shifted left_shifted /d/      NA                  6.77  exposure
#> 5 left_shifted left_shifted /d/      NA                -15.4   exposure
#> 6 left_shifted left_shifted /d/      NA                 -8.59  exposure
head(get_exposure_data(sfi))
#> # A tibble: 6 × 6
#>   group        group_unique category response_category     VOT Phase   
#>   <fct>        <fct>        <fct>    <chr>               <dbl> <fct>   
#> 1 left_shifted left_shifted /d/      NA                  1.06  exposure
#> 2 left_shifted left_shifted /d/      NA                 -0.317 exposure
#> 3 left_shifted left_shifted /d/      NA                -18.4   exposure
#> 4 left_shifted left_shifted /d/      NA                  6.77  exposure
#> 5 left_shifted left_shifted /d/      NA                -15.4   exposure
#> 6 left_shifted left_shifted /d/      NA                 -8.59  exposure
head(get_test_data(sfi))
#> # A tibble: 6 × 6
#>   group       group_unique category response_category   VOT Phase
#>   <fct>       <fct>        <fct>    <chr>             <dbl> <fct>
#> 1 no_exposure no_exposure  NA       /d/                13.8 test 
#> 2 no_exposure no_exposure  NA       /d/                13.8 test 
#> 3 no_exposure no_exposure  NA       /d/                13.8 test 
#> 4 no_exposure no_exposure  NA       /d/                13.8 test 
#> 5 no_exposure no_exposure  NA       /d/                13.8 test 
#> 6 no_exposure no_exposure  NA       /d/                13.8 test 
```
