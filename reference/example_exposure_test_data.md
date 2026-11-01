# Example exposure and test experimental data

Generates simulated exposure-test experiment datasets for phonetic
categorization and belief updating models. Shifts between exposure
conditions hold the mean distance between categories constant, while
variances remain unchanged. Shift magnitudes are approximately 50\\
distance between prior category means. A `"no_exposure"` baseline
condition tests the unshifted prior. Test locations form a regular grid
spanning the category continuum and the adaptation shift (forming a
parallelogram in 2D space), with simulated categorization responses.

## Usage

``` r
example_exposure_test_data(
  model_family = c("NIW", "NIX", "MNIX"),
  n_cues = NULL,
  categories = c("/d/", "/t/"),
  n_exposure_per_category = 20L,
  n_test_steps_category = 5L,
  n_test_steps_adaptation = 3L,
  n_test_tokens_per_location = 40L,
  seed = 42L
)
```

## Arguments

- model_family:

  Character string indicating the model family to simulate: `"NIW"`
  (default), `"NIX"`, or `"MNIX"`.

- n_cues:

  Integer number of acoustic cues (1 for NIX, 2-3 for MNIX, 1 to 3 for
  NIW). Default: 1 (or 2 for MNIX).

- categories:

  Character vector of length 2 naming the categories. Default:
  `c("/d/", "/t/")`.

- n_exposure_per_category:

  Integer number of exposure observations per category for each shifted
  condition. Default: 20L.

- n_test_steps_category:

  Integer number of test continuum steps between categories. Default:
  5L.

- n_test_steps_adaptation:

  Integer number of test continuum steps between prior and posterior
  means. Default: 3L.

- n_test_tokens_per_location:

  Integer number of test trials evaluated at each test location.
  Default: 40L.

- seed:

  Optional integer random seed for reproducibility. Default: 42L.

## Value

A list with two tibbles:

- exposure:

  Exposure trials containing `Phase`, `Condition`, `group`, `Subject`,
  cue columns, and `category`.

- test:

  Test trials containing `Phase`, `Condition`, `group`, `Subject`, cue
  columns, and `response_category`.

## See also

[`example_ideal_adaptor_stanfit_input`](https://hlplab.github.io/MVBeliefUpdatr/reference/example_ideal_adaptor_stanfit_input.md),
[`example_ideal_adaptor_staninput`](https://hlplab.github.io/MVBeliefUpdatr/reference/example_ideal_adaptor_stanfit_input.md),
[`example_model`](https://hlplab.github.io/MVBeliefUpdatr/reference/example-s7-objects.md)

## Examples

``` r
# Generate small exposure and test experiment datasets
exp_test <- example_exposure_test_data(
  model_family = "NIX",
  n_cues = 1,
  n_exposure_per_category = 5L,
  n_test_steps_category = 3L,
  n_test_steps_adaptation = 2L,
  n_test_tokens_per_location = 5L
)
head(exp_test$exposure)
#> # A tibble: 6 × 6
#>   Phase    Condition    group        Subject              VOT category
#>   <fct>    <fct>        <fct>        <chr>              <dbl> <fct>   
#> 1 exposure left_shifted left_shifted left_shifted_sub1   4.18 /d/     
#> 2 exposure left_shifted left_shifted left_shifted_sub1 -20.0  /d/     
#> 3 exposure left_shifted left_shifted left_shifted_sub1  -8.18 /d/     
#> 4 exposure left_shifted left_shifted left_shifted_sub1 -19.1  /d/     
#> 5 exposure left_shifted left_shifted left_shifted_sub1   4.69 /d/     
#> 6 exposure left_shifted left_shifted left_shifted_sub1  37.9  /t/     
head(exp_test$test)
#> # A tibble: 6 × 6
#>   Phase Condition   group       Subject            VOT response_category
#>   <fct> <fct>       <fct>       <chr>            <dbl> <fct>            
#> 1 test  no_exposure no_exposure no_exposure_sub1  13.8 /d/              
#> 2 test  no_exposure no_exposure no_exposure_sub1  13.8 /d/              
#> 3 test  no_exposure no_exposure no_exposure_sub1  13.8 /d/              
#> 4 test  no_exposure no_exposure no_exposure_sub1  13.8 /d/              
#> 5 test  no_exposure no_exposure no_exposure_sub1  13.8 /d/              
#> 6 test  no_exposure no_exposure no_exposure_sub1  37.1 /t/              
```
