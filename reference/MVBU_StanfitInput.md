# Specify control parameters for new_ideal_adaptor_stanfit_input()

This function is used to specify control parameters for the
[`new_ideal_adaptor_stanfit_input()`](https://hlplab.github.io/MVBeliefUpdatr/reference/new_ideal_adaptor_stanfit_input.md)
function, and to provide reasonable defaults for any of the unspecified
parameters.

An S7 Class for Ideal-Adaptor Fit-Input Objects

## Usage

``` r
control_staninput(
  split_loglik_per_observation = 0,
  transform_type = c("identity", "center", "standardize", "PCA whiten", "ZCA whiten")[3]
)

IdealAdaptorStanfitInput(
  data = list(),
  transform_information = NULL,
  metadata = list()
)
```

## Arguments

- split_loglik_per_observation:

  Split log likelihood per observation.

- transform_type:

  Affine transformation type applied to the data.

- data:

  A named `list` with `$original` and `$transformed` data frames, or a
  single `data.frame` containing input data. Default:
  [`list()`](https://rdrr.io/r/base/list.html).

- transform_information:

  An
  [`MVBU_TransformInformation`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-TransformInformation-class.md)
  object. Default: `NULL`.

- metadata:

  A named `list` containing `label_information`. Default:
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

A list of control parameters.

## Examples

``` r
# Inspect an example IdealAdaptorStanfitInput object
sfi <- example_ideal_adaptor_stanfit_input()
sfi
#> <IdealAdaptorStanfitInput>
#>   Categories (2): /d/, /t/
#>   Cues (1): VOT
#>   Groups (3): no_exposure, left_shifted, right_shifted
#>   Data: 1480 rows

# Control parameters for stanfit input preparation
ctrl <- control_staninput(transform_type = "standardize")
ctrl
#> $split_loglik_per_observation
#> [1] 0
#> 
#> $transform_type
#> [1] "standardize"
#> 
```
