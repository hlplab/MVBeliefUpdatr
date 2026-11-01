# S7 class for ordered collections of cognitive models

Represents an ordered collection or list of cognitive model objects
(which may represent sequential update histories, parameter variations,
or distinct model variants).

## Usage

``` r
MVBU_ModelList(
  models = list(),
  model_labels = character(0),
  metadata = list(parameters = list())
)

as_model_list(models, model_labels = NULL, metadata = list())

# S3 method for class 'MVBU_ModelList'
print(x, ...)

# S3 method for class 'MVBU_ModelList'
summary(object, ...)
```

## Arguments

- models:

  A list of S7 cognitive model objects or `MVBU_StanfitPosterior`
  objects. Default: [`list()`](https://rdrr.io/r/base/list.html).

- model_labels:

  Character vector of labels corresponding to each model. Default:
  `character(0)`.

- metadata:

  Optional list of metadata attributes. Default:
  [`list()`](https://rdrr.io/r/base/list.html).

- x, object:

  An `MVBU_ModelList` object.

- ...:

  Additional arguments passed to methods.

## Examples

``` r
m1 <- example_uvg_ideal_observer()
m2 <- example_uvg_ideal_observer()
ml <- as_model_list(list(m1, m2), model_labels = c("Step 0", "Step 1"))
print(ml)
#> <MVBU_ModelList>
summary(ml)
#> MVBeliefUpdatr Model List containing 2 model(s):
#>   [1] Step 0          : MVBeliefUpdatr::UVG_IdealObserver
#>   [2] Step 1          : MVBeliefUpdatr::UVG_IdealObserver
#> 
#> Model Moments Summary across Steps:
#>    step category expected_mu marginal_mu
#>  Step 0      /b/          NA          NA
#>  Step 0      /d/          NA          NA
#>  Step 0      /g/          NA          NA
#>  Step 0      /k/          NA          NA
#>  Step 0      /p/          NA          NA
#>  Step 0      /t/          NA          NA
#>  Step 1      /b/          NA          NA
#>  Step 1      /d/          NA          NA
#>  Step 1      /g/          NA          NA
#>  Step 1      /k/          NA          NA
#>  Step 1      /p/          NA          NA
#>  Step 1      /t/          NA          NA
```
