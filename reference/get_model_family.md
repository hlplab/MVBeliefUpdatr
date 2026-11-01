# Get the model family from an object

Get the model family from an object

## Usage

``` r
get_model_family(x, ...)

get_family(x, ...)

get_family_name(x, ...)

## S7 method for class <IdealAdaptorStaninput>
get_family_name(x, ...)

## S7 method for class <IdealAdaptorStaninput>
get_family(x, ...)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>
get_model_family(x, ...)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>
get_family_name(x, ...)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>
get_family(x, ...)
```

## Arguments

- x:

  An object.

- ...:

  Additional arguments passed to methods.

## Value

A character string representing the model family.

## See also

`get_family`, `get_family_name`

## Examples

``` r
mod <- example_mvg_ideal_observer(n_cues = 2)
get_model_family(mod)
#> [1] "MVBeliefUpdatr::MVG_IdealObserver"
get_family_name(mod)
#> [1] "MVG"
get_family(mod)
#> <MVBeliefUpdatr::MVBU_Family>
#>  @ family_names                  :List of 3
#>  .. $ name_short  : chr "MVG"
#>  .. $ name_long   : chr "Multivariate Gaussian Ideal Observer"
#>  .. $ name_details:List of 4
#>  ..  ..$ CategoryRepresentation_name: chr "MVG_CategoryRepresentation"
#>  ..  ..$ CategoryTemplate_name      : chr "MVG_CategoryRepresentationTemplate"
#>  ..  ..$ CognitiveModel_name        : chr "MVG_IdealObserver"
#>  ..  ..$ CognitiveModelList_name    : chr "MVG_CognitiveModelList"
#>  @ density_function              : function (x, rep, log = FALSE, noise_treatment = "no_noise", Sigma_noise = NULL)  
#>  @ sampling_function             : function (rep, n = 1L, with_replacement = TRUE, ...)  
#>  @ expected_moments              : function (rep)  
#>  @ marginal_moments              : NULL
#>  @ sufficient_statistics_function: NULL
#>  @ stanmodel                     : NULL
```
