# List all registered model families

Returns the short names of all currently registered model families in
`MVBeliefUpdatr`.

## Usage

``` r
list_model_families()
```

## Value

Character vector of registered family short names.

## See also

[`get_model_family()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_model_family.md),
[`register_model_family()`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Family.md)

## Examples

``` r
list_model_families()
#>  [1] "EXEMPLAR"                               
#>  [2] "MNIX"                                   
#>  [3] "MNIX_category_weighted_integrator"      
#>  [4] "MNIX_group_category_weighted_integrator"
#>  [5] "MNIX_group_weighted_integrator"         
#>  [6] "MNIX_ideal_integrator"                  
#>  [7] "MNIX_informativity_weighted_integrator" 
#>  [8] "MNIX_weighted_integrator"               
#>  [9] "MUVG"                                   
#> [10] "MVG"                                    
#> [11] "NIW"                                    
#> [12] "NIX"                                    
#> [13] "UVG"                                    
#> [14] "category_weighted_integrator_MNIX"      
#> [15] "group_category_weighted_integrator_MNIX"
#> [16] "group_weighted_integrator_MNIX"         
#> [17] "ideal_integrator_MNIX"                  
#> [18] "informativity_weighted_integrator_MNIX" 
#> [19] "weighted_integrator_MNIX"               
```
