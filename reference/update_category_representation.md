# Update a category representation from observations

Update a category representation from observations

## Usage

``` r
update_category_representation(x, n_exposure, mean_exposure, css_exposure, ...)

## S7 method for classes <MVBeliefUpdatr::NIW_CategoryRepresentation>, <any>, <any>, <any>
update_category_representation(x, n_exposure, mean_exposure, css_exposure, ...)

## S7 method for classes <MVBeliefUpdatr::NIX_CategoryRepresentation>, <any>, <any>, <any>
update_category_representation(x, n_exposure, mean_exposure, css_exposure, ...)

## S7 method for classes <MNIX_CategoryRepresentation>, <any>, <any>, <any>
update_category_representation(x, n_exposure, mean_exposure, css_exposure, ...)
```

## Arguments

- x:

  A category representation object.

- n_exposure:

  Number of observations represented by the sufficient statistics.

- mean_exposure:

  Mean vector of the represented observations.

- css_exposure:

  Centered sum-of-squares matrix or vector of the represented
  observations.

- ...:

  Additional arguments passed to methods.

## Value

An updated category representation.

## Examples

``` r
cat_rep <- example_niw_category_representation(n_cues = 2)
update_category_representation(
  cat_rep,
  n_exposure = 5L,
  mean_exposure = c(0.5, 0.5),
  css_exposure = matrix(c(1, 0.2, 0.2, 1), nrow = 2)
)
#> <MVBeliefUpdatr::NIW_CategoryRepresentation>
#>   Category: /b/
#>   Cues (2): VOT, f0_semitones
#>   kappa: 6
#>   nu: 9
#>   m: c(1.8296, 0.1887)
#>   S:
#>                    VOT f0_semitones
#> VOT           77.23727    -12.30740
#> f0_semitones -12.30740     18.59764
```
