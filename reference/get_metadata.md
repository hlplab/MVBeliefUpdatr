# Get metadata from an object

Get metadata from an object

## Usage

``` r
get_metadata(x, ...)

## S7 method for class <MVBU_ModelList>
get_metadata(x, ...)
```

## Arguments

- x:

  An object.

- ...:

  Additional arguments passed to methods.

## Value

A list containing metadata associated with the object.

## Examples

``` r
mod <- example_mvg_ideal_observer(n_cues = 2)
get_metadata(mod)
#> $label_information
#> $label_information$cue
#> [1] "VOT"          "f0_semitones"
#> 
#> $label_information$category
#> [1] "/b/" "/d/" "/g/" "/k/" "/p/" "/t/"
#> 
#> $label_information$response_category
#> [1] "/b/" "/d/" "/g/" "/k/" "/p/" "/t/"
#> 
#> $label_information$group
#> character(0)
#> 
#> 
#> $parameters
#> $parameters$`/b/`
#> $parameters$`/b/`$mu
#> [1]  8.477551 -1.368047
#> 
#> $parameters$`/b/`$Sigma
#>                      VOT f0_semitones
#> VOT          23.20282907  -0.08870048
#> f0_semitones -0.08870048  14.68964368
#> 
#> 
#> $parameters$`/d/`
#> $parameters$`/d/`$mu
#> [1] 13.917921 -0.580694
#> 
#> $parameters$`/d/`$Sigma
#>                    VOT f0_semitones
#> VOT          75.914604    -1.445236
#> f0_semitones -1.445236    16.963832
#> 
#> 
#> $parameters$`/g/`
#> $parameters$`/g/`$mu
#> [1] 16.9919098 -0.7449837
#> 
#> $parameters$`/g/`$Sigma
#>                      VOT f0_semitones
#> VOT          107.3474463   -0.1125058
#> f0_semitones  -0.1125058   20.0938994
#> 
#> 
#> $parameters$`/k/`
#> $parameters$`/k/`$mu
#> [1] 54.544491  0.576098
#> 
#> $parameters$`/k/`$Sigma
#>                    VOT f0_semitones
#> VOT          409.69309    -11.39129
#> f0_semitones -11.39129     13.89484
#> 
#> 
#> $parameters$`/p/`
#> $parameters$`/p/`$mu
#> [1] 50.6026289 -0.4758416
#> 
#> $parameters$`/p/`$Sigma
#>                     VOT f0_semitones
#> VOT          408.755971    -6.736239
#> f0_semitones  -6.736239    17.966702
#> 
#> 
#> $parameters$`/t/`
#> $parameters$`/t/`$mu
#> [1] 61.2998157 -0.2026473
#> 
#> $parameters$`/t/`$Sigma
#>                     VOT f0_semitones
#> VOT          436.650466    -9.659716
#> f0_semitones  -9.659716    15.159686
#> 
#> 
#> 
```
