# Get transform metadata from an S7 fit or fit-input object

Get transform metadata from an S7 fit or fit-input object

## Usage

``` r
get_transform_information(x, ...)

## S7 method for class <any>
get_transform_information(x, ...)

## S7 method for class <MVBU_Stanfit>/<MVBU_StanfitPosterior>/<IdealAdaptorStanfitInput>
get_transform_information(x, ...)
```

## Arguments

- x:

  An S7 fit or fit-input object.

- ...:

  Additional arguments passed to methods.

## Value

A transform-information object.

## Examples

``` r
sfi <- example_ideal_adaptor_stanfit_input()
get_transform_information(sfi)
#> <MVBU_TransformInformation>
```
