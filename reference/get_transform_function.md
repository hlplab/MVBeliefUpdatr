# Get the transform/untransform function from an object

Get the transform/untransform function from an object

## Usage

``` r
get_transform_function(x, ...)

get_untransform_function(x, ...)

## S7 method for class <any>
get_transform_function(x, ...)

## S7 method for class <any>
get_untransform_function(x, ...)
```

## Arguments

- x:

  An S7 fit, fit-input, or transform-information object.

- ...:

  Additional arguments.

## Value

A function.

## Examples

``` r
sfi <- example_ideal_adaptor_stanfit_input()
tf <- get_transform_function(sfi)
utf <- get_untransform_function(sfi)
```
