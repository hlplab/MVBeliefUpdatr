# Get the registration information for a model family

Get the registration information for a model family

## Usage

``` r
get_model_family_registration(family)
```

## Arguments

- family:

  Model family name (character string).

## Value

A list with the registered class names and family object.

## See also

[`get_model_family()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_model_family.md),
[`list_model_families()`](https://hlplab.github.io/MVBeliefUpdatr/reference/list_model_families.md)

## Examples

``` r
reg <- get_model_family_registration("UVG")
reg$names$name_short
#> [1] "UVG"
```
