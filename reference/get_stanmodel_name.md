# Get the Stan model name associated with a model family

Get the Stan model name associated with a model family

## Usage

``` r
get_stanmodel_name(family)
```

## Arguments

- family:

  Model family short name or
  [MVBU_Family](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Family.md)
  object.

## Value

Stan model name (character string) or `NULL` if the family has no Stan
model.

## See also

[`get_model_family()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_model_family.md),
[`list_model_families()`](https://hlplab.github.io/MVBeliefUpdatr/reference/list_model_families.md)

## Examples

``` r
get_stanmodel_name("NIW")
#> [1] "NIW_ideal_adaptor"
get_stanmodel_name("UVG")
#> NULL
```
