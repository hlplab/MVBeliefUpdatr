# MVBeliefUpdatr Stan Input Classes and Constructors

MVBeliefUpdatr Stan Input Classes and Constructors

## Usage

``` r
MVBU_Staninput(values = list(), family = character(0))

IdealAdaptorStaninput(values = list(), family = character(0))
```

## Arguments

- values:

  A named list containing the Stan input values. Default:
  [`list()`](https://rdrr.io/r/base/list.html).

- family:

  Character string indicating the model family name. Default:
  `character(0)`.

## Examples

``` r
# Instantiate an MVBU_Staninput container object
inp <- MVBU_Staninput(values = list(N = 10L), family = "NIW")
is(inp, "MVBU_Staninput")
#> [1] TRUE
```
