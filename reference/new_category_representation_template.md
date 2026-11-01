# Construct an MVBU Category Representation Template

Construct an MVBU Category Representation Template

## Usage

``` r
new_category_representation_template(representations, metadata = list())
```

## Arguments

- representations:

  A named or unnamed list of `MVBU_CategoryRepresentation` objects.

- metadata:

  Optional metadata list for template-level properties. Defaults to
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

An object of class `MVBU_CategoryRepresentationTemplate`.

## Examples

``` r
rep_b <- new_uvg_category_representation(
  category_labels = "/b/", cue_labels = "VOT",
  mu = 10, sigma2 = 25
)
rep_d <- new_uvg_category_representation(
  category_labels = "/d/", cue_labels = "VOT",
  mu = 50, sigma2 = 25
)
tmpl <- new_category_representation_template(list(rep_b, rep_d))
tmpl
#> <MVBU_CategoryRepresentationTemplate>
#>   Categories (2): /b/, /d/
#>   Cues (1): VOT
#>   Representations (2):
#>     /b/: UVG(mu = 10, sigma2 = 25)
#>     /d/: UVG(mu = 50, sigma2 = 25)
```
