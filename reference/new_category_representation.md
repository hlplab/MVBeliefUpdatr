# Construct a base category representation object

Construct a base category representation object

## Usage

``` r
new_category_representation(
  category_labels,
  cue_labels,
  category_likelihood_function = NULL,
  parameters = list(),
  metadata = list()
)
```

## Arguments

- category_labels:

  Character vector of category label(s).

- cue_labels:

  Character vector of cue label(s).

- category_likelihood_function:

  Optional category likelihood function. Default: `NULL`.

- parameters:

  Optional named list of representation parameters. Default:
  [`list()`](https://rdrr.io/r/base/list.html).

- metadata:

  Optional list of metadata. Default:
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

An
[MVBU_CategoryRepresentation](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md)
object.
