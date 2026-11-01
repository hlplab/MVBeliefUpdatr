# Convert S7 objects to legacy tibble format

Converts modern S7 cognitive models, category representation templates,
or individual category representations into legacy tibbles matching the
schemas from pre-S7 versions of MVBeliefUpdatr.

## Usage

``` r
# S3 method for class 'MVBU_Object'
as_tibble(x, ...)

# S3 method for class 'MVBU_CategoryRepresentation'
as_tibble(x, ...)

# S3 method for class 'MVBU_CategoryRepresentationTemplate'
as_tibble(x, ...)

# S3 method for class 'MVBU_CognitiveModel'
as_tibble(x, ...)
```

## Arguments

- x:

  An S7 object inheriting from
  [`MVBU_Object`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md).

- ...:

  Additional arguments (currently unused).

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
columns and classes adhering to the legacy format.
