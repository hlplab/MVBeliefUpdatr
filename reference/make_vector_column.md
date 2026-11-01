# Deprecated: Combine columns into a vector column

**\[deprecated\]** `make_vector_column()` is deprecated as of
MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0.

## Usage

``` r
make_vector_column(data, cols, vector_col, .keep = "all")
```

## Arguments

- data:

  A tibble or data.frame.

- cols:

  A character vector of variable names to combine.

- vector_col:

  Name of the new vector-valued column.

- .keep:

  A tidyselect option passed to dplyr::mutate.

## Value

Same as `data`.
