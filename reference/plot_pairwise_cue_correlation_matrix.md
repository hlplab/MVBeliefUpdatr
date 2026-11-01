# Deprecated: plot_pairwise_cue_correlation_matrix

**\[deprecated\]** `plot_pairwise_cue_correlation_matrix()` was
deprecated in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please
use
[`plot_parameters()`](https://hlplab.github.io/MVBeliefUpdatr/reference/plot_parameters.md)
instead.

## Usage

``` r
plot_pairwise_cue_correlation_matrix(
  data,
  cues,
  category = category,
  category.colors = seq_along(unique(data$category))
)
```

## Arguments

- data:

  Data frame.

- cues:

  Cue columns.

- category:

  Category column.

- category.colors:

  Category colors.

## Value

A `ggplot` object.

## See also

[`plot_parameters()`](https://hlplab.github.io/MVBeliefUpdatr/reference/plot_parameters.md)
