# Sample observations from category representations, templates, or cognitive models

Generates simulated cue observations by sampling from an individual
category representation, a multi-category template, or a complete
cognitive model.

- For single representations
  ([`MVBU_CategoryRepresentation`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md)),
  samples \\n\\ observations from the representation's distribution.

- For templates
  ([`MVBU_CategoryRepresentationTemplate`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md)),
  samples \\n\\ total observations distributed uniformly across
  categories.

- For cognitive models
  ([`MVBU_CognitiveModel`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md)),
  samples \\n\\ total observations distributed across categories
  proportionally to the model's `category_prior`.

## Usage

``` r
sample_observations(
  x,
  n = 1L,
  with_replacement = TRUE,
  randomize_order = TRUE,
  ...
)

## S7 method for class <MVBU_Stanfit>
sample_observations(
  x,
  n = 1L,
  with_replacement = TRUE,
  randomize_order = TRUE,
  ...
)
```

## Arguments

- x:

  An
  [`MVBU_CategoryRepresentation`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md),
  [`MVBU_CategoryRepresentationTemplate`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md),
  or
  [`MVBU_CognitiveModel`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md)
  object.

- n:

  Non-negative integer specifying the total number of observations to
  sample. Defaults to `1L`.

- with_replacement:

  Logical; for exemplar representations, specifies whether to sample
  exemplars with replacement (`TRUE`, default) or without replacement
  (`FALSE`). Ignored for parametric families.

- randomize_order:

  Logical; whether to randomize the row order of the returned
  observations. Defaults to `TRUE`.

- ...:

  Additional arguments passed to methods (e.g., `Ns` or
  `randomize.order` for backwards compatibility).

## Value

A tibble with \\n\\ rows containing a `category` factor column (with
levels matching category labels) and one numeric column per cue.

## See also

[`plot_categories`](https://hlplab.github.io/MVBeliefUpdatr/reference/plot_categories.md),
[`likelihood`](https://hlplab.github.io/MVBeliefUpdatr/reference/likelihood.md),
[`posterior`](https://hlplab.github.io/MVBeliefUpdatr/reference/posterior.md)

## Examples

``` r
mod <- example_mvg_ideal_observer(n_cues = 2)
sample_observations(mod, n = 5)
#> # A tibble: 5 × 3
#>   category   VOT f0_semitones
#>   <fct>    <dbl>        <dbl>
#> 1 /t/      95.0        -0.694
#> 2 /g/       5.54        1.68 
#> 3 /p/      62.5        -3.42 
#> 4 /p/      82.4        -5.93 
#> 5 /t/      64.9        -4.47 
```
