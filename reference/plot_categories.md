# Plot category representations and distributions

Generic plotting function for visualizing category representations
(densities, contour lines, confidence ellipses, exemplar clouds, and 3D
density surfaces/ellipsoids) across single category representations,
multi-category templates, cognitive decision models, and fitted Stan
models.

## Usage

``` r
plot_categories(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  resolution = 100,
  n_exemplars = 0L,
  ...
)

## S7 method for class <MVBU_CategoryRepresentation>
plot_categories(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  resolution = 100,
  n_exemplars = 0L,
  ...
)

## S7 method for class <MVBU_CategoryRepresentationTemplate>
plot_categories(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  resolution = 100,
  n_exemplars = 0L,
  ...
)

## S7 method for class <MVBU_CognitiveModel>
plot_categories(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  resolution = 100,
  n_exemplars = 0L,
  ...
)

## S7 method for class <MVBU_Stanfit>
plot_categories(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  resolution = 100,
  n_exemplars = 0L,
  ...
)

## S7 method for class <MVBU_ModelList>
plot_categories(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  resolution = 100,
  n_exemplars = 0L,
  ...
)

## S7 method for class <MVBU_StanfitPosterior>
plot_categories(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  resolution = 100,
  n_exemplars = 0L,
  ...
)
```

## Arguments

- x:

  An MVBU object
  ([`MVBU_CategoryRepresentation`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md),
  [`MVBU_CategoryRepresentationTemplate`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md),
  [`MVBU_CognitiveModel`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md),
  or
  [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)).

- cues:

  Character vector of cue names to plot (1, 2, or 3 cues). Defaults to
  all cues of the model up to 3 (or the first 3 cues if the model has
  more than 3 cues, analytically marginalizing out all remaining
  dimensions). At most 3 cue dimensions can be plotted simultaneously.

- categories:

  Character vector of category names to include. Defaults to `NULL` (all
  categories).

- groups:

  Character vector of group names to plot. Defaults to `NULL` (all
  groups).

- aes:

  Plot aesthetic: `"contour"`, `"fill"`, `"fill-gradient"`,
  `"fill-discrete"`, `"scatter"`, or combinations like
  `c("fill-gradient", "contour")`. For 1D and 2D categories, defaults to
  `c("fill-gradient", "contour")`. For 3D slices, defaults to
  `"contour"`. For 3D interactive parametric representations, defaults
  to `"fill-discrete"`. For exemplar representations, defaults to
  `"scatter"`.

- levels:

  Numeric vector or list specifying density / probability levels.
  Defaults to two-tailed central probabilities corresponding to 1, 2, 3
  sigma (`2 * stats::pnorm(1:3) - 1`) for 1D/2D and sliced plots, and 2
  sigma (`2 * stats::pnorm(2) - 1`) for 3D interactive ellipsoids. For a
  bivariate normal distribution \\\mathcal{N}(\boldsymbol{\mu},
  \boldsymbol{\Sigma})\\, the squared Mahalanobis distance follows
  \\\chi^2_2\\, meaning the ellipse enclosing mass \\p \in (0, 1)\\ has
  radius \\r = \sqrt{F\_{\chi^2_2}^{-1}(p)} = \sqrt{-2 \ln(1 - p)}\\.

- limits:

  Optional named list or numeric vector specifying axis / cue limits
  (e.g. `list(F1 = c(200, 800), F2 = c(800, 2500))`). Defaults to `NULL`
  (automatically computed from category distributions or exemplars).

- resolution:

  Integer specifying the grid resolution along continuous cue
  dimensions. Defaults to `100L` for 1D/2D category plots and `60L` for
  3D sliced category plots.

- n_exemplars:

  Integer specifying max number of exemplars to randomly sample and plot
  for exemplar representations. Defaults to `0L` for 1D/2D

- ...:

  Additional arguments passed to specific plotting methods:

  - `uncertainty_treatment`: Treatment of posterior parameter
    uncertainty in fitted Stan models (`"marginalize"` or `"sample"`).

  - `ndraws`: Integer specifying the number of posterior MCMC draws to
    use from
    [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
    (defaults to `20L`).

  - `show_exposure_data`, `show_test_data`: Logical; if `TRUE`, overlays
    observed exposure/test data points.

  - `parallel`, `n_cores`: Logical and integer controlling parallel grid
    evaluation.

  - `noise_treatment`: Optional character string specifying perceptual
    noise treatment (`"no_noise"`, `"sample"`, or `"marginalize"`).

  - `interactive`: Logical; if `TRUE`, produces interactive 3D WebGL
    scenes via Plotly.

  - `slices`, `slice_cue`: Configuration for 3D sliced category plots.

## Value

A `ggplot` object (or a `plotly` htmlwidget if `interactive = TRUE`).
Because ggplot output is a standard `ggplot` / patchwork object, it can
be further customized with additional layers, scales, themes, or facets.

## Details

In 3D interactive category plots (`plot_categories` with
`interactive = TRUE`), large diamond markers represent category means.
For exemplar representations, the mean is computed across all stored
exemplars.

## See also

[`plot_categorization_functions`](https://hlplab.github.io/MVBeliefUpdatr/reference/plot_categorization_functions.md),
[`plot_parameters`](https://hlplab.github.io/MVBeliefUpdatr/reference/plot_parameters.md),
[`likelihood`](https://hlplab.github.io/MVBeliefUpdatr/reference/likelihood.md),
[`posterior`](https://hlplab.github.io/MVBeliefUpdatr/reference/posterior.md)

## Examples

``` r
mod <- example_mvg_ideal_observer(n_cues = 2)
plot_categories(mod)
```
