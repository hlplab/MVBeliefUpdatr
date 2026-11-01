# Plot categorization decision functions and response boundaries

Visualizes predicted categorization response probabilities and decision
boundaries across cue dimensions for cognitive models and fitted Stan
models. Supports proportional probability matching, deterministic
criterion/MAP classification, and stochastic sampling decision rules,
fully accounting for perceptual noise and lapses.

## Usage

``` r
plot_categorization_functions(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  decision_rule = "proportional",
  resolution = 100,
  noise_treatment = NULL,
  lapse_treatment = NULL,
  ...
)

## S7 method for class <MVBU_CognitiveModel>
plot_categorization_functions(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  decision_rule = "proportional",
  resolution = 100,
  noise_treatment = NULL,
  lapse_treatment = NULL,
  ...
)

## S7 method for class <MVBU_Stanfit>
plot_categorization_functions(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  decision_rule = "proportional",
  resolution = 100,
  noise_treatment = NULL,
  lapse_treatment = NULL,
  ...
)

## S7 method for class <MVBU_ModelList>
plot_categorization_functions(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  decision_rule = "proportional",
  resolution = 100,
  noise_treatment = NULL,
  lapse_treatment = NULL,
  ...
)

## S7 method for class <MVBU_StanfitPosterior>
plot_categorization_functions(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  decision_rule = "proportional",
  resolution = 100,
  noise_treatment = NULL,
  lapse_treatment = NULL,
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

  Character vector of category names to evaluate and plot. Defaults to
  the first category (`get_category_labels(x)[1L]`).

- groups:

  Character vector of group names to plot. Defaults to `NULL` (all
  groups).

- aes:

  Plot aesthetic: `"contour"`, `"fill"`, `"fill-gradient"`, or
  `"fill-discrete"`. Defaults to `"contour"` for 1D/2D static plots, and
  `"fill-discrete"` for 2D interactive plots (`interactive = TRUE`).
  When `interactive = TRUE` (or `aes = "interactive"`), renders a 3D
  response probability surface via Plotly.

- levels:

  Numeric vector specifying response probability contour levels.
  Defaults to `c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99)`.

- limits:

  Optional named list or numeric vector specifying axis / cue limits
  (e.g. `list(F1 = c(200, 800), F2 = c(800, 2500))`). Defaults to `NULL`
  (automatically computed from category distributions or exemplars).

- decision_rule:

  Character string specifying the decision rule to use for
  categorization predictions (`"proportional"`, `"criterion"`, or
  `"sampling"`). Defaults to `"proportional"`.

- resolution:

  Integer specifying grid resolution along continuous cue dimensions.
  Defaults to `100L` for 1D and `60L` for 2D/3D.

- noise_treatment:

  Treatment of perceptual noise during evaluation: `"no_noise"`,
  `"sample"`, or `"marginalize"`. Defaults to `NULL` (uses the model's
  configured noise treatment).

- lapse_treatment:

  Optional character string specifying treatment of attentional lapses:
  `"no_lapses"`, `"sample"`, or `"marginalize"`.

  - `"no_lapses"`: Lapses are ignored; categorization reflects pure
    normative evidence evaluation (even if `lapse_rate > 0`).

  - `"sample"`: Simulates stochastic attentional lapses on individual
    observations. With probability `lapse_rate`, a lapse occurs and the
    posterior is replaced by `lapse_bias`.

  - `"marginalize"`: Evaluates categorization in the limit across
    infinitely many encounters (expected response probabilities), mixing
    normative posteriors with `lapse_bias`: \\(1 - \lambda) p(c
    \mid x) + \lambda b_c\\.

  Defaults to `NULL` (uses the model's configured lapse treatment).

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

## See also

[`plot_categories`](https://hlplab.github.io/MVBeliefUpdatr/reference/plot_categories.md),
[`plot_parameters`](https://hlplab.github.io/MVBeliefUpdatr/reference/plot_parameters.md),
[`categorize`](https://hlplab.github.io/MVBeliefUpdatr/reference/categorize.md),
[`posterior`](https://hlplab.github.io/MVBeliefUpdatr/reference/posterior.md)

## Examples

``` r
mod <- example_mvg_ideal_observer(n_cues = 2)
plot_categorization_functions(mod)
```
