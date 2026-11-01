# Plot model parameter estimates, correlations, and pairwise posterior distributions

Visualizes category parameter estimates (means, covariance ellipses,
SDs, correlations), parameter correlation matrices across MCMC draws,
and pairwise posterior parameter scatter/contour matrices for category
representations, templates, cognitive models, and fitted Stan models.

## Usage

``` r
plot_parameters(
  x,
  pars = NULL,
  categories = NULL,
  groups = NULL,
  ndraws = 100,
  combine_into_single_plot = TRUE,
  index_panels = FALSE,
  ...
)

plot_parameter_correlations(
  x,
  pars = NULL,
  categories = NULL,
  groups = NULL,
  cues = NULL,
  ndraws = 100,
  ...
)

plot_parameters_pairwise(
  x,
  pars = NULL,
  categories = NULL,
  groups = NULL,
  cues = NULL,
  ndraws = 100,
  ...
)

## S7 method for class <MVBU_CognitiveModel>
plot_parameters(
  x,
  pars = NULL,
  categories = NULL,
  groups = NULL,
  ndraws = 100,
  combine_into_single_plot = TRUE,
  index_panels = FALSE,
  ...
)

## S7 method for class <MVBU_Stanfit>
plot_parameters(
  x,
  pars = NULL,
  categories = NULL,
  groups = NULL,
  ndraws = 100,
  combine_into_single_plot = TRUE,
  index_panels = FALSE,
  ...
)

## S7 method for class <MVBU_Stanfit>
plot_parameter_correlations(
  x,
  pars = NULL,
  categories = NULL,
  groups = NULL,
  cues = NULL,
  ndraws = 100,
  ...
)

## S7 method for class <MVBU_Stanfit>
plot_parameters_pairwise(
  x,
  pars = NULL,
  categories = NULL,
  groups = NULL,
  cues = NULL,
  ndraws = 100,
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

- pars:

  Character vector of parameter names to filter in `plot_parameters`,
  `plot_parameter_correlations`, and `plot_parameters_pairwise`.
  Defaults to `NULL` (all parameters).

- categories:

  Character vector of category names to include. Defaults to `NULL` (all
  categories).

- groups:

  Character vector of group names to plot for
  [`MVBU_Stanfit`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Stanfit.md)
  objects. For `plot_parameters`, defaults to all groups including the
  prior. For `plot_parameters_pairwise` and
  `plot_parameter_correlations`, defaults to `"prior"`.

- ndraws:

  Integer specifying the number of posterior MCMC draws to sample for
  plotting parameter distributions. Defaults to `100L`.

- combine_into_single_plot:

  Logical; if `TRUE` (default), combines multi-panel parameter subplots
  into a single unified figure via patchwork. If `FALSE`, returns a
  named list of individual `ggplot` panel objects.

- index_panels:

  Logical; if `TRUE`, adds panel index labels (A, B, ...) to multi-panel
  figures. Defaults to `FALSE`.

- ...:

  Additional arguments passed to specific parameter plotting methods.

- cues:

  Character vector of cue names to plot (1, 2, or 3 cues). Defaults to
  all cues of the model up to 3 (or the first 3 cues if the model has
  more than 3 cues, analytically marginalizing out all remaining
  dimensions). At most 3 cue dimensions can be plotted simultaneously.

## Value

A `ggplot` object (or list of `ggplot` objects if
`combine_into_single_plot = FALSE`).

## Details

In parameter plots, scale is displayed as standard deviation \\\tau =
\text{SD} = \sqrt{\text{diag}(\boldsymbol{\Sigma})}\\ (or the expected
standard deviation \\\sqrt{\text{diag}(\mathbf{S}) / (\nu - D - 1)}\\
for NIW/NIX models), while off-diagonal covariances are plotted as
dimensionless correlation coefficients \\\rho \in \[-1, 1\]\\.

`plot_parameter_correlations` displays a correlation heatmap of Pearson
correlation coefficients (\\\rho\\) computed across MCMC draws for all
parameters across categories and selected groups. Parameters for
distinct groups form separate rows and columns without averaging across
groups. Parameter names indicate category, cue, and group index in their
subscripts (with subscript `0` denoting the prior group). Gray bounding
boxes along the diagonal group parameters belonging to each category.

`plot_parameters_pairwise` displays a matrix of pairwise parameter
relationships across categories and groups. The lower triangle shows
sample scatterplots with LOESS trendlines; the diagonal displays
marginal posterior auto-densities scaled to panel height; and the upper
triangle displays 2D posterior density contour lines (`geom_density_2d`
with normalized density contours). Points and curves are colored by
group.

## See also

[`plot_categories`](https://hlplab.github.io/MVBeliefUpdatr/reference/plot_categories.md),
[`plot_categorization_functions`](https://hlplab.github.io/MVBeliefUpdatr/reference/plot_categorization_functions.md),
[`get_parameters`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_parameters.md),
[`get_draws`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_draws.md)

## Examples

``` r
mod <- example_mvg_ideal_observer(n_cues = 2)
plot_parameters(mod)
```
