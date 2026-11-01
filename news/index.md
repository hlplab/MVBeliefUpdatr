# Changelog

## MVBeliefUpdatr 0.1.0

Initial release of **MVBeliefUpdatr**, a package for modeling, fitting,
and visualizing incremental Bayesian categorization and belief updating
for Gaussian ideal observers and ideal adaptors.

### Key Features

#### S7 Class Architecture

- Unified hierarchy for cognitive models based on R’s modern ‘S7’ class
  system.
- Separation of concerns across:
  - `MVBU_CategoryRepresentation`: category distribution representations
    for univariate Gaussian (`UVG`), multivariate Gaussian (`MVG`), and
    univariate Gaussian with category-specific variance (`MUVG`).
  - `MVBU_CategoryTemplate`: collections of category representations
    representing an observer’s prior or updated beliefs.
  - `MVBU_CognitiveModel`: ideal observers, ideal adaptors, and
    exemplar-based categorization models.

#### Bayesian Belief Updating & Conjugate Hyperpriors

- Analytical closed-form Bayesian belief updating under conjugate
  hyperprior distributions:
  - Normal-Inverse-Chi-squared (`NIX`) for 1D cues.
  - Normal-Inverse-Wishart (`NIW`) for multivariate cues.
  - Multivariate Normal-Inverse-Chi-squared (`MNIX`) for multivariate
    cues with diagonal covariance matrices.
- Support for four moment-updating regimes per category:
  - Unknown prior mean and unknown prior variance/covariance (Case 1).
  - Unknown prior mean and known prior variance/covariance (Case 2).
  - Known prior mean and unknown prior variance/covariance (Case 3).
  - Known prior mean and known prior variance/covariance (Case 4).
- Declarative prior specification with
  [`parameter_prior()`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Prior.md),
  supporting informative Normal, Cauchy, LKJ, Beta, and Dirichlet
  priors.
- Partial and global parameter fixing with
  [`fix_parameter()`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_FixedParameter.md).

#### Stan Integration & MCMC Parameter Inference

- Pre-compiled Stan models providing full Markov Chain Monte
  Carlo (MCMC) parameter estimation for univariate and multivariate
  Gaussian ideal adaptors.
- Integration with ‘loo’ for Leave-One-Out Cross-Validation model
  comparison via [`loo()`](https://mc-stan.org/loo/reference/loo.html),
  [`waic()`](https://mc-stan.org/loo/reference/waic.html), and
  [`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html).
- Tidy extraction methods
  ([`get_parameters()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_parameters.md),
  [`get_draws()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_draws.md),
  `get_predictions()`) supporting both natural cue space and
  standardized space.

#### Visualization Suite

- Rich visualization tools built with ‘ggplot2’ and ‘patchwork’:
  - [`plot_categories()`](https://hlplab.github.io/MVBeliefUpdatr/reference/plot_categories.md):
    1D density profiles, 2D ellipses, and 3D sliced hypersurfaces.
  - [`plot_categorization_functions()`](https://hlplab.github.io/MVBeliefUpdatr/reference/plot_categorization_functions.md):
    posterior categorization functions along continuous cue dimensions.
  - [`plot_parameters()`](https://hlplab.github.io/MVBeliefUpdatr/reference/plot_parameters.md):
    multi-panel posterior distributions and updating trajectories.
