# Multi-Cue Normal-Inverse-Chi-Squared (MNIX) Family: Representations, Templates, and Models

Constructs and manages multi-cue Normal-Inverse-\\\chi^2\\ (MNIX)
category representations, templates, and ideal adaptors under early and
late cue integration.

## Usage

``` r
new_mnix_category_representation(
  category_labels,
  cue_labels,
  m,
  kappa,
  nu,
  sigma2,
  weights = NULL,
  metadata = list()
)

new_mnix_category_representation_from_data(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = 3,
  ...
)

new_mnix_category_representation_template_from_data(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = 3,
  ...
)

new_mnix_ideal_adaptor_ideal_integrator(
  category_template = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
  metadata = list()
)

new_mnix_ideal_adaptor_weighted_integrator(
  category_template = NULL,
  weights = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
  metadata = list()
)

new_mnix_ideal_adaptor_category_weighted_integrator(
  category_template = NULL,
  weights = NULL,
  weights_matrix = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
  metadata = list()
)

new_mnix_ideal_adaptor_informativity_weighted_integrator(
  category_template = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
  metadata = list()
)

new_mnix_ideal_adaptor(
  category_template = NULL,
  weights = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
  metadata = list()
)

new_mnix_ideal_adaptor_from_data(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = 3,
  weights = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
  ...
)

new_ideal_integrator_mnix_ideal_adaptor(
  category_template = NULL,
  cues = NULL,
  input_cues = NULL,
  weights = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
  metadata = list(),
  ...
)

new_ideal_integrator_mnix_ideal_adaptor_from_data(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = 3,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate)
)

new_weighted_integrator_mnix_ideal_adaptor(
  category_template = NULL,
  weights = NULL,
  cues = NULL,
  input_cues = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
  metadata = list()
)

new_weighted_integrator_mnix_ideal_adaptor_from_data(
  data,
  category = "category",
  cues,
  weights = NULL,
  kappa = nu,
  nu = 3,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate)
)

new_category_weighted_integrator_mnix_ideal_adaptor(
  category_template = NULL,
  weights = NULL,
  weights_matrix = NULL,
  cues = NULL,
  input_cues = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
  metadata = list()
)

new_informativity_weighted_integrator_mnix_ideal_adaptor(
  category_template = NULL,
  cues = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
  metadata = list()
)

new_mnix_ideal_adaptor_adaptive_weighted_integrator(...)

new_mnix_ideal_adaptor_category_adaptive_weighted_integrator(...)
```

## Arguments

- metadata:

  Optional list of metadata. Defaults to
  [`list()`](https://rdrr.io/r/base/list.html).

- m:

  Numeric vector of per-cue location hyperparameters.

- kappa:

  Positive numeric scalar or vector of precision scalings (\\\kappa \>
  0\\). For late integration, a vector of length equal to cues (or
  scalar recycled across cues). For early integration, a scalar
  reflecting composite dimension belief strength.

- nu:

  Positive numeric scalar or vector of degrees of freedom (\\\nu \>
  2\\). For late integration, a vector of length equal to cues (or
  scalar recycled across cues). For early integration, a scalar
  reflecting composite dimension degrees of freedom.

- sigma2:

  Numeric vector of per-cue scale parameters.

- weights:

  Numeric vector or matrix of integration/attention weights. For
  representations and global weighted models, a numeric simplex vector
  summing to 1. For category-weighted models, an \\M \times K\\ matrix
  whose rows each sum to 1.

- category_labels:

  Character vector of category label(s).

- cue_labels:

  Character vector of cue labels (at least 2 cue dimensions).

- data:

  Data frame containing category and cue observations.

- category:

  Character string indicating the category column in `data`. Defaults to
  `"category"`.

- cues:

  Character vector indicating cue columns in `data` (at least 2 cues).

- ...:

  Additional arguments passed to methods or constructors.

- category_template:

  An
  [MVBU_CategoryRepresentationTemplate](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md)
  object containing MNIX category representations.

- decision_rule:

  Categorization decision rule: `"sampling"`, `"criterion"`, or
  `"proportional"`. Defaults to `"sampling"`.

- category_prior:

  Optional numeric vector of prior category probabilities summing to 1.
  Defaults to equal priors.

- lapse_rate:

  Numeric scalar lapse probability in `[0, 1]`. Defaults to 0.

- lapse_bias:

  Optional numeric vector of lapse category probabilities summing to 1.
  Defaults to uniform bias `1/M`.

- Sigma_noise:

  Optional perceptual/measurement noise covariance matrix or vector.

- noise_treatment:

  Treatment of perceptual noise during evaluation and updating:
  `"no_noise"`, `"sample"`, or `"marginalize"`.

  - `"no_noise"`: No perceptual noise is applied to input cues, and no
    noise is assumed during categorization.

  - `"sample"`: Simulates sensory noise on a specific observation. Noise
    drawn from \\\mathcal{N}(0, \Sigma\_{\text{noise}})\\ is added to
    input cues (and projected into composite cue space in early
    integration), and the observer accounts for their own noise
    uncertainty during categorization by adding
    \\\Sigma\_{\text{noise}}\\ to their likelihood scale.

  - `"marginalize"`: Evaluates categorization in the limit across
    infinitely many noisy encounters (expected consequences). No noise
    is sampled onto the input stimulus, but \\\Sigma\_{\text{noise}}\\
    is added to the category covariance scale.

  Defaults to `.infer_noise_treatment(Sigma_noise)`: `"marginalize"` if
  noise is present, `"no_noise"` otherwise.

- lapse_treatment:

  Treatment of attentional lapses during evaluation and updating:
  `"no_lapses"`, `"sample"`, or `"marginalize"`.

  - `"no_lapses"`: Lapses are ignored; categorization and updating
    reflect pure normative evidence evaluation (even if
    `lapse_rate > 0`).

  - `"sample"`: Simulates stochastic attentional lapses on individual
    observations. With probability `lapse_rate`, a lapse occurs: during
    categorization, the posterior is replaced by `lapse_bias`; during
    updating, the lapsing observation is treated as uninformative.

  - `"marginalize"`: Evaluates categorization in the limit across
    infinitely many encounters (expected response probabilities). The
    category posterior is the mixture \\(1 - \lambda) p(c \mid x) +
    \lambda b_c\\, and evidence during updating is downweighted by \\1 -
    \lambda\\.

  Defaults to `.infer_lapse_treatment(lapse_rate)`: `"marginalize"` if
  `lapse_rate > 0`, `"no_lapses"` otherwise.

- weights_matrix:

  Deprecated. Use `weights` instead.

- input_cues:

  Deprecated. Use `cues` instead.

## Value

An S7 object of the respective MNIX class.

## Details

The MNIX family unifies multi-cue belief updating under both late and
early integration regimes:

- **Late Integration**: Categories are represented by independent
  marginal Student-\\t\\ predictive components along each cue dimension
  (\\\kappa\\ and \\\nu\\ have length \\K\\). Cue log-likelihoods are
  combined at decision time via Bayes rule or attention weights.

- **Early Integration**: Multi-cue observations are projected to a
  composite dimension \\z = \mathbf{w}^\top \mathbf{x}\\ prior to
  updating (\\\kappa\\ and \\\nu\\ are scalar length 1).

## See also

[family-muvg](https://hlplab.github.io/MVBeliefUpdatr/reference/family-muvg.md),
[family-niw](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md),
[family-nix](https://hlplab.github.io/MVBeliefUpdatr/reference/family-nix.md)

## Examples

``` r
# Mixture / Multivariate NIX representation
mu_b <- c(VOT = 10, f0 = 200)
sigma2_b <- c(VOT = 25, f0 = 36)
rep_b <- new_mnix_category_representation(
  category_labels = "/b/", cue_labels = c("VOT", "f0"),
  m = mu_b, kappa = c(5, 5), nu = c(8, 8),
  sigma2 = sigma2_b
)
mu_d <- c(VOT = 50, f0 = 240)
sigma2_d <- c(VOT = 25, f0 = 36)
rep_d <- new_mnix_category_representation(
  category_labels = "/d/", cue_labels = c("VOT", "f0"),
  m = mu_d, kappa = c(5, 5), nu = c(8, 8),
  sigma2 = sigma2_d
)
tmpl <- new_category_representation_template(list(rep_b, rep_d))
m_mnix <- new_mnix_ideal_adaptor(category_template = tmpl)
m_mnix
#> <MVBeliefUpdatr::MNIX_IdealAdaptor_IdealIntegrator>
#>   Categories (2): /b/, /d/
#>   Cues (2): VOT, f0
#>   Decision rule: sampling
#>   Category prior: /b/: 0.5, /d/: 0.5
#>   Lapse rate: 0 (treatment: no_lapses)
#>   Perceptual noise: none (treatment: no_noise)
#>   Category representations (2):
#>     /b/: MNIX(kappa = vector(2), nu = vector(2), m = vector(2), sigma2 = vector(2))
#>       E[mu] = vector(2), E[Sigma] = matrix(2, 2), Sigma_marg = matrix(2, 2)
#>     /d/: MNIX(kappa = vector(2), nu = vector(2), m = vector(2), sigma2 = vector(2))
#>       E[mu] = vector(2), E[Sigma] = matrix(2, 2), Sigma_marg = matrix(2, 2)
```
