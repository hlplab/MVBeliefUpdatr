# Multi-Cue Univariate Gaussian (MUVG) Family: Representations, Templates, and Models

Constructs and manages multi-cue univariate Gaussian integration
representations, templates, and ideal observers.

## Usage

``` r
new_muvg_category_representation(
  category_labels,
  cue_labels,
  mu,
  sigma2,
  weights = NULL,
  metadata = list()
)

new_muvg_category_representation_from_data(data, category = "category", cues)

new_muvg_category_representation_template_from_data(
  data,
  category = "category",
  cues
)

new_muvg_ideal_observer(
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

new_muvg_ideal_observer_from_data(
  data,
  category = "category",
  cues,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate)
)
```

## Arguments

- metadata:

  Optional list of metadata.

- mu:

  Numeric vector of per-cue means.

- sigma2:

  Numeric vector of per-cue positive variances.

- weights:

  Optional numeric vector of integration weights in `[0, 1]` summing
  to 1. If `NULL`, defaults to precision weights \\w_i \propto
  1/\sigma_i^2\\.

- category_labels:

  Character vector of category label(s).

- cue_labels:

  Character vector of cue labels matching the length of cue components.

- data:

  Data frame containing category and cue observations.

- category:

  Character string indicating the category column in `data`. Defaults to
  `"category"`.

- cues:

  Character vector indicating cue columns in `data`.

- category_template:

  An
  [MVBU_CategoryRepresentationTemplate](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md)
  object containing MUVG category representations.

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

- Sigma_noise:

  Optional perceptual/measurement noise covariance matrix or vector.

- noise_treatment:

  Treatment of perceptual noise during evaluation: `"no_noise"`,
  `"sample"`, or `"marginalize"`.

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

  Treatment of attentional lapses during evaluation: `"no_lapses"`,
  `"sample"`, or `"marginalize"`.

  - `"no_lapses"`: Lapses are ignored; categorization reflects pure
    normative evidence evaluation (even if `lapse_rate > 0`).

  - `"sample"`: Simulates stochastic attentional lapses on individual
    observations. With probability `lapse_rate`, a lapse occurs: the
    posterior is replaced by `lapse_bias`.

  - `"marginalize"`: Evaluates categorization in the limit across
    infinitely many encounters (expected response probabilities). The
    category posterior is the mixture \\(1 - \lambda) p(c \mid x) +
    \lambda b_c\\.

  Defaults to `.infer_lapse_treatment(lapse_rate)`: `"marginalize"` if
  `lapse_rate > 0`, `"no_lapses"` otherwise.

## Value

An S7 object of the respective MUVG class.

## Details

Under conditional cue independence and integrated cue \\z = \sum_i w_i
x_i\\, the category likelihood is: \$\$p(z \mid c) =
\mathcal{N}\\\left(z;\\ \sum_i w_i \mu_i,\\ \sum_i w_i^2
\sigma_i^2\right)\$\$

## See also

[family-mnix](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mnix.md),
[family-mvg](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md),
[family-uvg](https://hlplab.github.io/MVBeliefUpdatr/reference/family-uvg.md)

## Examples

``` r
# MUVG: Multivariate independent Gaussian (diagonal covariance)
mu_b <- c(VOT = 10, f0 = 200)
sigma2_b <- c(VOT = 25, f0 = 36)
rep_b <- new_muvg_category_representation(
  category_labels = "/b/", cue_labels = c("VOT", "f0"),
  mu = mu_b, sigma2 = sigma2_b
)
mu_d <- c(VOT = 50, f0 = 240)
sigma2_d <- c(VOT = 25, f0 = 36)
rep_d <- new_muvg_category_representation(
  category_labels = "/d/", cue_labels = c("VOT", "f0"),
  mu = mu_d, sigma2 = sigma2_d
)
tmpl <- new_category_representation_template(list(rep_b, rep_d))
m_muvg <- new_muvg_ideal_observer(category_template = tmpl)
m_muvg
#> <MVBeliefUpdatr::MUVG_IdealObserver>
#>   Categories (2): /b/, /d/
#>   Cues (2): VOT, f0
#>   Decision rule: sampling
#>   Category prior: /b/: 0.5, /d/: 0.5
#>   Lapse rate: 0 (treatment: no_lapses)
#>   Perceptual noise: none (treatment: no_noise)
#>   Category representations (2):
#>     /b/: MUVG(mu = vector(2), sigma2 = vector(2))
#>     /d/: MUVG(mu = vector(2), sigma2 = vector(2))
```
