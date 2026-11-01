# Univariate Gaussian (UVG) Family: Representations, Templates, and Models

Constructs and manages univariate Gaussian category representations,
templates, and ideal observer cognitive models.

## Usage

``` r
new_uvg_category_representation(
  category_labels,
  cue_labels,
  mu,
  sigma2,
  metadata = list()
)

new_uvg_category_representation_from_data(data, category = "category", cues)

new_uvg_category_representation_template_from_data(
  data,
  category = "category",
  cues
)

new_uvg_ideal_observer(
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

new_uvg_ideal_observer_from_data(
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

  Numeric scalar mean.

- sigma2:

  Positive numeric scalar variance (\\\sigma^2 \> 0\\).

- category_labels:

  Character vector of category label(s).

- cue_labels:

  Character vector of cue label(s). For UVG, exactly one cue dimension
  is allowed.

- data:

  Data frame containing category and cue observations.

- category:

  Character string indicating the category column in `data`. Defaults to
  `"category"`.

- cues:

  Character string indicating the cue column in `data`. Defaults to
  `"cue"`.

- category_template:

  An
  [MVBU_CategoryRepresentationTemplate](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md)
  object containing UVG category representations.

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

  Optional perceptual/measurement noise covariance matrix or scalar
  variance.

- noise_treatment:

  Treatment of perceptual noise during evaluation: `"no_noise"`,
  `"sample"`, or `"marginalize"`.

  - `"no_noise"`: No perceptual noise is applied to input cues, and no
    noise is assumed during categorization.

  - `"sample"`: Simulates sensory noise on a specific observation. Noise
    drawn from \\\mathcal{N}(0, \Sigma\_{\text{noise}})\\ is added to
    input cues, and the observer accounts for their own noise
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

An S7 object of the respective UVG class.

## Details

The univariate Gaussian (`UVG`) family models category structure over a
single continuous cue dimension: \$\$p(x \mid c) = \mathcal{N}(x; \mu_c,
\sigma_c^2)\$\$

Category representations are constructed with
new_uvg_category_representation or
new_uvg_category_representation_from_data. Cognitive ideal observer
models are constructed with new_uvg_ideal_observer or
new_uvg_ideal_observer_from_data.

## See also

[family-nix](https://hlplab.github.io/MVBeliefUpdatr/reference/family-nix.md),
[family-mvg](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md),
[family-muvg](https://hlplab.github.io/MVBeliefUpdatr/reference/family-muvg.md),
[family-exemplar](https://hlplab.github.io/MVBeliefUpdatr/reference/family-exemplar.md)

## Examples

``` r
# Create UVG representations
rep_b <- new_uvg_category_representation(
  category_labels = "/b/", cue_labels = "VOT",
  mu = 10, sigma2 = 25
)
rep_d <- new_uvg_category_representation(
  category_labels = "/d/", cue_labels = "VOT",
  mu = 50, sigma2 = 25
)
# Combine into template and ideal observer
tmpl <- new_category_representation_template(list(rep_b, rep_d))
m_uvg <- new_uvg_ideal_observer(category_template = tmpl)
m_uvg
#> <MVBeliefUpdatr::UVG_IdealObserver>
#>   Categories (2): /b/, /d/
#>   Cues (1): VOT
#>   Decision rule: sampling
#>   Category prior: /b/: 0.5, /d/: 0.5
#>   Lapse rate: 0 (treatment: no_lapses)
#>   Perceptual noise: none (treatment: no_noise)
#>   Category representations (2):
#>     /b/: UVG(mu = 10, sigma2 = 25)
#>     /d/: UVG(mu = 50, sigma2 = 25)
```
