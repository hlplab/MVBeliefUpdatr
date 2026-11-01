# Normal-Inverse-Chi-Squared (NIX) Family: Representations, Templates, and Models

Constructs and manages Normal-Inverse-\\\chi^2\\ (NIX) category
representations, templates, and ideal adaptor cognitive models for 1D
cues under belief updating.

## Usage

``` r
new_nix_category_representation(
  category_labels,
  cue_labels,
  m,
  kappa,
  nu,
  sigma2,
  metadata = list()
)

new_nix_category_representation_from_data(
  data,
  category = "category",
  cues = "cue",
  kappa = nu,
  nu = 3,
  ...
)

new_nix_category_representation_template_from_data(
  data,
  category = "category",
  cues = "cue",
  kappa = nu,
  nu = 3,
  ...
)

new_nix_ideal_adaptor(
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

new_nix_ideal_adaptor_from_data(
  data,
  category = "category",
  cues = "cue",
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
```

## Arguments

- metadata:

  Optional list of metadata. Defaults to
  [`list()`](https://rdrr.io/r/base/list.html).

- m:

  Numeric scalar location hyperparameter.

- kappa:

  Positive numeric scalar precision scaling (\\\kappa \> 0\\).

- nu:

  Positive numeric scalar degrees of freedom (\\\nu \> 0\\).

- sigma2:

  Positive numeric scalar scale parameter (\\\sigma^2 \> 0\\).

- category_labels:

  Character vector of category label(s).

- cue_labels:

  Character vector of cue label(s). For NIX, exactly one cue dimension
  is allowed.

- data:

  Data frame containing category and cue observations.

- category:

  Character string indicating the category column in `data`. Defaults to
  `"category"`.

- cues:

  Character string indicating the cue column in `data`. Defaults to
  `"cue"`.

- ...:

  Additional arguments passed to methods or constructors.

- category_template:

  An
  [MVBU_CategoryRepresentationTemplate](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU-core-classes.md)
  object containing NIX category representations.

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

  Optional perceptual/measurement noise covariance matrix or scalar
  variance.

- noise_treatment:

  Treatment of perceptual noise during evaluation and updating:
  `"no_noise"`, `"sample"`, or `"marginalize"`.

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

## Value

An S7 object of the respective NIX class.

## Details

The Normal-Inverse-\\\chi^2\\ (`NIX`) family models adaptor-style
uncertainty over univariate Gaussian category structure. The predictive
likelihood is a Student-\\t\\ distribution: With \\s_c =
\sqrt{\sigma_c^2(\kappa_c+1)/\kappa_c}\\, \$\$p(x \mid c) =
\frac{1}{s_c} \\ t\_{\nu_c}\\\left(\frac{x-m_c}{s_c}\right)\$\$

## See also

[family-uvg](https://hlplab.github.io/MVBeliefUpdatr/reference/family-uvg.md),
[family-niw](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md),
[family-mnix](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mnix.md)

## Examples

``` r
# Normal-Inverse-Chi-squared (NIX) representation for 1D cues
rep_b <- new_nix_category_representation(
  category_labels = "/b/", cue_labels = "VOT",
  m = 10, kappa = 5, nu = 8, sigma2 = 25
)
rep_d <- new_nix_category_representation(
  category_labels = "/d/", cue_labels = "VOT",
  m = 50, kappa = 5, nu = 8, sigma2 = 25
)
tmpl <- new_category_representation_template(list(rep_b, rep_d))
m_nix <- new_nix_ideal_adaptor(category_template = tmpl)
m_nix
#> <MVBeliefUpdatr::NIX_IdealAdaptor>
#>   Categories (2): /b/, /d/
#>   Cues (1): VOT
#>   Decision rule: sampling
#>   Category prior: /b/: 0.5, /d/: 0.5
#>   Lapse rate: 0 (treatment: no_lapses)
#>   Perceptual noise: none (treatment: no_noise)
#>   Category representations (2):
#>     /b/: NIX(kappa = 5, nu = 8, m = 10, sigma2 = 25)
#>       E[mu] = 10, E[sigma2] = 33.33, sigma2_marg = 40
#>     /d/: NIX(kappa = 5, nu = 8, m = 50, sigma2 = 25)
#>       E[mu] = 50, E[sigma2] = 33.33, sigma2_marg = 40
```
