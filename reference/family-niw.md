# Normal-Inverse-Wishart (NIW) Family: Representations, Templates, and Models

Constructs and manages Normal-Inverse-Wishart (NIW) category
representations, templates, and ideal adaptor cognitive models.

## Usage

``` r
new_niw_category_representation(
  category_labels,
  cue_labels,
  m,
  kappa,
  nu,
  S,
  metadata = list()
)

new_niw_category_representation_from_data(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = length(cues) + 2
)

new_niw_category_representation_template_from_data(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = length(cues) + 2
)

new_niw_ideal_adaptor(
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

new_niw_ideal_adaptor_from_data(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = length(cues) + 2,
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

- m:

  Numeric vector of location hyperparameters.

- kappa:

  Positive numeric scalar precision scaling (\\\kappa \> 0\\).

- nu:

  Numeric scalar degrees of freedom (\\\nu \> D - 1\\).

- S:

  Numeric square symmetric scale matrix.

- category_labels:

  Character vector of category label(s).

- cue_labels:

  Character vector of cue labels.

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
  object containing NIW category representations.

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

  Optional perceptual/measurement noise covariance matrix.

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

## Value

An S7 object of the respective NIW class.

## Details

The Normal-Inverse-Wishart (`NIW`) family induces a multivariate
Student-\\t\\ posterior predictive distribution: For cue dimensionality
\\D\\, \\\nu_t = \nu - D + 1\\ and \\\mathbf{\Sigma}\_t = \frac{\kappa +
1}{\kappa \nu_t} \mathbf{S}\\: \$\$p(\mathbf{x} \mid c) =
t\_{\nu_t}(\mathbf{x}; \mathbf{m}, \mathbf{\Sigma}\_t)\$\$

## See also

[family-mvg](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md),
[family-nix](https://hlplab.github.io/MVBeliefUpdatr/reference/family-nix.md),
[family-mnix](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mnix.md)

## Examples

``` r
# Normal-Inverse-Wishart (NIW) representation for multidimensional cues
m_b <- c(VOT = 10, f0 = 200)
S_b <- matrix(c(25, 5, 5, 36), nrow = 2,
  dimnames = list(names(m_b), names(m_b)))
rep_b <- new_niw_category_representation(
  category_labels = "/b/", cue_labels = c("VOT", "f0"),
  m = m_b, kappa = 5, nu = 10, S = S_b
)
m_d <- c(VOT = 50, f0 = 240)
S_d <- matrix(c(25, -2, -2, 36), nrow = 2,
  dimnames = list(names(m_d), names(m_d)))
rep_d <- new_niw_category_representation(
  category_labels = "/d/", cue_labels = c("VOT", "f0"),
  m = m_d, kappa = 5, nu = 10, S = S_d
)
tmpl <- new_category_representation_template(list(rep_b, rep_d))
m_niw <- new_niw_ideal_adaptor(category_template = tmpl)
m_niw
#> <MVBeliefUpdatr::NIW_IdealAdaptor>
#>   Categories (2): /b/, /d/
#>   Cues (2): VOT, f0
#>   Decision rule: sampling
#>   Category prior: /b/: 0.5, /d/: 0.5
#>   Lapse rate: 0 (treatment: no_lapses)
#>   Perceptual noise: none (treatment: no_noise)
#>   Category representations (2):
#>     /b/: NIW(kappa = 5, nu = 10, m = vector(2), S = matrix(2, 2))
#>       E[mu] = vector(2), E[Sigma] = matrix(2, 2), Sigma_marg = matrix(2, 2)
#>     /d/: NIW(kappa = 5, nu = 10, m = vector(2), S = matrix(2, 2))
#>       E[mu] = vector(2), E[Sigma] = matrix(2, 2), Sigma_marg = matrix(2, 2)
```
