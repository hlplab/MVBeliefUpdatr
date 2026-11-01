# Exemplar Model Family: Representations, Templates, and Models

Constructs and manages non-parametric exemplar-based category
representations, templates, and cognitive models, providing a direct
mathematical bridge between Kernel Density Estimation (KDE) and the
Generalized Context Model (GCM; (Nosofsky 1986; Ashby and Alfonso-Reese
1995) ).

## Usage

``` r
new_exemplar_category_representation(
  category_labels,
  cue_labels,
  exemplars,
  exemplar_weights = NULL,
  c = NULL,
  metadata = list()
)

new_exemplar_category_representation_from_data(
  data,
  category = "category",
  cues,
  c = NULL
)

new_exemplar_category_representation_template_from_data(
  data,
  category = "category",
  cues,
  c = NULL
)

new_exemplar_model(
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

new_exemplar_model_from_data(
  data,
  category = "category",
  cues,
  c = NULL,
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

- exemplars:

  Numeric matrix or data frame of exemplars (rows are exemplars; columns
  correspond to `cue_labels`).

- exemplar_weights:

  Optional numeric vector of exemplar weights in `[0, 1]` summing to 1.
  Defaults to equal weights (`1 / n_exemplars`).

- c:

  Optional positive numeric scalar specifying the sensitivity / inverse
  kernel bandwidth scale. If `NULL` (default), `c` is estimated
  adaptively via Silverman's (1986) multivariate rule of thumb: \$\$h^2
  = \left(\frac{4}{(d + 2) n}\right)^{\frac{2}{d + 4}}, \quad c =
  \frac{1}{2 h^2}\$\$

- category_labels:

  Character vector of category label(s).

- cue_labels:

  Character vector of cue dimension label(s).

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
  object containing Exemplar category representations.

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

An S7 object of class Exemplar_CategoryRepresentation.

## Details

### Relation to the Generalized Context Model (GCM)

In the Generalized Context Model (Nosofsky 1986; Kruschke 1992) ,
psychological similarity between an observed stimulus \\\mathbf{x}\\ and
a stored exemplar \\\mathbf{e}\_i\\ is given by:
\$\$\operatorname{sim}(\mathbf{x}, \mathbf{e}\_i) = \exp\left(-c \cdot
d(\mathbf{x}, \mathbf{e}\_i)^p\right)\$\$ where \\c \> 0\\ is the
perceptual sensitivity parameter and \\d(\mathbf{x}, \mathbf{e}\_i)\\ is
the attention-weighted Minkowski distance: \$\$d(\mathbf{x},
\mathbf{e}\_i) = \left(\sum\_{k=1}^d w_k \|x_k -
e\_{i,k}\|^\tau\right)^{\frac{1}{\tau}}\$\$ with attention weights
\\w_k\\ (\\\sum w_k = 1\\).

Under the Luce choice rule with category response biases \\b_C\\, the
categorization probability in GCM is: \$\$P(C \mid \mathbf{x}) =
\frac{b_C \sum\_{j \in C} \operatorname{sim}(\mathbf{x},
\mathbf{e}\_{C,j})}{ \sum_K b_K \sum\_{j \in K}
\operatorname{sim}(\mathbf{x}, \mathbf{e}\_{K,j})}\$\$

### Density Estimation (KDE) Formulation

In `MVBeliefUpdatr`, category representations are framed in terms of
probability density estimation (Ashby and Alfonso-Reese 1995; Shi et al.
2010) : \$\$p(\mathbf{x} \mid C) = \sum\_{i=1}^n \pi\_{C,i} \\
\mathcal{N}\left(\mathbf{x} \mid \mathbf{e}\_{C,i},
\mathbf{\Sigma}\_{\text{kernel}}\right)\$\$

When the kernel covariance is matched to GCM sensitivity and attention:
\$\$\mathbf{\Sigma}\_{\text{kernel}}^{-1} = 2c \mathbf{W} \iff
\mathbf{\Sigma}\_{\text{kernel}} = \frac{1}{2c} \mathbf{W}^{-1}\$\$ the
Gaussian likelihood is directly proportional to GCM similarity:
\$\$\mathcal{N}\left(\mathbf{x} \mid \mathbf{e}\_i,
\mathbf{\Sigma}\_{\text{kernel}}\right) \propto
\operatorname{sim}(\mathbf{x}, \mathbf{e}\_i)\$\$ When plugged into
Bayes' rule \\P(C \mid \mathbf{x}) = \frac{P(C) p(\mathbf{x} \mid
C)}{\sum_K P(K) p(\mathbf{x} \mid K)}\\, the normalizing constant
cancels out, resulting in an exact mathematical equivalence to the GCM
(Ashby and Alfonso-Reese 1995) .

### Implicit Assumptions on \\\tau\\ and \\p\\

This density-based formulation corresponds to the GCM under two specific
geometric assumptions:

1.  **Euclidean distance metric (\\\tau = 2\\)**: Assumes integral
    (holistic) perceptual dimensions where dimensions combine
    Euclidean-wise (Shepard 1987; Nosofsky 1986) .

2.  **Gaussian similarity gradient (\\p = 2\\)**: Assumes a Gaussian
    similarity decay profile \\\exp(-c \cdot d^2)\\ rather than an
    exponential profile \\\exp(-c \cdot d)\\, characteristic of integral
    dimensions and Gaussian kernel density estimators.

In addition, `MVBeliefUpdatr` generalizes standard GCM by allowing a
full covariance matrix \\\mathbf{\Sigma}\_{\text{kernel}} =
\frac{1}{2c}\operatorname{Cov}(\mathbf{X})\\, accommodating correlated
cue dimensions (e.g., acoustic formant cues).

### Bandwidth Scale \\c\\ and Silverman's Rule

If `c` is not provided (`NULL`), `c` defaults to Silverman's (1986)
multivariate rule-of-thumb bandwidth factor: \$\$h^2 =
\left(\frac{4}{(d + 2) n}\right)^{\frac{2}{d + 4}}, \quad c = \frac{1}{2
h^2}\$\$ which asymptotically minimizes the Mean Integrated Squared
Error (MISE) to the true category distribution.

## References

Nosofsky RM (1986). “Attention, similarity, and the
identification–categorization relationship.” *Journal of Experimental
Psychology: General*, **115**(1), 39.

Shepard RN (1987). “Toward a universal law of generalization for
psychological science.” *Science*, **237**(4820), 1317–1323.

Ashby FG, Alfonso-Reese PC (1995). “Categorization as probability
density estimation.” *Journal of Mathematical Psychology*, **39**(2),
216–233.

Kruschke JK (1992). “ALCOVE: An exemplar-based connectionist model of
category learning.” *Psychological Review*, **99**(1), 22–44.

## See also

[family-mvg](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md),
[family-uvg](https://hlplab.github.io/MVBeliefUpdatr/reference/family-uvg.md)

## Examples

``` r
# Create an exemplar category representation from numeric points
pts <- matrix(c(1, 2, 1.5, 2.5), nrow = 2, byrow = TRUE)
cat_rep <- new_exemplar_category_representation(
  category_labels = "A",
  cue_labels = c("cue1", "cue2"),
  exemplars = pts,
  c = 1.0
)
cat_rep
#> <MVBeliefUpdatr::Exemplar_CategoryRepresentation>
#>   Category: A
#>   Cues (2): cue1, cue2
#>   Exemplars (2 points):
#>     [1] (1, 2)
#>     [2] (1.5, 2.5)

# Create an exemplar model from data
df <- data.frame(
  category = c("A", "A", "B", "B"),
  x = c(1, 2, 5, 6),
  y = c(2, 1, 6, 5)
)
mod <- new_exemplar_model_from_data(
  df,
  category = "category",
  cues = c("x", "y")
)
mod
#> <MVBeliefUpdatr::Exemplar_Model>
#>   Categories (2): A, B
#>   Cues (2): x, y
#>   Decision rule: sampling
#>   Category prior: A: 0.5, B: 0.5
#>   Lapse rate: 0 (treatment: no_lapses)
#>   Perceptual noise: none (treatment: no_noise)
#>   Category representations (2):
#>     A: Exemplar(N = 2, (1, 2), (2, 1))
#>     B: Exemplar(N = 2, (5, 6), (6, 5))
```
