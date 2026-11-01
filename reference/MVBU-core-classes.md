# MVBU Core S7 Class Architecture

The MVBeliefUpdatr S7 organizes model objects around an explicit
compositional structure. These classes define the canonical S7 API; any
older list/data-frame-style inputs are supported only through
transitional compatibility methods and should be phased out over time.

Constructs a cognitive model binding a category representation template
to model-level decision and uncertainty parameters.

## Usage

``` r
new_cognitive_model(
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
```

## Arguments

- metadata:

  Optional list of metadata. Defaults to
  [`list()`](https://rdrr.io/r/base/list.html).

- category_template:

  An MVBU_CategoryRepresentationTemplate object containing category
  representations.

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
  Defaults to uniform bias across categories.

- Sigma_noise:

  Optional perceptual/measurement noise covariance matrix or scalar
  variance. Defaults to `NULL`.

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

An MVBU_CognitiveModel object.

## Details

- MVBU_CategoryRepresentation defines the structure for one
  category-level representational object. This could be a uni- or
  multivariate Gaussian, a mixture of Gaussians, or an exemplar-based
  representation.

- MVBU_CategoryRepresentationTemplate collects one or more
  category-representation objects into a validated set used by a model.
  This corresponds to the notion of "templates" in e.g., (Nearey and
  Assmann 2007) .

- MVBU_CognitiveModel combines category-template structure with
  model-level decision behavior (`decision_rule`, `category_prior`,
  `lapse_rate`, `lapse_bias`).

Family-level semantics are layered on top of these base classes. Ideal
observer/adaptor pairings follow:

- [UVG_IdealObserver](https://hlplab.github.io/MVBeliefUpdatr/reference/family-uvg.md)
  `UVG` \<-\>
  [NIX_IdealAdaptor](https://hlplab.github.io/MVBeliefUpdatr/reference/family-nix.md)
  `NIX`

- [MUVG_IdealObserver](https://hlplab.github.io/MVBeliefUpdatr/reference/family-muvg.md)
  `MUVG` \<-\>
  [MNIX_IdealAdaptor](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mnix.md)
  `MNIX`

- [MVG_IdealObserver](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md)
  `MVG` \<-\>
  [NIW_IdealAdaptor](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md)
  `NIW`

while
[Exemplar_Model](https://hlplab.github.io/MVBeliefUpdatr/reference/family-exemplar.md)
`EXEMPLAR` is treated as a standalone family.

## MVBU_Object

Root S7 base class for MVBeliefUpdatr objects.

## MVBU_CategoryRepresentation

Abstract category-level representational class.

Expected properties:

- `category_labels`: one or more category labels associated with the
  object.

- `cue_labels`: one or more cue-dimension labels.

- `category_likelihood_function`: family-specific likelihood function

## MVBU_CategoryRepresentationTemplate

Container class for per-category representation objects.

Expected properties:

- `representations`: list of objects inheriting from
  `MVBU_CategoryRepresentation`.

- `metadata`: optional template-level metadata.

Validators enforce at least one representation and (if named) unique,
non-empty representation names.

## MVBU_CognitiveModel

Abstract model class binding a category-representation template to
model-level decision and uncertainty parameters.

Expected properties:

- `category_template`: `MVBU_CategoryRepresentationTemplate`.

- `category_posterior_functions`: named list of treatment-specific
  posterior functions.

- `decision_rule`: scalar character rule label.

- `category_prior`: numeric probability vector over represented
  categories.

- `lapse_behavior`: list bundling attentional lapse parameters:
  `lapse_rate` (scalar probability in `[0, 1]`), `lapse_bias` (vector of
  lapse response probabilities across categories), and `lapse_treatment`
  (`"no_lapses"`, `"sample"`, or `"marginalize"`).

- `noise_behavior`: list bundling perceptual noise parameters:
  `Sigma_noise` (perceptual noise covariance matrix) and
  `noise_treatment` (`"no_noise"`, `"sample"`, or `"marginalize"`).
  `Sigma_noise` is stored as a square matrix and may be supplied as a
  vector (converted to a diagonal matrix) or as a matrix.

- `metadata`: optional model metadata.

## References

Nearey TM, Assmann PF (2007). “Probabilistic 'sliding template' models
for indirect vowel normalization.” In Solé M, Beddor PS, Ohala M (eds.),
*Experimental Approaches to Phonology*, 246–269. Oxford University
Press, Oxford.

Xie X, Jaeger TF, Kurumada C (2023). “What we do (not) know about the
mechanisms underlying adaptive speech perception: A computational
framework and review.” *Cortex*, **166**, 377–424.

Xie, X., Jaeger, T. F., & Kurumada, C. (2021). What we do (not) know
about the mechanisms underlying adaptive speech perception: A
computational review.
