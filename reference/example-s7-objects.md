# Example S7 objects built from the bundled Chodroff-Wilson data

Helper functions to construct example S7 representations, templates, and
cognitive models from the bundled `mixer6` dataset.

## Usage

``` r
example_uvg_category_representation(n_cues = 1, category = "/b/")

example_nix_category_representation(
  n_cues = 1,
  category = "/b/",
  kappa = 1,
  nu = n_cues + 2
)

example_muvg_category_representation(n_cues = 2, category = "/b/")

example_mnix_category_representation(
  n_cues = 2,
  category = "/b/",
  kappa = 1,
  nu = n_cues + 2
)

example_mvg_category_representation(n_cues = 1, category = "/b/")

example_niw_category_representation(
  n_cues = 1,
  category = "/b/",
  kappa = 1,
  nu = n_cues + 2
)

example_exemplar_category_representation(n_cues = 1, category = "/b/")

example_category_representation(type, n_cues = NULL, category = "/b/", ...)

example_uvg_category_representation_template(n_cues = 1, categories = NULL)

example_nix_category_representation_template(
  n_cues = 1,
  categories = NULL,
  kappa = 1,
  nu = n_cues + 2
)

example_muvg_category_representation_template(n_cues = 2, categories = NULL)

example_mnix_category_representation_template(
  n_cues = 2,
  categories = NULL,
  kappa = 1,
  nu = n_cues + 2
)

example_mvg_category_representation_template(n_cues = 1, categories = NULL)

example_niw_category_representation_template(
  n_cues = 1,
  categories = NULL,
  kappa = 1,
  nu = n_cues + 2
)

example_exemplar_category_representation_template(
  n_cues = 1,
  categories = NULL
)

example_category_representation_template(
  type,
  n_cues = NULL,
  categories = NULL,
  ...
)

example_uvg_ideal_observer(n_cues = 1, categories = NULL, ...)

example_nix_ideal_adaptor(
  n_cues = 1,
  categories = NULL,
  kappa = 1,
  nu = n_cues + 2,
  ...
)

example_muvg_ideal_observer(n_cues = 2, categories = NULL, ...)

example_mnix_ideal_adaptor(
  n_cues = 2,
  categories = NULL,
  kappa = 1,
  nu = n_cues + 2,
  ...
)

example_mvg_ideal_observer(n_cues = 1, categories = NULL, ...)

example_niw_ideal_adaptor(
  n_cues = 1,
  categories = NULL,
  kappa = 1,
  nu = n_cues + 2,
  ...
)

example_exemplar_model(n_cues = 1, categories = NULL, ...)

example_model(type, n_cues = NULL, categories = NULL, ...)
```

## Arguments

- n_cues:

  Integer number of cues: 1 (`"VOT"`), 2 (`"VOT"`, `"f0_semitones"`), or
  3 (`"VOT"`, `"f0_semitones"`, `"vowel_duration"`; maximum is 3).
  Default: 1.

- category:

  Character string giving the category label.

- kappa:

  Strength of belief (pseudocount) about the category mean \\\mu\\.
  Default: 1.

- nu:

  Degrees of freedom (pseudocount) for covariance/variance beliefs
  \\\Sigma\\. Default: `n_cues + 2`.

- type:

  Character string indicating the model types. Valid choices are:

  - `"UVG"`: Univariate Gaussian.

  - `"NIX"`: Normal-Inverse-chi^2.

  - `"MUVG"`: Multivariate Univariate Gaussian.

  - `"MNIX"`: Multivariate Normal-Inverse-chi^2.

  - `"MVG"`: Multivariate Gaussian.

  - `"NIW"`: Normal-Inverse-Wishart.

  - `"EXEMPLAR"`: Exemplar.

- ...:

  Additional arguments passed to underlying constructors.

- categories:

  Character vector of category labels.

## Details

The parameter `n_cues` specifies the number of acoustic cues to include
(up to a maximum of 3 cues):

- `n_cues = 1`: Voice onset time (`"VOT"`). Required for univariate
  models (UVG and NIX).

- `n_cues = 2`: Voice onset time and fundamental frequency (`"VOT"` and
  `"f0_semitones"`).

- `n_cues = 3` (maximum): Voice onset time, fundamental frequency, and
  vowel duration (`"VOT"`, `"f0_semitones"`, and `"vowel_duration"`).

## Examples

``` r
# 1D UVG category representation from bundled data
uvg_rep <- example_uvg_category_representation()
uvg_rep
#> <MVBeliefUpdatr::UVG_CategoryRepresentation>
#>   Category: /b/
#>   Cues (1): VOT
#>   mu: 8.397
#>   sigma2: 23.86

# 2D MVG ideal observer model from bundled data
mvg_mod <- example_mvg_ideal_observer(n_cues = 2)
mvg_mod
#> <MVBeliefUpdatr::MVG_IdealObserver>
#>   Categories (6): /b/, /d/, /g/, /k/, /p/, /t/
#>   Cues (2): VOT, f0_semitones
#>   Decision rule: sampling
#>   Category prior: /b/: 0.167, /d/: 0.167, /g/: 0.167, /k/: 0.167, /p/: 0.167, /t/: 0.167
#>   Lapse rate: 0 (treatment: no_lapses)
#>   Perceptual noise: none (treatment: no_noise)
#>   Category representations (6):
#>     /b/: MVG(mu = vector(2), Sigma = matrix(2, 2))
#>     /d/: MVG(mu = vector(2), Sigma = matrix(2, 2))
#>     /g/: MVG(mu = vector(2), Sigma = matrix(2, 2))
#>     /k/: MVG(mu = vector(2), Sigma = matrix(2, 2))
#>     /p/: MVG(mu = vector(2), Sigma = matrix(2, 2))
#>     /t/: MVG(mu = vector(2), Sigma = matrix(2, 2))
```
