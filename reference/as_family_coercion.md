# Coerce representations, templates, and cognitive models between families

Coerces between point-estimate families (UVG, MUVG, MVG), uncertainty
families (NIX, MNIX, NIW), and the non-parametric Exemplar family.
Supported coercions are: same-family (returned unchanged), matched pairs
(UVG\<-\>NIX, MUVG\<-\>MNIX, MVG\<-\>NIW), and any family \<-\>
EXEMPLAR. Coercions into NIX/MNIX/NIW require `kappa` and `nu`;
coercions into EXEMPLAR require `n`, the number of exemplars to sample.

## Usage

``` r
as_uvg_category_representation(x)

as_nix_category_representation(x, kappa, nu)

as_muvg_category_representation(x, weights = NULL)

as_mnix_category_representation(x, kappa, nu, weights = NULL)

as_mvg_category_representation(x)

as_niw_category_representation(x, kappa, nu)

as_exemplar_category_representation(x, n)

as_uvg_category_representation_template(x)

as_nix_category_representation_template(x, kappa, nu)

as_muvg_category_representation_template(x, weights = NULL)

as_mnix_category_representation_template(x, kappa, nu, weights = NULL)

as_mvg_category_representation_template(x)

as_niw_category_representation_template(x, kappa, nu)

as_exemplar_category_representation_template(x, n)

as_uvg_ideal_observer(x)

as_nix_ideal_adaptor(x, kappa, nu)

as_muvg_ideal_observer(x, weights = NULL)

as_mnix_ideal_adaptor(x, kappa, nu, weights = NULL)

as_mvg_ideal_observer(x)

as_niw_ideal_adaptor(x, kappa, nu)

as_exemplar_model(x, n)

as_ideal_adaptor(x, model_family = NULL, ...)

as_ideal_observer(x, model_family = NULL, ...)
```

## Arguments

- x:

  An `MVBU_CategoryRepresentation`,
  `MVBU_CategoryRepresentationTemplate`, or `MVBU_CognitiveModel`
  object.

- kappa:

  Strength of belief (pseudocount) about the category mean \\\mu\\.

- nu:

  Degrees of freedom (pseudocount) for covariance/variance beliefs
  \\\Sigma\\.

- weights:

  Optional per-cue weights for MUVG/MNIX targets. (default: preserved
  from `x` if available, otherwise the family constructor's default.)

- n:

  Number of exemplars to sample when coercing into the EXEMPLAR family.

- model_family:

  Optional character string specifying target model family ("NIX",
  "MNIX", "NIW" for adaptors; "UVG", "MUVG", "MVG" for observers).

- ...:

  Additional arguments passed to constructors.

## Value

An object of the same kind as `x` (representation, template, or model),
coerced to the target family.

## Details

Also provides generic lifters `as_ideal_adaptor()` and
`as_ideal_observer()` to lift representations or templates into full
cognitive models.

## Examples

``` r
# Coerce between UVG and NIX representations
uvg_rep <- new_uvg_category_representation(
  category_labels = "A",
  cue_labels = "cue1",
  mu = 0,
  sigma2 = 1
)
nix_rep <- as_nix_category_representation(uvg_rep, kappa = 10, nu = 15)
nix_rep
#> <MVBeliefUpdatr::NIX_CategoryRepresentation>
#>   Category: A
#>   Cues (1): cue1
#>   kappa: 10
#>   nu: 15
#>   m: 0
#>   sigma2: 0.8667

# Coerce back to UVG
as_uvg_category_representation(nix_rep)
#> <MVBeliefUpdatr::UVG_CategoryRepresentation>
#>   Category: A
#>   Cues (1): cue1
#>   mu: 0
#>   sigma2: 1

# Sample exemplars from a parametric representation
ex_rep <- as_exemplar_category_representation(uvg_rep, n = 20)
ex_rep
#> <MVBeliefUpdatr::Exemplar_CategoryRepresentation>
#>   Category: A
#>   Cues (1): cue1
#>   Exemplars (20 points):
#>     [1] (-0.5428)
#>     [2] (0.581)
#>     [3] (0.7682)
#>     [4] (0.4638)
#>     [5] (-0.8858)
#>     ...

# Lift representation to an ideal observer model
obs <- as_ideal_observer(uvg_rep)
obs
#> <MVBeliefUpdatr::UVG_IdealObserver>
#>   Categories (1): A
#>   Cues (1): cue1
#>   Decision rule: sampling
#>   Category prior: A: 1
#>   Lapse rate: 0 (treatment: no_lapses)
#>   Perceptual noise: none (treatment: no_noise)
#>   Category representations (1):
#>     A: UVG(mu = 0, sigma2 = 1)
```
