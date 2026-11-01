# Aggregate representations, templates, or cognitive models

Aggregates multiple S7 category representations, category representation
templates, or cognitive models of the same family by computing weighted
averages of their category distributions and model parameters. Accepts
an
[MVBU_ModelList](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_ModelList.md)
container object or a list of cognitive models.

## Usage

``` r
aggregate(x, ...)

aggregate_models(models, weights = NULL, ...)
```

## Arguments

- x:

  An S7 representation, template, cognitive model, or list of such
  objects. When `x` is not an S7 object or list of S7 objects, dispatch
  falls back to [`aggregate`](https://rdrr.io/r/stats/aggregate.html).

- ...:

  Additional S7 objects of the same class (for variadic usage of
  `aggregate`) or arguments passed to methods.

- models:

  An
  [MVBU_ModelList](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_ModelList.md)
  container object or non-empty list of cognitive models.

- weights:

  Optional numeric vector of positive weights with the same length as
  the number of objects being aggregated. If `NULL` (default), uniform
  weights \\1/N\\ are used. Weights are automatically normalized to sum
  to 1.

## Value

An aggregated S7 object of the same class.

An S7 object of the same class as the inputs with aggregated parameters.

## Details

Aggregation operates hierarchically depending on the class of the input:

### Parameter Aggregation Functions

|  |  |  |
|----|----|----|
| **Object Level** | **Parameter / Slot** | **Aggregation Function** |
| *Category Representation (MVG)* | `mu` (mean vector) | Weighted arithmetic mean: \\\sum_i w_i \mu_i\\ |
| *Category Representation (MVG)* | `Sigma` (covariance matrix) | Weighted arithmetic mean: \\\sum_i w_i \Sigma_i\\ |
| *Category Representation (UVG)* | `mu` (mean scalar) | Weighted arithmetic mean: \\\sum_i w_i \mu_i\\ |
| *Category Representation (UVG)* | `sigma2` (variance scalar) | Weighted arithmetic mean: \\\sum_i w_i \sigma^2_i\\ |
| *Category Representation (MUVG)* | `mu` (mean vector) | Weighted arithmetic mean: \\\sum_i w_i \mu_i\\ |
| *Category Representation (MUVG)* | `sigma2` (variance vector) | Weighted arithmetic mean: \\\sum_i w_i \sigma^2_i\\ |
| *Category Representation (NIW)* | `m` (location vector) | Weighted arithmetic mean: \\\sum_i w_i m_i\\ |
| *Category Representation (NIW)* | `S` (scatter matrix) | Weighted arithmetic mean: \\\sum_i w_i S_i\\ |
| *Category Representation (NIW)* | `kappa` (degrees of freedom) | Weighted arithmetic mean: \\\sum_i w_i \kappa_i\\ |
| *Category Representation (NIW)* | `nu` (degrees of freedom) | Weighted arithmetic mean: \\\sum_i w_i \nu_i\\ |
| *Category Representation (NIX)* | `m` (location scalar) | Weighted arithmetic mean: \\\sum_i w_i m_i\\ |
| *Category Representation (NIX)* | `S` (scatter scalar) | Weighted arithmetic mean: \\\sum_i w_i S_i\\ |
| *Category Representation (NIX)* | `kappa` (degrees of freedom) | Weighted arithmetic mean: \\\sum_i w_i \kappa_i\\ |
| *Category Representation (NIX)* | `nu` (degrees of freedom) | Weighted arithmetic mean: \\\sum_i w_i \nu_i\\ |
| *Category Representation (MNIX)* | `m`, `S`, `kappa`, `nu` | Weighted arithmetic means across cues |
| *Category Representation (EXEMPLAR)* | `exemplars` | Row-binding (pooling) of all exemplars: `rbind` |
| *Category Representation (EXEMPLAR)* | `exemplar_weights` | Re-weighted by model weights \\w_i W_i\\ and normalized |
| *Category Representation (EXEMPLAR)* | `c` (sensitivity scalar) | Weighted arithmetic mean: \\\sum_i w_i c_i\\ |
| *Category Template* | `representations` | Aggregates each representation category-by-category |
| *Cognitive Model* | `category_template` | Template aggregation across models |
| *Cognitive Model* | `category_prior` | Weighted arithmetic mean: \\\sum_i w_i P(c)\_i\\ |
| *Cognitive Model* | `lapse_rate` | Weighted arithmetic mean: \\\sum_i w_i \lambda_i\\ |
| *Cognitive Model* | `lapse_bias` | Weighted arithmetic mean: \\\sum_i w_i \beta_i\\ |
| *Cognitive Model* | `Sigma_noise` | Weighted arithmetic mean: \\\sum_i w_i \Sigma\_{\text{noise}, i}\\ |
| *Cognitive Model* | `decision_rule` | Inherited from first model |
| *Cognitive Model* | `noise_treatment` | Inherited from first model |
| *Cognitive Model* | `lapse_treatment` | Inherited from first model |

## See also

[`new_mvg_ideal_observer`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md),
[`new_niw_ideal_adaptor`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md),
[`new_exemplar_model`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-exemplar.md)

## Examples

``` r
m1 <- example_mvg_ideal_observer(n_cues = 2)
m2 <- example_mvg_ideal_observer(n_cues = 2)
aggregate(m1, m2, weights = c(0.5, 0.5))
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
# Create two simple UVG models
rep1 <- new_uvg_category_representation("A", "cue1", mu = 1, sigma2 = 1)
rep2 <- new_uvg_category_representation("A", "cue1", mu = 3, sigma2 = 1)
tpl1 <- new_category_representation_template(list(rep1))
tpl2 <- new_category_representation_template(list(rep2))
mod1 <- new_uvg_ideal_observer(tpl1)
mod2 <- new_uvg_ideal_observer(tpl2)

# Aggregate the models with equal weights
agg_mod <- aggregate_models(list(mod1, mod2), weights = c(0.5, 0.5))
agg_mod
#> <MVBeliefUpdatr::UVG_IdealObserver>
#>   Categories (1): A
#>   Cues (1): cue1
#>   Decision rule: sampling
#>   Category prior: A: 1
#>   Lapse rate: 0 (treatment: no_lapses)
#>   Perceptual noise: none (treatment: no_noise)
#>   Category representations (1):
#>     A: UVG(mu = 2, sigma2 = 1)
```
