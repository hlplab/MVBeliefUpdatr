# Backward Compatibility: Working with Legacy Code

## Introduction and Architectural Shift

Earlier versions (pre-v1.0 / pre-v0.1.0) of **MVBeliefUpdatr**
represented multidimensional category distributions, category templates,
and cognitive decision models as nested `tibble` data frames with
list-columns of vectors and parameter matrices (`mu`, `Sigma`, `m`, `S`,
etc.). While convenient for simple tabular wrangling, weakly typed
tabular models suffered from several limitations:

- Lack of parameter dimension and positive-definiteness validation at
  object creation.
- Loosely coupled model metadata, cue labels, decision rules, and noise
  parameters.
- Dispersed and ad-hoc method dispatch across functions and packages.

Version 1.0 rebuilds **MVBeliefUpdatr** around R’s modern **S7** object
system, providing formal class definitions, robust property validators,
unified generics, and a decoupled model family registry (`MVBU_Family`).

To provide a smooth, non-breaking transition path for existing analysis
pipelines and user scripts, **MVBeliefUpdatr** provides an extensive
compatibility bridge implemented in
`R/S7-to-legacy-tibble-compatibility.R` and deprecation mappings.

``` r

library(MVBeliefUpdatr)
library(tibble)
library(dplyr)
```

------------------------------------------------------------------------

## 1. Coercing S7 Objects with `as_tibble()`

Calling
[`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
on any modern S7 category representation, category template, or
cognitive model produces a tibble matching the pre-S7 column schema,
variable names, types, and column ordering:

``` r

# Create an S7 cognitive model
observer_model <- example_mvg_ideal_observer()

# Coerce an S7 cognitive model to a legacy tibble
model_tibble <- as_tibble(observer_model)
model_tibble
#> # A tibble: 6 × 7
#>   category mu           Sigma         prior lapse_rate lapse_bias Sigma_noise
#>   <fct>    <named list> <named list>  <dbl>      <dbl>      <dbl> <list>     
#> 1 /b/      <dbl [1]>    <dbl [1 × 1]> 0.167          0      0.167 <NULL>     
#> 2 /d/      <dbl [1]>    <dbl [1 × 1]> 0.167          0      0.167 <NULL>     
#> 3 /g/      <dbl [1]>    <dbl [1 × 1]> 0.167          0      0.167 <NULL>     
#> 4 /k/      <dbl [1]>    <dbl [1 × 1]> 0.167          0      0.167 <NULL>     
#> 5 /p/      <dbl [1]>    <dbl [1 × 1]> 0.167          0      0.167 <NULL>     
#> 6 /t/      <dbl [1]>    <dbl [1 × 1]> 0.167          0      0.167 <NULL>

# Inspect column names and list-column contents
names(model_tibble)
#> [1] "category"    "mu"          "Sigma"       "prior"       "lapse_rate" 
#> [6] "lapse_bias"  "Sigma_noise"
model_tibble$mu
#> $`/b/`
#>     VOT 
#> 8.39663 
#> 
#> $`/d/`
#>      VOT 
#> 13.81306 
#> 
#> $`/g/`
#>      VOT 
#> 16.71445 
#> 
#> $`/k/`
#>      VOT 
#> 54.04223 
#> 
#> $`/p/`
#>      VOT 
#> 50.23845 
#> 
#> $`/t/`
#>      VOT 
#> 60.48317
model_tibble$Sigma
#> $`/b/`
#>          VOT
#> VOT 23.85665
#> 
#> $`/d/`
#>          VOT
#> VOT 77.63595
#> 
#> $`/g/`
#>          VOT
#> VOT 108.7265
#> 
#> $`/k/`
#>          VOT
#> VOT 418.7947
#> 
#> $`/p/`
#>          VOT
#> VOT 445.1705
#> 
#> $`/t/`
#>          VOT
#> VOT 476.4007
```

When called,
[`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
issues an informative deprecation warning noting that S7 objects are the
canonical format and that behavioral metadata (such as decision noise
parameters, lapses, and internal validators) are not fully represented
in tabular form.

------------------------------------------------------------------------

## 2. Deprecated `dplyr` Verb Forwarding

For legacy workflows that previously piped models directly into `dplyr`
verbs (e.g., `model %>% mutate(...)` or `model %>% filter(...)`), S3
method forwarding intercepts these calls, converts the object to a
tibble, applies the requested verb, and alerts the user:

``` r

# Legacy workflow: calling filter() or select() on an S7 model
filtered_df <- suppressWarnings(filter(observer_model, category == "/i/"))
filtered_df
#> # A tibble: 0 × 7
#> # ℹ 7 variables: category <fct>, mu <named list>, Sigma <named list>,
#> #   prior <dbl>, lapse_rate <dbl>, lapse_bias <dbl>, Sigma_noise <list>

selected_df <- suppressWarnings(select(observer_model, category, mu, Sigma))
selected_df
#> # A tibble: 6 × 3
#>   category mu           Sigma        
#>   <fct>    <named list> <named list> 
#> 1 /b/      <dbl [1]>    <dbl [1 × 1]>
#> 2 /d/      <dbl [1]>    <dbl [1 × 1]>
#> 3 /g/      <dbl [1]>    <dbl [1 × 1]>
#> 4 /k/      <dbl [1]>    <dbl [1 × 1]>
#> 5 /p/      <dbl [1]>    <dbl [1 × 1]>
#> 6 /t/      <dbl [1]>    <dbl [1 × 1]>
```

------------------------------------------------------------------------

## 3. Parameter & Moment Accessors

In pre-v0.1.0 code, parameters were typically extracted by indexing
list-columns:

``` r

# Legacy workflow (pre-v0.1.0):
mu_cat1 <- model_df$mu[[1]]
Sigma_cat1 <- model_df$Sigma[[1]]
```

In the modern S7 architecture, parameters and distributional moments are
standardized via uniform S7 generics, described next.

#### Extracting Model Parameters (`get_parameters`)

The
[`get_parameters()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_parameters.md)
generic retrieves the structured named list of parameters guaranteed by
the `metadata$parameters` property across all four architectural tiers:

``` r

# Extract parameters across a cognitive model
params_model <- get_parameters(observer_model)
params_model
#> $`/b/`
#> $`/b/`$mu
#> [1] 8.39663
#> 
#> $`/b/`$Sigma
#>          VOT
#> VOT 23.85665
#> 
#> 
#> $`/d/`
#> $`/d/`$mu
#> [1] 13.81306
#> 
#> $`/d/`$Sigma
#>          VOT
#> VOT 77.63595
#> 
#> 
#> $`/g/`
#> $`/g/`$mu
#> [1] 16.71445
#> 
#> $`/g/`$Sigma
#>          VOT
#> VOT 108.7265
#> 
#> 
#> $`/k/`
#> $`/k/`$mu
#> [1] 54.04223
#> 
#> $`/k/`$Sigma
#>          VOT
#> VOT 418.7947
#> 
#> 
#> $`/p/`
#> $`/p/`$mu
#> [1] 50.23845
#> 
#> $`/p/`$Sigma
#>          VOT
#> VOT 445.1705
#> 
#> 
#> $`/t/`
#> $`/t/`$mu
#> [1] 60.48317
#> 
#> $`/t/`$Sigma
#>          VOT
#> VOT 476.4007

# Extract parameters from an individual category representation
rep_first <- observer_model@category_template@representations[[1]]
get_parameters(rep_first)
#> $mu
#> [1] 8.39663
#> 
#> $Sigma
#>          VOT
#> VOT 23.85665
```

#### Accessing Expected and Marginal Moments

For Bayesian belief distributions (such as NIW, NIX, and MNIX ideal
adaptors), distributional moments are accessed directly through
dedicated generics rather than manually computed from hyperparameters
($`\kappa, \nu, m, S`$):

``` r

# Instantiate an NIW ideal adaptor model
adaptor_model <- example_niw_ideal_adaptor()

# Expected moments conditional on hyperparameter beliefs:
get_expected_mu(adaptor_model)
#> $`/b/`
#> [1] 8.39663
#> 
#> $`/d/`
#> [1] 13.81306
#> 
#> $`/g/`
#> [1] 16.71445
#> 
#> $`/k/`
#> [1] 54.04223
#> 
#> $`/p/`
#> [1] 50.23845
#> 
#> $`/t/`
#> [1] 60.48317
get_expected_sigma(adaptor_model)
#> $`/b/`
#>          [,1]
#> [1,] 23.85665
#> 
#> $`/d/`
#>          [,1]
#> [1,] 77.63595
#> 
#> $`/g/`
#>          [,1]
#> [1,] 108.7265
#> 
#> $`/k/`
#>          [,1]
#> [1,] 418.7947
#> 
#> $`/p/`
#>          [,1]
#> [1,] 445.1705
#> 
#> $`/t/`
#>          [,1]
#> [1,] 476.4007

# Marginal moments incorporating estimation uncertainty:
get_marginal_mu(adaptor_model)
#> $`/b/`
#> [1] 8.39663
#> 
#> $`/d/`
#> [1] 13.81306
#> 
#> $`/g/`
#> [1] 16.71445
#> 
#> $`/k/`
#> [1] 54.04223
#> 
#> $`/p/`
#> [1] 50.23845
#> 
#> $`/t/`
#> [1] 60.48317
get_marginal_sigma(adaptor_model)
#> $`/b/`
#>         [,1]
#> [1,] 47.7133
#> 
#> $`/d/`
#>          [,1]
#> [1,] 155.2719
#> 
#> $`/g/`
#>         [,1]
#> [1,] 217.453
#> 
#> $`/k/`
#>          [,1]
#> [1,] 837.5894
#> 
#> $`/p/`
#>         [,1]
#> [1,] 890.341
#> 
#> $`/t/`
#>          [,1]
#> [1,] 952.8014
```

------------------------------------------------------------------------

## 4. Transitioning Collections of Models to `MVBU_ModelList`

In legacy workflows, researchers often tracked model adaptation across
trials, blocks, or experimental conditions by storing models inside
tibbles:

``` r

# Legacy workflow:
history_tbl <- tibble(
  step = 1:2,
  model = list(model_step1, model_step2)
)
```

In modern **MVBeliefUpdatr**, collections of models are represented
using the typed `MVBU_ModelList` class:

``` r

# Construct a typed collection of cognitive models
model_seq <- MVBU_ModelList(
  models = list(observer_model, adaptor_model),
  model_labels = c("Initial_Observer", "Adapted_Model")
)
model_seq
#> <MVBU_ModelList>

# Subsetting returns an MVBU_ModelList or individual model:
model_seq[1]
#> <MVBU_ModelList>
model_seq[[1]]
#> <MVBeliefUpdatr::MVG_IdealObserver>
#>   Categories (6): /b/, /d/, /g/, /k/, /p/, /t/
#>   Cues (1): VOT
#>   Decision rule: sampling
#>   Category prior: /b/: 0.167, /d/: 0.167, /g/: 0.167, /k/: 0.167, /p/: 0.167, /t/: 0.167
#>   Lapse rate: 0 (treatment: no_lapses)
#>   Perceptual noise: none (treatment: no_noise)
#>   Category representations (6):
#>     /b/: MVG(mu = 8.397, Sigma = 23.86)
#>     /d/: MVG(mu = 13.81, Sigma = 77.64)
#>     /g/: MVG(mu = 16.71, Sigma = 108.7)
#>     /k/: MVG(mu = 54.04, Sigma = 418.8)
#>     /p/: MVG(mu = 50.24, Sigma = 445.2)
#>     /t/: MVG(mu = 60.48, Sigma = 476.4)
model_seq[["Adapted_Model"]]
#> NULL

# Batch parameter extraction across all models:
get_parameters(model_seq)
#> list()
```

`MVBU_ModelList` objects integrate directly with
[`plot_categories()`](https://hlplab.github.io/MVBeliefUpdatr/reference/plot_categories.md)
and
[`plot_categorization_functions()`](https://hlplab.github.io/MVBeliefUpdatr/reference/plot_categorization_functions.md)
for visualizing incremental belief updating trajectories.

------------------------------------------------------------------------

## 5. API Migration Mapping Table

The table below maps legacy (pre-v0.1.0) functions and naming
conventions to their modern v0.1.0+ S7 API replacements. All legacy
functions issue deprecation warnings in v0.1.0 and are scheduled for
removal in v0.2.0.

| Legacy API Concept / Function (pre-v0.1.0) | Modern S7 API Replacement (v0.1.0+) | Description / Purpose |
|:---|:---|:---|
| [`make_MVG_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/make_MVG_from_data.md) | [`new_mvg_category_representation_template_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md) | Estimate MVG category template from data |
| [`make_MVG_ideal_observer_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/make_MVG_ideal_observer_from_data.md) | [`new_mvg_ideal_observer_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md) | Construct MVG cognitive model from data |
| [`make_NIW_belief_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/make_NIW_belief_from_data.md) | [`new_niw_category_representation_template_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md) | Estimate NIW belief template from data |
| [`make_NIW_ideal_adaptor_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/make_NIW_ideal_adaptor_from_data.md) | [`new_niw_ideal_adaptor_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md) | Construct NIW ideal adaptor model from data |
| `make_NIX_belief_from_data()` | [`new_nix_category_representation_template_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-nix.md) | Estimate NIX belief template from data |
| `make_NIX_ideal_adaptor_from_data()` | [`new_nix_ideal_adaptor_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-nix.md) | Construct NIX ideal adaptor model from data |
| `make_MNIX_belief_from_data()` | [`new_mnix_category_representation_template_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mnix.md) | Estimate MNIX belief template from data |
| `make_MNIX_ideal_adaptor_from_data()` | [`new_mnix_ideal_adaptor_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mnix.md) | Construct MNIX ideal adaptor model from data |
| [`make_exemplars_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/make_exemplars_from_data.md) | [`new_exemplar_category_representation_template_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-exemplar.md) | Estimate Exemplar category template from data |
| [`make_exemplar_model_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/make_exemplar_model_from_data.md) | [`new_exemplar_model_from_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-exemplar.md) | Construct Exemplar model from data |
| [`lift_MVG_to_MVG_ideal_observer()`](https://hlplab.github.io/MVBeliefUpdatr/reference/lift_likelihood_to_model.md) | `lift_representation_to_model()` / [`new_mvg_ideal_observer()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mvg.md) | Wrap MVG category template in cognitive model |
| [`lift_NIW_belief_to_NIW_ideal_adaptor()`](https://hlplab.github.io/MVBeliefUpdatr/reference/lift_likelihood_to_model.md) | `lift_representation_to_model()` / [`new_niw_ideal_adaptor()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md) | Wrap NIW template in ideal adaptor model |
| `lift_NIX_belief_to_NIX_ideal_adaptor()` | `lift_representation_to_model()` / [`new_nix_ideal_adaptor()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-nix.md) | Wrap NIX template in ideal adaptor model |
| `lift_MNIX_belief_to_MNIX_ideal_adaptor()` | `lift_representation_to_model()` / [`new_mnix_ideal_adaptor()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-mnix.md) | Wrap MNIX template in ideal adaptor model |
| [`lift_exemplars_to_exemplar_model()`](https://hlplab.github.io/MVBeliefUpdatr/reference/lift_likelihood_to_model.md) | `lift_representation_to_model()` / [`new_exemplar_model()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-exemplar.md) | Wrap Exemplar template in cognitive model |
| [`is.ideal_adaptor_stanfit()`](https://hlplab.github.io/MVBeliefUpdatr/reference/deprecated-functions.md) | `S7::S7_inherits(x, IdealAdaptorStanfit)` | Type check for S7 Stanfit objects |
| `read_ideal_adaptor_stanfit()` | [`read_stanfit()`](https://hlplab.github.io/MVBeliefUpdatr/reference/read_stanfit.md) | Read saved Stanfit object from RDS |
| `write_ideal_adaptor_stanfit()` | [`write_stanfit()`](https://hlplab.github.io/MVBeliefUpdatr/reference/read_stanfit.md) | Write Stanfit object to RDS |
| `get_info_from_model()` | [`get_parameters()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_parameters.md), [`summary()`](https://rdrr.io/r/base/summary.html), or [`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html) | Extract model parameters / metadata |
| [`get_model_family()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_model_family.md) | [`get_family_name()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_model_family.md) / [`get_family()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_model_family.md) | Retrieve model family name or `MVBU_Family` |
| `get_representation_type()` | [`get_family_name()`](https://hlplab.github.io/MVBeliefUpdatr/reference/get_model_family.md) | Retrieve representation family name |
| [`sample_observation()`](https://hlplab.github.io/MVBeliefUpdatr/reference/sample_observation.md) | [`sample_observations()`](https://hlplab.github.io/MVBeliefUpdatr/reference/sample_observations.md) | Draw random samples from category representations |
| `prepare_staninput()` | [`new_ideal_adaptor_stanfit_input()`](https://hlplab.github.io/MVBeliefUpdatr/reference/new_ideal_adaptor_stanfit_input.md) | Prepare typed input for ideal adaptor Stan models |
| `x_ss` / `c_ss` | `x_css` / `c_css` | Centered sum of squares in sufficient statistics |
| Nested model tibbles | `MVBU_ModelList` | Typed collections of cognitive models |

------------------------------------------------------------------------

## 6. Related Documentation

- [S7 Class Architecture, Constructors, and
  Workflows](https://hlplab.github.io/MVBeliefUpdatr/articles/s7-class-structure-and-workflows.md):
  Detailed technical guide to S7 classes, generics, and constructors.
- [Fitting and Working with MVBeliefUpdatr Stanfit
  Models](https://hlplab.github.io/MVBeliefUpdatr/articles/fitting-and-working-with-stanfit-models.md):
  Fitting, evaluating, and diagnosing Bayesian belief updating models.
- [Visualizing Categories, Categorization Functions, and Belief
  Updating](https://hlplab.github.io/MVBeliefUpdatr/articles/visualizing-models-and-categories.md):
  Plotting 1D, 2D, and 3D category distributions and decision
  boundaries.
- [Extending MVBeliefUpdatr: Adding Model
  Families](https://hlplab.github.io/MVBeliefUpdatr/articles/extending-mvbeliefupdatr-adding-model-families.md):
  Registering custom model families using `MVBU_Family`.
