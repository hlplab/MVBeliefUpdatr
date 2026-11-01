# Contributor Guide: Adding Permanent Model Families

This guide outlines the systematic procedure for contributing a new
model family as a permanent, first-class citizen in the
**MVBeliefUpdatr** core.

------------------------------------------------------------------------

## 1. Motivation: Family Registration vs. First-Class Core Classes

While **MVBeliefUpdatr** supports lightweight custom model families via
[`new_mvbu_family()`](https://hlplab.github.io/MVBeliefUpdatr/reference/MVBU_Family.md)
and general `MVBU_CategoryRepresentation` instances, promoting
frequently used or foundational families to core built-in classes
provides significant engineering advantages:

1.  **IDE Tab Completion and Discoverability**: Dedicated constructor
    functions (such as
    [`new_uvg_category_representation()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-uvg.md)
    and
    [`new_niw_category_representation()`](https://hlplab.github.io/MVBeliefUpdatr/reference/family-niw.md))
    feature explicit, typed arguments that autocomplete seamlessly in
    IDEs.
2.  **Strict Property Validation at Construction**: Concrete S7 classes
    validate domain constraints immediately upon instantiation (e.g.,
    positive variance, symmetric covariance matrices, probability
    distributions summing to 1).
3.  **High-Performance S7 Method Dispatch**: Direct S7 class hierarchy
    dispatch avoids reflective introspection of metadata or property
    names at runtime.
4.  **Integration with Compiled Stan Generative Models**: Built-in
    classes can interface directly with precompiled Stan models and C++
    routines for Bayesian parameter estimation.

------------------------------------------------------------------------

## 2. Architecture and Design Principles

Every built-in model family adheres to the following structural
conventions:

1.  **Family Descriptors (`MVBU_Family`)**: Model families are
    registered in `.build_builtin_families()` within `R/S7-class.R`.
    Each family encapsulates:

    - `density_function`: Pure mathematical density function supporting
      perceptual noise treatments (`"no_noise"`, `"sample"`,
      `"marginalize"`).
    - `sampling_function`: Vectorized random number generator.
    - `expected_moments`: Extraction of first and second moments (`mu`,
      `Sigma`).
    - `default_prior`: Standard prior configuration constructor.
    - `stanmodel`: Optional list specifying `stanmodel_name` (e.g.,
      `"NIX_ideal_adaptor"`) and `stanfit_class` (e.g.,
      `"NIX_IdealAdaptorStanfit"`), or `NULL` if no Stan model is
      associated.

2.  **Internal Likelihood Closure Parity**: Concrete category
    representation constructors do *not* handwrite custom closure math.
    Instead, they instantiate the S7 class and uniformly derive the
    category likelihood closure directly from the family descriptor:

    ``` r

    rep <- MyFamily_CategoryRepresentation(...)
    rep@category_likelihood_function <- .make_family_likelihood_function(rep)
    ```

3.  **Base-Class-First Logic**: Package methods and algorithms target
    the base class `MVBU_CategoryRepresentation`. Specialized methods
    are only implemented where mathematical optimizations (such as
    Cholesky prefactorization) warrant concrete overrides.

------------------------------------------------------------------------

## 3. Step-by-Step Implementation Procedure

### Step 1: Author S7 Class and Constructors (`R/S7-class-{family}.R`)

Create a dedicated file `R/S7-class-{family}.R` containing: 1. S7 class
definition inheriting from `MVBU_CategoryRepresentation`. 2. S7
validator enforcing strict domain constraints. 3. Constructor
`new_{family}_category_representation()` using
`.make_family_likelihood_function(rep)`. 4. Empirical constructor
`new_{family}_category_representation_from_data()`.

``` r

#' @rdname family-myfamily
#' @export
new_myfamily_category_representation <- function(
  category_labels,
  cue_labels,
  param1,
  param2,
  metadata = list()
) {
  # Rationale for retaining typed class constructor:
  # Concrete typed constructors provide strict property validation, IDE tab
  # completion, and fast S7 method dispatch, while obtaining the likelihood
  # closure uniformly from the registered model family density function.
  rep <- MyFamily_CategoryRepresentation(
    category_likelihood_function = NULL,
    metadata = set_labels(
      metadata,
      cue = as.character(cue_labels),
      category = as.character(category_labels),
      response_category = as.character(category_labels),
      group = if (!is.null(metadata$group)) {
        as.character(metadata$group)
      } else {
        character(0)
      }
    ),
    param1 = as.numeric(param1),
    param2 = as.numeric(param2)
  )
  rep@category_likelihood_function <- .make_family_likelihood_function(rep)
  rep
}
```

### Step 2: Register Family Descriptor (`R/S7-class.R`)

In `R/S7-class.R`, add the family registration inside
`.build_builtin_families()`:

``` r

"MYFAMILY" = new_mvbu_family(
  family_names = list(
    name_short = "MYFAMILY",
    name_full = "My Family Distribution",
    name_category_representation = "MyFamily_CategoryRepresentation",
    name_category_representation_template =
      "MyFamily_CategoryRepresentationTemplate",
    name_cognitive_model = "MyFamily_Model"
  ),
  description = "Description of family.",
  density_function = function(
    x, rep, log = FALSE,
    noise_treatment = "no_noise",
    Sigma_noise = NULL
  ) {
    # Vectorized density evaluation supporting noise treatment
  },
  sampling_function = function(rep, n = 1L) {
    # Vectorized random draw generation
  },
  expected_moments = function(rep) {
    list(mu = ..., Sigma = ...)
  },
  stanmodel = list(
    stanmodel_name = "MYFAMILY_ideal_adaptor",
    stanfit_class = "MYFAMILY_IdealAdaptorStanfit"
  ) # or NULL if analytic
)
```

Ensure `.get_family_from_class_name()` maps
`MyFamily_CategoryRepresentation` to `"MYFAMILY"`.

### Step 3: Implement S7 Methods (`R/S7-methods.R`)

Register specific S7 methods if specialized behavior is required:

- `get_parameters(MyFamily_CategoryRepresentation)`
- `sample_observations(MyFamily_CategoryRepresentation)`
- If dynamic adaptor: `update_template(MyFamily_Model)`

### Step 4: Develop Unit Tests (`tests/testthat/`)

1.  Create `tests/testthat/test-XX-S7-{family}.R` covering:

    - Object creation and validator edge cases.
    - Likelihood evaluation against reference distributions.
    - Noise treatment modes (`"no_noise"`, `"sample"`, `"marginalize"`).
    - Expected moments and sampling fidelity.

2.  Verify cross-family parity by running:

    ``` r

    testthat::test_file("tests/testthat/test-17-S7-model-cross-family-parity.R")
    ```

------------------------------------------------------------------------

## 4. Documentation & Packaging Standards

The package documentation is strictly governed by `roxygen2`
conventions:

### 4.1 Roxygen Rd Grouping & Naming Conventions

1.  **Grouping via `@rdname`**: Group semantically related functions,
    classes, and methods onto a single unified manual topic. For
    instance, all constructors and accessors for `MyFamily` should
    share:

    ``` r

    #' @rdname family-myfamily
    #' @family model families
    ```

2.  **Cross-Linking with `@seealso`**: Explicitly link related model
    families, constructors, and mathematical methods. Example:

    ``` r

    #' @seealso \code{\link{new_mvg_category_representation}},
    #'   \code{\link{get_family}}, \code{\link{evaluate_model}}
    ```

3.  **Internal Helpers**: Internal functions must start with a leading
    dot (`.my_internal_helper`), specify `@keywords internal`, and omit
    export tags (or use `@noRd` if no manual page is generated).

### 4.2 S7 Class Documentation

All S7 class definitions must include comprehensive documentation: -
Mark class topics with `#' @docType class` and `#' @export`. - Document
every property/slot using `#' @property <name> <description>`. - In
constructor functions, document all parameters with their types,
constraints, and explicit defaults (e.g., `(default: 0)`).

### 4.3 File Ordering and `@include` Dependencies

Because `roxygen2` manages `Collate` ordering for S7 class hierarchies,
declare inter-file dependencies explicitly at the top of each R file:

``` r

#' @include S7-class.R
#' @include S7-generics.R
#' @include S7-class-myfamily.R
NULL
```

Always verify that parent classes and generics precede their child
classes and method definitions in the `@include` statements.

------------------------------------------------------------------------

## 5. Coercion Architecture (`.mvbu_coerce_category_representation`)

**MVBeliefUpdatr** provides seamless cross-family conversion via
`.mvbu_coerce_category_representation()` and its template counterpart
`.mvbu_coerce_category_representation_template()` in
`R/S7-family-coercion.R`.

### 5.1 Conversion Pathways

1.  **Identity (Same Family)**: If `from == to`, the object returns
    unaltered.
2.  **Matched Point/Uncertain Pairs**:
    - `UVG` \<-\> `NIX`
    - `MUVG` \<-\> `MNIX`
    - `MVG` \<-\> `NIW` Coercing from point (e.g., `MVG`) to belief
      distribution (e.g., `NIW`) requires explicit confidence parameters
      `kappa` and `nu`:
    - Mean center: `m = mu`
    - Scale matrix: `S = Sigma * (nu - D - 1)` Coercing from belief to
      point converts expected moments:
    - `mu = m`
    - `Sigma = S / (nu - D - 1)`
3.  **Arbitrary Family \<-\> `EXEMPLAR`**:
    - **Parametric to Exemplar**: Generates `n` random points via the
      family’s vectorized `sampling_function()`.
    - **Exemplar to Parametric**: Computes sample moments (`x_mean`,
      `x_css`) from stored exemplar points.

### 5.2 Integrating New Families into Coercion

When introducing a permanent family: 1. Define whether it has a matched
point/uncertain counterpart. 2. Register the family token in the
`allowed` vector in `R/S7-family-coercion.R`. 3. Provide dedicated
coercion wrappers:
`r as_myfamily_category_representation <- function(x, ...) { .mvbu_coerce_category_representation(x, to = "MYFAMILY", ...) }`
4. Verify coercion parity in `tests/testthat/test-07-S7-as-family.R`.

------------------------------------------------------------------------

## 6. Stanfit Ecosystem Integration

Permanent cognitive models intended for Bayesian estimation integrate
with the compiled Stan infrastructure:

### 6.1 Stan Generative Programs (`inst/stan/`)

- Place Stan programs in `inst/stan/{family}_ideal_adaptor.stan`.
- Ensure the family’s `stanmodel$stanmodel_name` property matches the
  Stan program basename without `.stan` extension (e.g.,
  `"NIX_ideal_adaptor"`), and `stanmodel$stanfit_class` matches the
  Stanfit S7 class name (e.g., `"NIX_IdealAdaptorStanfit"`).
- Registering the family dynamically populates bidirectional lookups
  (`.to_stanmodel()` and `.to_family()`) in `.mvbu_family_registry`.
- Use shared licensing/data chunks via `#include` statements in
  `inst/stan/include/`.
- Ensure headers compile automatically into C++ exports in
  `src/stanExports_{family}_ideal_adaptor.h` via `rstantools`.

### 6.2 Data Builders and Input Structures

1.  **Raw Frame Preparation**: Implement preprocessing within
    `.prepare_staninput_frame()` in `R/S7-stanfit-input.R`, ensuring
    cues, categories, and group factors are strictly validated and
    aligned.
2.  **Typed Staninput Class**: Define `MyFamily_IdealAdaptorStaninput`
    in `R/S7-staninput.R` inheriting from `IdealAdaptorStaninput`.
3.  **Data Summarization**: Implement `.summarize_exposure()` and
    `.summarize_test()` blocks mapping observations into numeric arrays,
    matrices, and sufficient statistics.

### 6.3 Parameter Extraction and Translation Tables

Stan parameter draws must translate seamlessly into cognitive model
moments:

1.  Map Stan parameter vectors (e.g., `m`, `S`, `kappa`, `nu`) into
    standardized moments in `.compute_exposure_category_statistic()` in
    `R/S7-stanfit-methods.R`.
2.  Implement model reconstruction hooks in `.reconstruct_model_list()`
    in `R/S7-stanfit-postprocess.R` to reconstruct `MVBU_CognitiveModel`
    instances from posterior draws.

------------------------------------------------------------------------

## 7. Testing & Quality Assurance Standards

### 7.1 Standard Test Fixtures

To guarantee fast, deterministic tests without duplicating data
fixtures:

- **`pb52`**: Standard Peterson & Barney (1952) vowel formant dataset
  (1D and 2D formant cues `F1`, `F2`).
- **`mixer6`**: 6-talker acoustic exposure dataset for multi-group and
  talker adaptation scenarios.
- **Package Generators**: Use standard generators rather than ad-hoc
  setups:
  - [`example_exposure_test_data()`](https://hlplab.github.io/MVBeliefUpdatr/reference/example_exposure_test_data.md)
  - [`example_ideal_adaptor_stanfit_input()`](https://hlplab.github.io/MVBeliefUpdatr/reference/example_ideal_adaptor_stanfit_input.md)
  - [`example_ideal_adaptor_staninput()`](https://hlplab.github.io/MVBeliefUpdatr/reference/example_ideal_adaptor_stanfit_input.md)
  - [`example_mvg_ideal_observer()`](https://hlplab.github.io/MVBeliefUpdatr/reference/example-s7-objects.md)
  - [`example_niw_ideal_adaptor()`](https://hlplab.github.io/MVBeliefUpdatr/reference/example-s7-objects.md)
- **Pre-fitted MCMC Models**: Store serialized test fits strictly in
  `tests/testthat/models/`. Never run long MCMC chains during unit test
  runs.

### 7.2 Parameter Recovery Testing

Every new estimator or Stan model must demonstrate parameter recovery:

1.  Simulate known true parameters using the family’s
    `sampling_function()`.
2.  Fit the model or calculate analytic updates.
3.  Assert that estimated parameters fall within nominal error
    tolerances (e.g., 95% posterior credible intervals or analytical
    equivalence).

### 7.3 Parity Snapshots and Cross-Family Parity

Verify that the new family produces expected likelihoods and
classifications consistent with analytical baselines by adding test
expectations to: -
`tests/testthat/test-17-S7-model-cross-family-parity.R` -
`tests/testthat/test-01-S7-classes.R`

------------------------------------------------------------------------

## 8. Quality Assurance & Phase-End Verification Checklist

Before opening a pull request, ensure every item on this checklist is
satisfied:

All functions, classes, and methods have complete roxygen documentation.

No documentation warnings or unresolved Rd links during
`devtools::document()`.

Code formatting strictly complies with the line length limit (all lines
\<= 80 characters).

All unit tests pass cleanly:

``` r

devtools::test()
```

All package vignettes compile and render without errors:

``` r

rmarkdown::render(
  "vignettes/extending-mvbeliefupdatr-adding-model-families.Rmd"
)
```

Cross-family parity tests pass with 0 failures:

``` r

testthat::test_file("tests/testthat/test-17-S7-model-cross-family-parity.R")
```
