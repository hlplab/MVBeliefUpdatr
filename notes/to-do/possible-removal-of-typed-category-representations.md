# Audit and Roadmap: Possible Removal of Typed Category Representations

## Executive Summary

Up to the completion of Phase 7 of the S7 class migration, the **MVBeliefUpdatr** package defined separate S7 subclasses for each model family across all four hierarchy tiers:
- **Category Representations**: `UVG_CategoryRepresentation`, `MVG_CategoryRepresentation`, `MUVG_CategoryRepresentation`, `NIX_CategoryRepresentation`, `MNIX_CategoryRepresentation`, `NIW_CategoryRepresentation`, `Exemplar_CategoryRepresentation`.
- **Category Templates**: `UVG_CategoryRepresentationTemplate`, etc.
- **Cognitive Models**: `UVG_IdealObserver`, `MVG_IdealObserver`, `MUVG_IdealObserver`, `NIX_IdealAdaptor`, `MNIX_IdealAdaptor`, `NIW_IdealAdaptor`, `Exemplar_Model`.
- **Model Lists**: `UVG_CognitiveModelList`, etc.

With the introduction of the unified **`MVBU_Family`** registry (`register_model_family()`), all family-specific computational logic (density evaluation, observation sampling, expected moments, marginal moments) has been centralized into first-class family objects. 

This document provides:
1. A **complete quantitative audit** of every code location across the package that remains conditional on or dispatches upon typed representation/model subclasses.
2. A **design plan** to ensure existing functions deploy cleanly on base classes (`MVBU_CategoryRepresentation`, `MVBU_CognitiveModel`, etc.) without branching complications.
3. A **step-by-step roadmap** of all further steps required if the package were to become 100% independent of typed subclasses (relying exclusively on generic base classes parameterized by `MVBU_Family`).
4. Documentation of why typed core classes are retained today (IDE tab-completion, strict slot typing).

---

## 1. Quantitative Audit of Typed Subclasses

A full scan of active codebase files (`R/*.R`, excluding deprecated files) identifies **130 `S7_inherits` call sites** and **28 typed S7 method dispatches** across 13 core package files.

### File-by-File Breakdown

| File | Typed References | Primary Role / Responsibility |
| :--- | :---: | :--- |
| `R/S7-methods.R` | 46 | Parameter formatting, summaries, parameter extraction generics (`get_parameters`, `get_expected_*`, `get_marginal_*`). |
| `R/S7-plot-methods.R` | 45 | Cue slicing (`.subset_representation_cues`), family name labelling (`.format_model_family_name`), exemplar point-vs-ellipse plotting. |
| `R/S7-class.R` | 32 | Default family registrations, class definitions, base constructors. |
| `R/S7-update-model.R` | 18 | Analytic Bayesian updating methods (`update_template`, `update_model`) for NIW, NIX, MNIX. |
| `R/basics.R` | 15 | Parameter validation helpers (`.validate_category_representation_parameters`). |
| `R/S7-make-objects.R` | 14 | Model dispatch in `new_cognitive_model()` mapping templates to observer/adaptor subclasses. |
| `R/S7-to-legacy-tibble-compatibility.R` | 12 | Serialization back to legacy S3/tibble data formats. |
| `R/S7-staninput.R` | 9 | Mapping model hyperpriors into Stan data structures. |
| `R/S7-class-*.R` (7 files) | 49 | Typed class declarations, slot definitions, and constructors. |
| `R/S7-stanfit-methods.R` | 6 | S7 method dispatch for Stan fitting routines. |
| `R/S7-plot-engine.R` | 5 | Visual element rendering (distinguishing exemplar scatter from ellipses). |
| `R/S7-stanfit-fitting.R` | 3 | Fitting dispatch for Stan adaptors. |
| `R/S7-stanfit-input.R` | 3 | S7 class inheritance assertions for Stan adaptors. |
| `R/S7-generics.R` | 2 | Generic declarations with typed signatures. |
| **Total** | **259** | |

---

## 2. Categorization by Subsystem

### 2.1 Plotting & Visualization (`R/S7-plot-methods.R`, `R/S7-plot-engine.R`)
- **Total occurrences**: 50
- **Current usage**:
  1. `.format_model_family_name(x)` (14 sites): Translates representation/model types to strings like `"1D univariate Gaussian category"`.
  2. `.subset_representation_cues(r, cues)` (20 sites): Slices a multidimensional representation down to 1D or 2D for marginal plotting. Different representations slice differently (e.g. slicing mean vector and covariance submatrix for MVG; slicing `nu`, `kappa`, `m`, and `S` for NIW; slicing columns of matrix for Exemplar).
  3. Exemplar scatter vs. parametric ellipse (16 sites): Checking `S7_inherits(r, Exemplar_CategoryRepresentation)` to decide whether to plot discrete points (`geom_point`) or density contour ellipses (`stat_ellipse` / polygon paths).
- **Remediation / General Route**:
  - `.format_model_family_name(x)` can query `get_family(x)@name_long` and dimension directly on base classes.
  - Slicing cues can be handled via a family slot `slice_cues(rep, cue_indices)` on `MVBU_Family`.
  - Exemplar plotting check can query `identical(get_family(r)@name_short, "EXEMPLAR")` or a family property `is_parametric`.

### 2.2 Formatting, Summary, & Parameter Extraction (`R/S7-methods.R`)
- **Total occurrences**: 46
- **Current usage**:
  1. `get_parameters` (7 typed method dispatches): Returns named list of parameters for each typed class.
  2. `format` for representations (7 sites): Prints formatted representation string like `UVG(mu = 0, sigma2 = 1)`.
  3. `format` for models (7 sites): Prints model name and summary.
  4. `.format_rep_params` (7 sites): Produces inline parameter strings.
  5. `.rep_to_param_df` (7 sites): Converts parameters to row of tibble for `summary(template)`.
  6. Adaptor expectation extractors (11 sites): `get_expected_mu`, `get_expected_sigma`, `get_marginal_mu`, `get_marginal_sigma` for NIX, MNIX, NIW.
- **Remediation / General Route**:
  - Implement a generic `get_parameters(MVBU_CategoryRepresentation)` that inspects non-metadata slots (`names(S7::props(x))`) or a `parameters` property.
  - Implement a generic `format(MVBU_CategoryRepresentation)` that formats `sprintf("%s(%s)", get_family_name(x), .format_rep_params(x))`.
  - Delegate expectation extractors directly to `get_family(x)@expected_moments(x)` and `get_family(x)@marginal_moments(x)`.

### 2.3 Object Construction & Model Dispatch (`R/S7-make-objects.R`)
- **Total occurrences**: 14
- **Current usage**:
  - `new_cognitive_model(template, ...)` branches on template representation class (`UVG`, `MVG`, `NIW`, etc.) to construct the corresponding typed model (`UVG_IdealObserver`, `NIW_IdealAdaptor`, etc.).
- **Remediation / General Route**:
  - If representations are untyped, `new_cognitive_model` queries `get_family(template)`. The family object determines whether it is an observer or adaptor (e.g. via `fam@is_adaptor` or registry), and instantiates a generic `MVBU_CognitiveModel`.

### 2.4 Analytic Belief Updating (`R/S7-update-model.R`)
- **Total occurrences**: 18
- **Current usage**:
  - `update_template` and `update_model` dispatch on `NIW_IdealAdaptor`, `NIX_IdealAdaptor`, `MNIX_IdealAdaptor`.
- **Remediation / General Route**:
  - Add an optional `update_function` slot to `MVBU_Family`:
    `update_function = function(rep, observations, ...)`
  - Base method `update_model(MVBU_CognitiveModel, observations, ...)` simply delegates:
    `new_rep <- fam@update_function(rep, obs)`

### 2.5 Input Validation (`R/basics.R`)
- **Total occurrences**: 15
- **Current usage**:
  - `.validate_category_representation_parameters` checks parameter validity (positive variance, positive definite covariance, positive degrees of freedom) conditional on family name.
- **Remediation / General Route**:
  - Add an optional `validate_parameters` closure to `MVBU_Family`.

### 2.6 Stan Fitting & Serialization (`R/S7-stan*.R`, `R/S7-to-legacy-tibble-compatibility.R`)
- **Total occurrences**: 39
- **Current usage**:
  - Mapping model representations into Stan data structures (`NIX_IdealAdaptorStaninput`, etc.) and converting S7 objects back into legacy tibble formats.
- **Remediation / General Route**:
  - Stan routines can query family name or a Stan template mapping closure registered with the family.

---

## 3. Deployment Strategy: Base Classes First

To prevent existing functions from becoming overly complicated by dual-handling branches:
1. **Always deploy primary logic on base classes**:
   - Write methods for `MVBU_CategoryRepresentation`, `MVBU_CategoryRepresentationTemplate`, `MVBU_CognitiveModel`, and `MVBU_CognitiveModelList`.
   - Base methods rely exclusively on polymorphic getters (`get_family()`, `get_labels()`, `get_parameters()`, `get_category_likelihood_function()`, `expected_moments()`).
2. **Typed subclasses delegate to base methods**:
   - Typed subclasses should not duplicate logic. If a typed subclass requires specialization, it calls the base method and augments the result.
3. **No dual-branching spaghetti**:
   - Avoid `if (inherits(x, "UVG_CategoryRepresentation") || inherits(x, "MVBU_CategoryRepresentation"))`.
   - Since all typed subclasses inherit from `MVBU_CategoryRepresentation`, an S7 method dispatched on `MVBU_CategoryRepresentation` automatically covers all typed subclasses as well as custom untyped instances.

---

## 4. Complete Roadmap for Eventual Removal of Typed Representations

If future maintenance desires the complete removal of typed subclasses in favor of pure `MVBU_CategoryRepresentation(family = ...)` objects, the following steps must be completed:

### Phase 1: Property Storage Standardization
- Replace typed slots (`rep@mu`, `rep@sigma2`, `rep@S`, etc.) with a uniform `parameters = S7::class_list` property on `MVBU_CategoryRepresentation`, or allow dynamic properties via a shared dictionary.
- Update `get_parameters(rep)` to return `rep@parameters`.

### Phase 2: Family-Level Cue Slicing
- Add `slice_cues = S7::class_function` to `MVBU_Family`.
- Migrate `.subset_representation_cues()` logic from `R/S7-plot-methods.R` into family definitions in `R/S7-class.R`.

### Phase 3: Family-Level Analytic Updating
- Add `update_function = S7::class_function` to `MVBU_Family` for adaptor families.
- Replace typed `update_model` / `update_template` methods in `R/S7-update-model.R` with a single base dispatch on `MVBU_CognitiveModel`.

### Phase 4: Family-Level Validation Closures
- Move parameter validation rules from `R/basics.R` into a `validate_parameters` closure in `MVBU_Family`.

### Phase 5: Stan Integration Generalization
- Register Stan data preparation closures within the family objects (`fam@stan_input_function`).

### Phase 6: Subclass Deprecation and Removal
- Mark `UVG_CategoryRepresentation`, `MVG_IdealObserver`, etc., as deprecated wrappers around `new_category_representation()` and `new_cognitive_model()`.
- After one release cycle, remove the 7 `R/S7-class-*.R` files and associated typed tests.

---

## 5. Rationale for Retaining Core Typed Subclasses Today

Despite the architectural elegance of pure family objects, the core typed classes are intentionally preserved in **MVBeliefUpdatr** for two critical developer and user experience reasons:

1. **Interactive IDE Tab-Completion**:
   - In modern R development environments (RStudio, VSCode, Positron), typed S7 classes provide instant autocompletion for object properties (e.g. `rep@mu`, `rep@sigma2`, `rep@kappa`).
   - Dynamic parameter lists in generic objects do not expose named slots to static analysis tools or autocomplete engines.

2. **Strict S7 Slot Typing and Early Validation**:
   - Dedicated S7 class definitions enforce strict slot constraints (e.g. `mu = S7::class_numeric`, `Sigma = S7::class_matrix`) at instantiation time.
   - This catches dimensional or type mismatches immediately at construction, long before complex numerical computations, plotting pipelines, or Stan MCMC sampling routines are invoked.
