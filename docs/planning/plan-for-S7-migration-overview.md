## MVBeliefUpdatr v1.0 Overhaul Plan (Editable)

### Scope Summary
- Big-bang v1.0 redesign
- Aggressive renaming in new code
- New migration aliases isolated in deprecated-2026.R for 2-3 minor releases (legacy deprecated-2025.R remains historical/out-of-scope)
- Explicit S7 composition: model contains representation
- Stan expansion/debugging deferred until architecture stabilizes
- Extensible family architecture: support multiple ideal observer/adaptor families (including MUVG, MNIX, and other non-Gaussian representations)

### Phase Status
| Phase | Status | Gate | Notes |
| --- | --- | --- | --- |
| 0 Architecture Contract and Freeze | Completed | Closed | Contract approved; interop bridge expansion continues in Phases 1-3 via watchlist (not a blocker). |
| 1 S7 Hierarchy Foundation | Completed | Closed | Base classes, validators, core generic aliases, extension hooks, and baseline dispatch tests are in place. |
| 2 Concrete Class Migration | Completed | Closed | Legacy->S7 adapters implemented for NIW/MVG/Exemplar and prototype MUVG/MNIX families; constructor/default normalization and S7-only adapter-output gates are covered in migration tests. |
| 3 API Unification and Method Coverage | Completed | Closed | S7 unification landed for categorization/prediction, accessors (`get_model_type`, `get_representation_type`, `get_noise_treatment`, `get_lapse_treatment`), Stanfit/Staninput classes, model aggregation (`aggregate_models()`), and unified plotting engine (`plot_categories`, `plot_categorization_function`). Vignettes and tests verified. |
| 4 Compatibility Shell | Completed | Closed | Package version updated to 0.1.0 (deprecation retirement target 0.2.0); standardized roxygen documentation, lifecycle warnings, and seealso links across all 22 deprecated files; as_tibble/dplyr compatibility verified. |
| 5 Data Model and Print Strategy | Completed | Closed | S7 print and summary methods implemented; as_tibble legacy bridge implemented with pre-S7 column fidelity. |
| 6 Performance and Caching Framework | Completed | Closed | Cached posterior closures implemented in plot engine; performance benchmarks in plotting vignette. |
| 7 Consistency and Quality Hardening | Completed | Closed | Dynamic/unified plotting framework, data extraction consistency, dataset ingestion, repeated chunk consolidation, test helper simplification, worked decision-rule examples, parameter naming harmonization, metadata unification, and roxygen tag standardization completed. All 2,319 unit tests passing with zero failures. |
| 7B Test Suite Architecture Cleanup | Completed | Closed | Reorganized into sequential numbering: active tests (`test-00-` to `test-39-`), deprecated tests (`test-80-` to `test-90-`). All 2,905+ tests passing. |
| 8 Documentation and Vignettes | Completed | Closed | All five vignettes authored and verified; Developer Guide for permanent model families authored (`docs/developer/contributor-guide-adding-permanent-model-families.md`). |
| 9 Stan Expansion Readiness & Prior Refactoring | Completed | Closed | 14 Stan models unified with customizable priors (`prior()`), partial parameter fixing (`fix_parameter()`), Cauchy/Beta/Dirichlet hyperpriors, and 4-case conjugate updating for known moments (μ₀, Σ₀). Full test suite (2,905 tests) passing. |

### Locked Decisions
- [x] Release strategy: big-bang v1.0
- [x] Naming direction: aggressive + aliases in deprecated-2026.R
- [x] Accessor/introspection syntax: prefer get_* for consistency
- [x] Backward wrappers: keep for 2-3 minor releases
- [x] Core composition: explicit model -> representation
- [x] Stan timeline: after class/API stabilization
- [x] Finalize cache profile defaults (standard profile default, with
  per-call override)
- [x] Documentation policy: all new code and all code integrated into the new S7
  scaffold must be roxygen documented
- [x] Documentation QA policy: each phase-end cleanup includes checks for broken
  links and roxygen/Rd issues

### Model Relationship Map (Locked Semantics)
- Family pairings encode observer/adaptor uncertainty relationships:
  - UVG (Univariate Gaussian, observer) <-> NIX (Normal-Inverse-chi^2, adaptor over UVG parameters)
  - MUVG (Multi-cue Univariate Gaussian integration, observer) <-> MNIX (Mixture of Normal-Inverse-chi^2, adaptor)
  - MVG (Multivariate Gaussian, observer) <-> NIW (Normal-Inverse-Wishart, adaptor over MVG parameters)
- Exemplar is intentionally treated as a standalone family (not an observer/adaptor conjugate pair in current scaffold).
- Core class composition semantics:
  - CategoryRepresentation: one category-level parametric/nonparametric representation object
  - CategoryRepresentationTemplate: a validated set of per-category representations
  - CognitiveModel: decision/lapse/prior machinery composed with a CategoryRepresentationTemplate

### Decision Blocks

#### Decision: Cache Semantics
- Global default: pure-by-default methods
- Persistent cache candidates (inside Stanfit objects):
  - extracted posterior draws
  - reused posterior summaries
  - bounded-size item-level posterior predictions
- Non-persistent candidates (outside core model objects):
  - large plotting grids
  - dense ad hoc evaluation surfaces
- Implementation idea: cache profiles via package options
- Final choice:
  - Default profile: [ ] minimal  [x] standard  [ ] eager
  - Allow per-call override: [x] yes  [ ] no

## Phase Plan

### Phase 0: Architecture Contract and Freeze
**Goal:** lock contract before coding

Checklist:
- [x] Finalize class family contract: Representation, CognitiveModel
- [x] Confirm that all user-facing core entities are model objects that compose representations
- [x] Finalize family-extension naming pattern for future model families (e.g., <Family>_Representation, <Family>_IdealObserverModel, <Family>_IdealAdaptorModel, <Family>_IdealAdaptorFit)
- [x] Lock naming conventions and migration vocabulary
- [x] Lock wrapper/deprecation policy text and timeline
- [x] Define phase gate criteria and rollback criteria

Phase gate (exit criteria):
- [x] Written architecture contract approved
- [x] Written naming/deprecation policy approved
- [x] Phase gate checklist approved

### Phase 1: S7 Hierarchy Foundation
**Goal:** establish base class and generic system

Checklist:
- [x] Implement abstract S7 base classes (Representation, CognitiveModel)
- [x] Add schema and semantic validators
- [x] Add composition structure (CognitiveModel includes representation slot)
- [x] Create core S7 generics: construct, validate, summarize, print, categorize/predict, posterior, plot-prep
- [x] Define extension hooks for future Stan families
- [x] Add explicit extension hooks for non-Gaussian representation families
- [x] Initialize interop bridge watchlist with baseline forwarded methods and dependency rationale

Phase gate:
- [x] Class instantiation and validator tests pass
- [x] Core generic dispatch tests pass

### Phase 2: Concrete Class Migration
**Goal:** move NIW/MVG/exemplar families to S7 and unify inferred classes

Checklist:
- [x] Migrate representation classes (NIW belief, MVG representation, exemplar representation)
- [x] Migrate cognitive model classes (NIW adaptor, MVG observer, exemplar model)
- [x] Standardize constructors and defaults
- [x] Validate migration pattern can be reused by at least one future-family prototype (MUVG, MNIX, or another non-Gaussian family)

Phase gate:
- [x] Constructor normalization tests pass
- [x] No mixed S4/S7 construction paths for migrated classes

### Phase 3: API Unification and Method Coverage
**Goal:** normalize user-facing behavior and signatures

Checklist:
- [x] Replace mixed S3/S4/direct dispatch with unified S7 API surface (core model, representation, plotting, and stanfit methods unified; legacy periphery functions in `update-decision-model.R` and wrappers in `deprecated-*` scheduled for Phase 5 cleanup).
- [x] Update handling of stanfit-related classes (completed in `R/S7-stanfit.R`, `R/S7-stanfit-input.R`, `R/S7-stanfit-methods.R`, and `R/S7-stanfit-utils.R`).
- [x] Normalize categorization/prediction signatures across model families
- [x] Support both single and list/batch forms consistently
- [x] Unify high-level plotting entry points with internal dimension specialization (`plot_categories()`, `plot_categorization_function()`, `plot_parameters()`, `plot_cue_correlations()`, `plot_cue_densities()`)
- [x] Add family-agnostic accessors: `get_model_type()`, `get_representation_type()`, `get_noise_treatment()`, `get_lapse_treatment()`
- [x] Clean up obsolete aliases and arguments (`sample_observation` -> `sample_observations()`, retired `wide` from `get_draws()`)
- [x] Fill method coverage gaps: S7 likelihood, categorization, evaluation, model aggregation (`aggregate_models()`), and legacy info consumers unified.
- [x] Ensure generic signatures remain family-agnostic for Gaussian and non-Gaussian model families
- [x] Close interop bridge watchlist items based on workflow tests and dependency budget (see `docs/planning/interop-bridge-watchlist.md`)

Phase gate:
- [x] Cross-family API parity tests pass
- [x] Signature consistency tests pass

### Phase 4: Compatibility Shell
**Goal:** keep old API available but isolated

Checklist:
- [x] Consistency Check: Standardize all deprecated functions in `deprecated-*.R` to use `lifecycle::deprecate_warn("0.1.0", ...)` (retiring in `0.2.0`), uniform Title format `Deprecated: <func>`, Description format with lifecycle badge, `@keywords internal`, `@export`, and `@seealso` pointing to exported public functions.
- [x] Ensure wrappers are thin adapters only (all deprecated functions issue lifecycle warnings and delegate directly to S7 constructors, generics, or helper methods).
- [x] Add migration mapping table old -> new API to the vignette for backward compatibility (documented in `vignettes/backward-compatibility-working-with-old-code.Rmd` Section 3).
- [x] Implement `as_tibble()` and deprecated `dplyr` methods (`mutate`, `filter`, etc.) for S7 objects to provide smooth backward compatibility for legacy tibble workflows (implemented in `R/S7-to-legacy-tibble-compatibility.R` and tested in `test-27`)

Phase gate:
- [x] Wrapper output-equivalence tests pass
- [x] Deprecation warnings fire as expected

### Phase 5: Data Model and Print Strategy
**Goal:** remove tibble-as-object while preserving readable printing

Checklist:
- [x] Store canonical fields once in S7 slots (no constant duplication across rows)
- [x] Keep human-readable print and summary methods for S7 models and templates
- [x] Add conversion helpers (`to_tibble()`, `as_tibble()`, `as_draws_df()`, `as_matrix` where needed)

Phase gate:
- [x] Object-size regression checks pass
- [x] Print snapshot tests pass

### Phase 6: Performance and Caching Framework
**Goal:** fix high-impact bottlenecks with controlled memory use

Planning note: see [docs/planning/phase6-performance-caching-framework.md](docs/planning/phase6-performance-caching-framework.md) for the detailed design and implementation outline for the posterior-kernel caching strategy.

Checklist:
- [x] Implement operation-aware caching for decision rules and likelihood/posterior closures in cognitive models (`get_category_likelihood_function`, `get_category_posterior_function`)
- [x] Functional S7 copy-on-modify cache semantics (no state leakage; models start with clean `@cache`)
- [x] Optional post-processing helper `add_parameter_draws()` on `MVBU_Stanfit` to pre-extract and cache draw matrices, eliminating duplicate `get_draws()` extraction work
- [x] Expanded model criteria helper `add_criterion()` supporting `loo`, `waic`, `kfold`, `loo_subsample`, `bayes_R2`, `loo_R2`, and `marglik` (via `bridgesampling`), following `brms` argument conventions
- [x] Model comparison helper `loo_compare()` / `compare_models()` on `MVBU_Stanfit` objects
- [x] Posterior predictive checks via `pp_check.MVBU_Stanfit()` leveraging `bayesplot`
- [x] Deduplication of likelihood evaluation paths by dispatching `likelihood()` methods through `get_category_likelihood_function()`

Priority hotspots to benchmark:
- [x] density plotting recomputation loops (optimized with cached closures)
- [x] repeated draw extraction paths (optimized with `add_parameter_draws`)
- [x] categorization function generation pipelines (optimized with cached vectorized closures)

Phase gate:
- [x] Cache test suite (`test-29-S7-caching.R`) passes cleanly (100% pass)
- [x] Full parity with `brms` model criteria, comparison, and posterior predictive checking verified

### Phase 7: Consistency and Quality Hardening
**Goal:** lock reliability, unify interfaces, and eliminate redundancy across the codebase

Checklist:
- [x] Dynamic & unified plotting framework:
  - Sequential belief updates across exposure blocks (`plot_model_updates()`)
  - Prior-to-posterior parameter transitions across MCMC posterior draws from `MVBU_Stanfit`
  - Unified `plot_sample()`, `plot_exposure_sample()`, and `plot_test_sample()` supporting 1D, 2D, and 3D observations with analytical slicing/marginalization
- [x] Data extraction & schema consistency:
  - Unified `get_data()`, `get_exposure_data()`, and `get_test_data()` with filtering, subsampling, and column name mapping
  - Standardized internal `group_unique` column tracking across inputs and fits
  - Consolidated roxygen documentation under `man/get_data.Rd`
- [x] Plotting deduplication:
  - Audit and consolidate overlapping grid evaluation, contour, and mesh rendering pipelines across `plot_categories()`, `plot_categorization_functions()`, and `plot_sample()`
  - Renamed generic and methods from `plot_categorization_function` to `plot_categorization_functions` following plural convention without keeping singular wrapper
  - Deduplicated cue limit computation (`.compute_default_cue_limits()`) and Plotly surface traces (`.add_2d_surface_trace()`) in `R/S7-plot-engine.R`
- [x] Phonetic datasets integration and standardization:
  - Ingested and compressed datasets in `data/`: `h95.rda`, `pb52.rda`, `swehvd.rda`, `mixer6.rda` via reproducible script `data-raw/import_phonetic_datasets.R`
  - Standardized schema across all datasets: consistent column ordering (`speaker` → demographics → category → context/sub-features → trial/task → acoustic cues → duration/rate → quality flags), `sex` column (`"female"`, `"male"`), standardized IPA notation in phonemic slashes (`/.../`), and `stop_poa`
  - Full roxygen documentation in `R/data.R`, BibTeX entries in `inst/REFERENCES.bib` (`hillenbrand1995`, `peterson-barney1952`, `barreda2015`, `persson2021`)
  - Deprecated `ChodroffWilson2018` in favor of `mixer6`
  - Added dedicated test suite `tests/testthat/test-00-datasets.R` (45 assertions) and updated `.example_data()` and `test-15-S7-from-data.R`
- [x] Codebase repeated chunk consolidation:
  - Undeprecated conjugate conversions (`get_expected_Sigma_from_S`, `get_S_from_expected_Sigma`, `get_expected_mu_from_m`, `get_m_from_expected_mu`) and moved to `R/S7-niw-nix-conversions.R`
  - Moved legacy posterior predictive and `get_D` to `R/deprecated-NIX-NIW-basics.R` with `@keywords internal`
  - Consolidated sufficient statistics computation into `get_sufficient_category_statistics()` and integrated into batch `update_template()` across NIW, NIX, and MNIX models
  - Deprecated single-table SS functions (`get_sum_of_squares_from_df`, aliases) and `make_vector_column` with `@keywords internal`
  - Added `.summarise_category_parameter_draws()` in `R/internal-utils-imported.R` and centralized default cue limits in `R/S7-plot-engine.R`
- [x] Test helper simplification:
  - Eliminated legacy `helper-vowel-data.R` / `make_vowel_test_data` in favor of package datasets (`pb52`, `h95`, `swehvd`, `mixer6`) and standard fixtures
- [x] Vignette enhancement:
  - Added explicit worked example demonstrating the output and return structure of the 3 decision rules (`"criterion"`, `"proportional"`, `"sampling"`) and the `simplify` argument in `vignettes/s7-class-structure-and-workflows.Rmd`
- [x] make extensions of plot_categories and plot_categorization_functions functions that animate model updates (both for histories of update_model or for stanfit ideal adaptors by stepping from prior to posterior for a number of equally-spaced draws) — *Completed via `plot_model_updates()`, `reconstruct_update_history()`, and `plot_categories(..., mode = "animate")`.*
+ [x] In `vignettes/s7-class-structure-and-workflows.Rmd`, reformatted `categorize_decision_rules` chunk and long comments to comply with ≤ 80 characters per line.
+ [x] In `add_posterior_latents`, replaced direct Gaussian density evaluations with model-type specific predictive densities (`.evaluate_stan_category_density`) dispatching on NIW (multivariate Student-t via `.dmvt_density`), NIX (univariate Student-t via `stats::dt`), and MNIX (product of independent univariate Student-t densities).
+ [x] Refactored `categorize()`:
  - Returns `response_category` (first column) and `response_probability` (second column).
  - `response_probability` is 1.0 for deterministic (`"criterion"`) or sampled (`"sample"`) decisions, and equal to the posterior probability under `"proportional"`.
  - Added argument `simplify = if (identical(decision_rule, "proportional")) FALSE else TRUE` returning a character vector of chosen categories when `TRUE`.
  - Updated downstream tests, methods, and vignettes.
+ [x] In `summary` and `print` outputs of `MVBU_Stanfit` and `MVBU_StanfitPosterior`, sorted `Sigma_marg` under `Sigma_exp` for each group and category combination using `factor(Parameter, levels = c("mu", "Sigma_exp", "Sigma_marg"))`.
+ [x] Cleaned up cache structure in `MVBU_Stanfit`:
  - `model@cache` now strictly holds `$draws` and `$summary`.
  - Removed redundant top-level `$category_moments`, `$expected_moments`, and `$marginal_moments`.
  - Updated getters (`get_expected_category_statistic`, `get_marginal_category_statistic`, etc.) to read/write cleanly to `x@cache$summary`.
+ [x] Updated documentation for `add_posterior_latents` and `MVBU_Stanfit`:
  - Documented all quantities stored in `cache$draws` and `cache$summary`.
  - Clarified that point estimates and quantiles in `summary` are obtained across MCMC posterior draws.
  - Linked to `vignette("fitting-and-working-with-stanfit-models")`.
+ [x] Fixed 3D Plotly animated output in `plot_categories()` (`.build_3d_model_list_mesh_plot`):
  - Removed dummy surface vertex marker (`x = c(x_all[1L]), ...`). Added category center markers at μ with `frame = factor(...)` matching the animation.
  - Fixed mesh colorscale (`colorscale = list(c(0, col_k), c(1, col_k))` with `intensity = 0`) to match category colors.
  - Enforced chronological slider and frame ordering (`Prior`, `Exposure_Step1`, ...) instead of Plotly's default alphabetical sorting.
+ [x] Standardized `reconstruct_update_history()` parameter naming:
  - Replaced legacy `full_posterior = TRUE/FALSE` with `uncertainty_treatment = "marginalize"` / `"discard"`.
  - Standardized treatments across all functions (`noise_treatment`, `lapse_treatment`, `uncertainty_treatment`).
  - Updated dynamic quantiles in `plot_categories()` / `plot_categorization_functions()` based on `levels`.
- [x] Singular / plural naming audit:
  - Systematically reviewed all exported functions and argument names for consistent singular vs. plural conventions (e.g., column selectors vs. vector inputs) and documented all intentional exceptions in `docs/conventions/naming-conventions.md`
- [x] Deprecation documentation & warning standardization:
  - Ensured all 27 deprecated files in `R/deprecated-*.R` have consistent roxygen documentation (using `@description \lifecycle{deprecated}`, `@keywords internal`, and exclusion from the main TOC index)
  - Ensured deprecation warnings consistently use `lifecycle::deprecate_warn("0.1.0", ...)` with standard session-level warning frequency and clean delegation to S7 generics
- [x] Cross-family method parity test matrix (in `test-17-S7-model-cross-family-parity.R`)
- [x] Performance regression checks and cache/memory bloat safeguards (in `test-18-S7-model-performance-cache-safeguards.R`)
- [x] Parameter naming harmonization and clean up `component_*` prefix:
  - Harmonized MUVG and MNIX constructors, validators, and representations to use `kappa` and `nu` consistently across families (retiring `component_kappa`, `component_nu`).
  - Removed `untransform_cues` parameter from `get_draws()` and upstream callers (preserving it strictly in `get_exposure_category_statistics()`).
  - Unified `add_category_representation()` parameters across all methods to accept either a scalar for the added category (rescaling existing categories) or a full probability vector.
- [x] Metadata architecture unification between `MVBU_Stanfit` and `MVBU_StanfitPosterior`:
  - Enforced identical `metadata$label_information` schema across both classes, allowing seamless metadata transfer during posterior extraction without separate `groups`/`categories`/`cues` slot branching.
- [x] Documentation tag placement and standardization:
  - Grouped metadata/routing tags (`@docType`, `@name`, `@rdname`, `@usage`, `@export`) systematically at the bottom of roxygen blocks, adhering to standard roxygen2 layout.
  - Eliminated redundant `@name` tags on single generic definitions.
  - Verified and confirmed that all 4 package vignettes build cleanly with 0 errors after parameter renamings.

Phase gate:
- [x] All functions and arguments adhere to documented singular/plural rules
- [x] Deprecated functions cleanly documented with `@keywords internal` and standard lifecycle warnings
- [x] All unit tests (2,314 tests across 40 test files) pass with 0 failures
- [x] Vignette examples run end-to-end without warnings

### Phase 7B: Test Suite Architecture Cleanup
**Goal:** systematically clean and modernize tests folder structure and test code quality

Checklist:
- [x] Reorganize test suite into sequential numeric naming (`test-XX-*.R`) with `S7-` and `deprecated-` prefixes
- [x] Sequence all active tests (`test-00-` to `test-37-`) before all deprecated tests (`test-80-` to `test-90-`)
- [x] Remove obsolete tests and dead aliases (e.g. `sample_observation`)
- [x] Verify all 2,314 tests across the entire suite pass cleanly with 0 failures
- [x] Standardize test inputs via built-in datasets (`pb52`, `mixer6`) and package generators (`example_exposure_test_data()`, `example_ideal_adaptor_stanfit_input()`, `example_ideal_adaptor_staninput()`); keep pre-fitted MCMC models strictly in `tests/testthat/models/`
- [x] Replace ad hoc test setups in legacy helpers with library datasets (`pb52`) and standardized package functions (`write_stanfit()`, `read_stanfit()`)
- [x] Add explicit parity snapshots/checks across NIW/NIX/MNIX/MVG/UVG/MUVG/exemplar families
- [x] Streamline the organization of tests into test files, so that things that are related to each other are tested in the same test file, and so that the name of the test file clearly communicates the content of the test file (grouped under Core `test-00` to `test-07`, Model operations `test-10` to `test-18`, Visualization `test-20` to `test-23`, Stanfit ecosystem `test-30` to `test-37`, and Deprecated `test-80` to `test-90`).

Phase gate:
- [x] New test folder structure documented and used consistently
- [x] All test files pass cleanly (2,314 tests, 0 failures)

### Phase 8A: Documentation, Vignettes, and Extensibility Polish
**Goal:** explain architecture and workflows clearly, modernize primitives, and enable runtime model extension

Checklist:
- [x] Add class-level reference docs for all S7 class families (core classes, representation classes, cognitive model classes)
- [x] v1.0 S7 Class Architecture and Workflows vignette (`vignettes/s7-class-structure-and-workflows.Rmd`)
- [x] Visualizing Models and Categories vignette (`vignettes/visualizing-models-and-categories.Rmd`)
- [x] Fitting and Working with MVBeliefUpdatr Stanfit Models vignette (`vignettes/fitting-and-working-with-stanfit-models.Rmd`)
- [x] Migration vignette with worked examples and API translation table (`vignettes/backward-compatibility-working-with-legacy-code.Rmd`)
- [x] Workflow vignettes (prediction/categorization, plotting, interop)
- [x] Include explicit observer/adaptor pairing map and standalone-family rationale (Exemplar) in architecture-facing docs
- [x] Contributor guide for adding model families vignette (`vignettes/extending-mvbeliefupdatr-adding-model-families.Rmd`)
- [x] Implement `MVBU_Family` descriptor and built-in family registry (density, sampling, moments)
- [x] Centralize `sample_observations()` on representation-level family dispatch (retiring `.mvbu_sample_from_representation()` switch logic)
- [x] Eliminate `.mvbu_point_summary_from_representation()` switch logic in favor of family moments dispatch
- [x] Implement `sample_observations()` for `MVBU_ModelList` allocating total `n` across models and preserving `model`/`group` metadata
- [x] Relocate legacy argument translation (`Ns`, `randomize.order`) exclusively into deprecated `sample_observation()`
- [x] Standardize sufficient category statistics on `x_css` and return `NA_integer_` for representations without sample counts
- [x] Fix uniform likelihood closure bug in extension vignette (`log = log`) restoring 100% response probability plateaus
- [x] Document distributional assumptions (e.g. Gaussian ellipses in `plot_categories()`) for custom category representations
- [x] Does the deprecated `dplyr` verb forwarding described in the backward-compatibility-... vignette actually work? If not implement.

Phase gate:
- [x] Vignette examples run end-to-end (all 5 vignettes render with 0 errors)
- [x] Targeted tests (`test-01`, `test-10`, `test-11`, `test-17`, `test-23`) pass with 0 failures
- [x] All 2,314+ package tests pass with 0 failures

### Phase 8B: Contributor Guide & Permanent Extensibility Infrastructure
**Goal:** author comprehensive guide and infrastructure for adding permanent model families

Checklist:
- [x] Author Developer Guide for permanent model families (`docs/developer/contributor-guide-adding-permanent-model-families.md`)
- [x] Document Stanfit ecosystem integration (`inst/stan/`, C++ headers, data builders, parameter translation tables)
- [x] Document Roxygen Rd standards, `@include` dependencies, and S7 class documentation rules
- [x] Document test fixtures, parameter recovery tests, and parity snapshots
- [x] Formalize representation family property/registration hooks (density, sampling, moments)
- [x] Audit coercion architecture (`.mvbu_coerce_category_representation`) for custom family registration hooks
- [x] Cross-link Developer Guide from the extension vignette

Phase gate:
- [x] Developer guide complete, reviewed, and linked from extension vignette
- [x] Extensibility hooks verified against test fixtures
- [x] Full package test suite passes cleanly with 0 failures

### Phase 9: Stan Expansion Readiness & Prior Refactoring
**Goal:** start Stan additions on stable foundations and unify priors / fixing

Checklist:
- [x] Implement analytical forward-updating functions (`update_template()` /
  `update_category_representation()`) for NIX, MNIX, and NIW models
- [x] Add analytical parity unit tests verifying that 1-cue NIX and 1-cue NIW
  forward updating yield identical belief parameters given identical prior
  beliefs and exposure data
- [x] Define template for new Stan-backed inferred models (documented in
  Developer Guide)
- [x] Begin Stan debugging and model expansion work (verify MNIX and NIX
  multi-cue integration architectures in Stan vs R alignment)
- [x] Add early and late multi-cue integration families (MNIX unweighted /
  weighted / category-weighted; NIX early precision / weighted /
  category-weighted; all 14 Stan models compiling and verified)
- [x] Implement declarative S7 `prior()` constructor in `R/S7-prior.R` with
  Cauchy priors on κ₀, ν₀, τ, Beta on lapse_rate, and Dirichlet on
  category_prior / lapse_bias
- [x] Implement declarative S7 `fix_parameter()` constructor in `R/S7-fixed.R`
  supporting global, per-category, and per-cue specification
- [x] Standardize hyperparameter naming (`prior_*_of_*`) and data-adaptive
  defaults across all 14 Stan models
- [x] Implement partial category fixing via index partitioning in Stan
  (sampling unknown categories while injecting known constants from data)
- [x] Implement 4-case conjugate belief updating in
  `inst/stan/include/belief_updating.stan` and Normal vs. Student-t
  likelihood branching for truly known category moments (μ₀, Σ₀) across all
  14 Stan models
- [x] Update R S7 stanfit input pipelines (`R/S7-staninput.R`) and fitting
  interfaces (`R/S7-stanfit-fitting.R`) to pass known moment data flags
- [x] Document conjugate updating behaviors and parameter fixing in
  `vignettes/fitting-and-working-with-stanfit-models.Rmd` and roxygen Rd files
- [x] Add parameter recovery test suite verifying analytical forward-updated
  beliefs against Stan posterior estimates (recovering prior beliefs given
  sufficient exposure data)

Phase gate:
- [x] All 14 Stan models compile and pass interface contract tests
- [x] Full test suite passes cleanly (over 3,040 tests across 45 test files)

### Phase-End Cleanup Standard (Applies to Every Phase)
Checklist:
- [x] All newly added functions/classes in phase scope have roxygen
  documentation.
- [x] All legacy code integrated into the S7 scaffold in phase scope is
  brought up to roxygen documentation standards.
- [x] Roxygen generation runs for the phase branch without introducing new
  unresolved-link warnings for phase-touched files.
- [x] Generated Rd output for phase-touched topics is checked for malformed
  markup/macros.
- [x] Finalize S7 class registration strategy: retain `package = NULL` for
  clean, unprefixed class names across the package build.

Phase gate:
- [x] Vignette examples run end-to-end
- [x] Migration guide covers all renamed/removed APIs

## Verification Matrix
- [x] Structural validation checks for all classes (`test-01-S7-classes.R`,
  `test-30-S7-stanfit-input.R`, `test-34-S7-stanfit-methods.R`,
  `test-38-S7-prior.R`, `test-39-S7-fixed.R`)
- [x] API parity checks across NIW/MVG/exemplar families
  (`test-17-S7-model-cross-family-parity.R`, `test-10-` through `test-14-`)
- [x] Family-extension checks: adding a new family does not require changing
  core generics (`test-02-S7-family.R`, developer guide and vignette)
- [x] Wrapper compatibility equivalence checks (`test-80-` through `test-90-`,
  `vignettes/backward-compatibility-working-with-legacy-code.Rmd`)
- [x] rstan/tidybayes interoperability checks (`test-31-`, `test-34-`,
  `test-40-`)
- [x] Performance benchmark checks (cached kernels in `R/S7-plot-methods.R`,
  plotting vignette benchmarks)
- [x] Documentation reproducibility checks (all 5 vignettes build end-to-end,
  `devtools::document()` passes cleanly with zero warnings)
- [x] tests folder architecture checks (45 sequentially numbered files, layout,
  fixtures, helper reuse, dead-file cleanup)

## Implementation Anchors
- R/class.R
- R/class-NIW-IA-stanfit.R
- R/class-NIW-IA-staninput.R
- R/make-objects.R
- R/get-info-from-model.R
- R/get-info-from-NIW-IA.R
- R/get-info-from-MVG-IO.R
- R/get-info-from-exemplar-model.R
- R/get-info-from-NIW-IA-stanfit.R
- R/plot-expected-categories.R
- R/methods-NIW-IA-stanfit.R
- R/update-NIW-IA.R
- R/deprecated-2025.R
- R/deprecated-2026.R

## Notes for Refinement
- Keep phases independently shippable where possible
- Attach explicit go/no-go criteria to each phase before coding starts
- If phase gate fails, resolve blockers before entering next phase
- Revisit grouped-model containers: group labels currently identify model instances in model combinations; formal grouped-model class/container design should be addressed in later phases.


# Later planned extensions

+ expand the vignette on fitting and working with stanfit models with a
  case that compares an MVBeliefUpdatr stanfit model against e.g., a logistic
  regression. This is a common use case. This might require thinking about
  how the log posterior densities are stored for MVBeliefUpdatr stanfit models:
  we might have to either give observation level log posterior densities or
  tell users how they'd need to fit their alternative model to have log
  posterior densities for the same observations.

+ Hybrid / Hierarchical Integration Architectures:
  + A subset of cues undergoes early perceptual fusion (e.g. spectral cues F1,
    F2 fused into an intermediate vowel quality percept z), which is then
    integrated late with remaining cues (e.g. duration, f0).
  + This family requires additional Stan models merging the early linear
    projection logic into the NIW or MNIX modeling code.

+ Adaptive Cue Weighting Models (`adaptive_weighted` and
  `adaptive_category_weighted`):
  + Models where attention weights evolve dynamically across exposure trials
    rather than being held static.
  + Theoretical space:
    + Global dynamic weights:
      `adaptive_weighted_integrator_NIX_ideal_adaptor` and
      `MNIX_ideal_adaptor_adaptive_weighted_integrator`.
    + Category-specific dynamic weights:
      `category_adaptive_weighted_integrator_NIX_ideal_adaptor` and
      `MNIX_ideal_adaptor_category_adaptive_weighted_integrator`.
  + Principled formulation of the adaptive function:
    + Error-driven attention shifting: delta-rule updating based on prediction
      error δₜ = yₜ - ŷₜ, where cue weights adjust in the direction reducing
      prediction error, governed by attention learning rate parameter α_w.
    + Attention-gated Bayesian updating: cues with higher informativity /
      lower predictive variance gain attention weight.
  + Stan Requirements:
    + Fitting requires trial-by-trial sequence data rather than batch sufficient
      statistics, and tracks dynamic simplex states across time.
  + To-dos for extension phase:
    + Implement Stan samplers and evaluate identification power and MCMC
      efficiency given the high parameter count.
    + Provide placeholder R functions/classes and documentation describing
      the theoretical formulation.

