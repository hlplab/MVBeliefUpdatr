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
| 7 Consistency and Quality Hardening | Not Started | Open | Follows implementation phases. |
| 7B Test Suite Architecture Cleanup | Completed | Closed | Reorganized into sequential numbering: active tests (`test-01-` to `test-26-`), deprecated tests (`test-80-` to `test-88-`) running strictly after active tests. 1,414 tests passing. |
| 8 Documentation and Vignettes | In Progress | Open | Comprehensive plotting vignette completed; S7 architecture & workflow vignette planned next. |
| 9 Stan Expansion Readiness | Not Started | Open | Deferred until architecture and API stabilize. |

### Locked Decisions
- [x] Release strategy: big-bang v1.0
- [x] Naming direction: aggressive + aliases in deprecated-2026.R
- [x] Accessor/introspection syntax: prefer get_* for consistency
- [x] Backward wrappers: keep for 2-3 minor releases
- [x] Core composition: explicit model -> representation
- [x] Stan timeline: after class/API stabilization
- [ ] Finalize cache profile defaults (minimal/standard/eager)
- [x] Documentation policy: all new code and all code integrated into the new S7 scaffold must be roxygen documented
- [x] Documentation QA policy: each phase-end cleanup includes checks for broken links and roxygen/Rd issues

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
- [ ] Singular / plural naming audit:
  - Systematically review all exported functions and argument names for consistent singular vs. plural conventions (e.g., column selectors vs. vector inputs) and document all intentional exceptions
- [ ] Test helper simplification:
  - Reduce and simplify test helpers (e.g. evaluating whether `helper-vowel-data.R` / `make_vowel_test_data` can be replaced by standard fixtures while maintaining full coverage)
- [ ] Deprecation documentation & warning standardization:
  - Ensure all deprecated functions have consistent roxygen documentation (using `@description \lifecycle{deprecated}`, `@keywords internal`, and exclusion from the main TOC index)
  - Ensure deprecation warnings consistently use `lifecycle::deprecate_warn("0.1.0", ...)` with standard session-level warning frequency (and active warning signaling during tests)
- [ ] Vignette enhancement:
  - Add explicit worked example demonstrating the output of the 3 decision rules (`"criterion"`, `"proportional"`, `"sampling"`) in the vignette covering `categorize()`
- [ ] Cross-family method parity test matrix
- [ ] Performance regression checks and cache/memory bloat safeguards

Phase gate:
- [ ] All functions and arguments adhere to documented singular/plural rules
- [ ] Deprecated functions cleanly documented with `@keywords internal` and standard lifecycle warnings
- [ ] All unit tests (1,680+) pass with 0 failures and 0 unexpected warnings
- [ ] Vignette examples run end-to-end without warnings

### Phase 7B: Test Suite Architecture Cleanup
**Goal:** systematically clean and modernize tests folder structure and test code quality

Checklist:
- [x] Reorganize test suite into sequential numeric naming (`test-XX-*.R`) with `S7-` and `deprecated-` prefixes
- [x] Sequence all active tests (`test-01-` to `test-26-`) before all deprecated tests (`test-80-` to `test-88-`)
- [x] Remove obsolete tests and dead aliases (e.g. `sample_observation`)
- [x] Verify all 1,414 tests across the entire suite pass cleanly with 0 failures
- [ ] Introduce a systematic `tests/testthat/data/` layout for reusable fixtures and generated test inputs
- [ ] Replace remaining ad hoc test setups in legacy helpers with shared fixtures; reduce ad-hoc code (use code from package instead)
- [ ] Add explicit parity snapshots/checks across NIW/MVG/exemplar families

Phase gate:
- [x] New test folder structure documented and used consistently
- [x] All 35 test files pass cleanly (1,414 tests, 0 failures, 1 skip)

### Phase 8: Documentation and Vignettes
**Goal:** explain architecture and workflows clearly

Checklist:
- [x] Add class-level reference docs for all S7 class families (core classes, representation classes, cognitive model classes classes)
- [ ] Include explicit observer/adaptor pairing map and standalone-family rationale (Exemplar) in architecture-facing docs
- [x] v1.0 S7 Class Architecture and Workflows vignette (`vignettes/s7-class-structure-and-workflows.Rmd`)
- [x] Visualizing Models and Categories vignette (`vignettes/visualizing-models-and-categories.Rmd`)
- [x] Fitting and Working with MVBeliefUpdatr Stanfit Models vignette (`vignettes/fitting-and-working-with-stanfit-models.Rmd`)
- [ ] migration vignette with worked examples
- [ ] workflow vignettes (prediction/categorization, interop)
- [ ] contributor guide for adding model families

### Phase-End Cleanup Standard (Applies to Every Phase)
Checklist:
- [ ] All newly added functions/classes in phase scope have roxygen documentation.
- [ ] All legacy code integrated into the S7 scaffold in phase scope is brought up to roxygen documentation standards.
- [ ] Roxygen generation runs for the phase branch without introducing new unresolved-link warnings for phase-touched files.
- [ ] Generated Rd output for phase-touched topics is checked for malformed markup/macros.
- [ ] Remove temporary S7 scaffold-only `package = NULL` overrides once the class registration strategy is finalized for the packaged build.

Phase gate:
- [ ] Vignette examples run end-to-end
- [ ] Migration guide covers all renamed/removed APIs

### Phase 9: Stan Expansion Readiness
**Goal:** start Stan additions on stable foundations

Checklist:
- [ ] Implement analytical forward-updating functions (`forward_update()`) for NIX, MNIX, and NIW models
- [ ] Add parameter recovery test suite verifying analytical forward-updated beliefs against Stan posterior estimates (recovering prior beliefs given sufficient exposure data)
- [ ] Define template for new Stan-backed inferred models
- [ ] Reuse diagnostics/extraction interfaces across Stan families
- [ ] Begin Stan debugging and model expansion work
- [ ] Add at least one additional family plan item (e.g., MOG or MNIX ideal observer/adaptor) using the same extension template

Phase gate:
- [ ] First new Stan-family prototype passes interface contract tests

## Verification Matrix
- [ ] Structural validation checks for all classes
- [ ] API parity checks across NIW/MVG/exemplar families
- [ ] Family-extension checks: adding a new family does not require changing core generics
- [ ] Wrapper compatibility equivalence checks
- [ ] rstan/tidybayes interoperability checks
- [ ] Performance benchmark checks
- [ ] Documentation reproducibility checks
- [ ] tests folder architecture checks (layout, fixtures, helper reuse, dead-file cleanup)

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

## Later planned extensions

### Issues to fix
+ Check MNIX cue weighting: there are at least two possible implementations:
  + adapt after integration (single nix sitting on top of integrated cue dimension)
  + adapt before integration (one nix for each cue; then integration)
  + make sure that whatever is implemented for Stan is also what is implemented in R
  + make sure that these choices are documented verbosely in roxygen
+ MNIX fitting is currently commented out. AI comments for MNIX failure:
  > "The root cause is visible now: the MNIX Stan program expects a covariance-style summary array, but the current builder is providing a sum-of-squares vector array instead. I’m aligning that data structure with the model’s declared interface before I verify again."

### Extensions
+ handling of tau_scale and hyper-priors more generally is really non-transparent atm. rather than handing known mu, sigma, switch to allowing specification of informative priors for all parameters. for m, S, etc. create tools that translate the prior from an intuitive space (mu, Sigma) to the relevant underlying parameter space (m, S, nu, kappa).
+ Inverse MUVG/MNIX models (both as Stanfit model and S7 representation/template/model): a model that accepts multiple cues as input but maintains and updates representations/templates/models over the single integrated cue dimension (e.g. UVG-I / NIX-I).
+ Write forward-updating functions for NIX, MNIX, NIW, revising existing NIW forward-updating.
  + Include equivalence checks verifying that 1-cue NIX and 1-cue NIW forward updating yield identical analytical belief parameters given identical prior beliefs and exposure data.
+ Add specification of category prior and lapse bias via `fixed_parameters` and extend Stan program to include category prior.
+ Expand Stan programs to allow users to specify prior $m, S$ for each category (not fixed point estimates; evaluate if/how this differs from handing $\mu, \Sigma$ and inferring $\kappa, \nu$).
+ Refine the pre-compiled 2D example `stanfit` object so category distributions have less overlap, creating a clearer categorization surface.
+ For stanfit fitting, allow specification of fixed parameters
  + for *some* mu and some sigma. that will require changes to the stan code (and might require making multiple versions of each stan model to maintain efficiency for the most common case in which none or all of the mu, sigma's are fixed.)
  + category priors and lapse bias 
+ make extensions of plot_categories and plot_categorization_functions functions that animate model updates (both for histories of update_model or for stanfit ideal adaptors by stepping from prior to posterior for a number of equally-spaced draws)

### Efficiency considerations
+ Consider making separate versions of stan models for 0, 1, or more observations. Functions could be shared between them to ease maintenance.