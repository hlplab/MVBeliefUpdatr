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
| 6 Performance and Caching Framework | In Progress | Open | Cached posterior closures implemented in plot engine; performance benchmarks in plotting vignette. |
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
- [ ] Object-size regression checks pass
- [x] Print snapshot tests pass

### Phase 6: Performance and Caching Framework
**Goal:** fix high-impact bottlenecks with controlled memory use

Planning note: see [docs/planning/phase6-performance-caching-framework.md](docs/planning/phase6-performance-caching-framework.md) for the detailed design and implementation outline for the posterior-kernel caching strategy.

Checklist:
- [x] Implement operation-aware caching for decision rules and posterior closures in plotting and categorization pipelines
- [ ] Add cache keys/invalidation (version, transform settings, draw filters, prediction options)
- [ ] Keep large plot-grid computations in dedicated grid/result objects

Priority hotspots to benchmark:
- [x] density plotting recomputation loops (optimized with cached closures)
- [ ] repeated draw extraction paths
- [ ] array-to-tibble loop conversions
- [x] categorization function generation pipelines (optimized with cached vectorized closures)

Phase gate:
- [ ] Benchmark thresholds met
- [ ] Memory ceiling checks pass

### Phase 7: Consistency and Quality Hardening
**Goal:** lock reliability and prevent regressions

Checklist:
- [ ] Add cross-family method parity test matrix
- [ ] Add compatibility and deprecation regression tests
- [ ] Add performance regression suite
- [ ] Add cache/memory bloat safeguards

Phase gate:
- [ ] CI matrix green
- [ ] No critical regressions vs baseline

### Phase 7B: Test Suite Architecture Cleanup
**Goal:** systematically clean and modernize tests folder structure and test code quality

Checklist:
- [x] Reorganize test suite into sequential numeric naming (`test-XX-*.R`) with `S7-` and `deprecated-` prefixes
- [x] Sequence all active tests (`test-01-` to `test-26-`) before all deprecated tests (`test-80-` to `test-88-`)
- [x] Remove obsolete tests and dead aliases (e.g. `sample_observation`)
- [x] Verify all 1,414 tests across the entire suite pass cleanly with 0 failures
- [ ] Introduce a systematic `tests/testthat/data/` layout for reusable fixtures and generated test inputs
- [ ] Replace remaining ad hoc test setups in legacy helpers with shared fixtures
- [ ] Add explicit parity snapshots/checks across NIW/MVG/exemplar families

Phase gate:
- [x] New test folder structure documented and used consistently
- [x] All 35 test files pass cleanly (1,414 tests, 0 failures, 1 skip)

### Phase 8: Documentation and Vignettes
**Goal:** explain architecture and workflows clearly

Checklist:
- [x] Add class-level reference docs for all S7 class families (core classes, representation classes, cognitive model classes classes)
- [ ] Include explicit observer/adaptor pairing map and standalone-family rationale (Exemplar) in architecture-facing docs
- [ ] v1.0 architecture vignette
- [ ] migration vignette with worked examples
- [ ] workflow vignettes (fitting, prediction/categorization, plotting, interop)
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