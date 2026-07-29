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
| 3 API Unification and Method Coverage | In Progress | Open | Initial S7 unification landed for categorization/prediction/posterior-prediction with consistent single and list/batch handling. |
| 4 Compatibility Shell | Not Started | Open | Can start late in Phase 3. |
| 5 Data Model and Print Strategy | Not Started | Open | Starts after class migration baseline. |
| 6 Performance and Caching Framework | Not Started | Open | Starts after API baseline stabilizes. |
| 7 Consistency and Quality Hardening | Not Started | Open | Follows implementation phases. |
| 7B Test Suite Architecture Cleanup | Not Started | Open | Scheduled later in process by request. |
| 8 Documentation and Vignettes | Not Started | Open | Runs continuously, final hardening late. |
| 9 Stan Expansion Readiness | Not Started | Open | Deferred until architecture and API stabilize. |

### Locked Decisions
- [x] Release strategy: big-bang v1.0
- [x] Naming direction: aggressive + aliases in deprecated-2026.R
- [x] Accessor/introspection syntax: prefer get_* for consistency
- [x] Backward wrappers: keep for 2-3 minor releases
- [x] Core composition: explicit model -> representation
- [x] Stan timeline: after class/API stabilization
- [ ] Finalize ModelDistribution interop strategy (A/B/C below)
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
  - ModelDistribution: inferred posterior/distribution-level package object associated with a model family

### Decision Blocks

#### Decision: ModelDistribution Interop
- Option A: strict wrapper with get_stanfit()/as_stanfit()
- Option B: wrapper + S3 forwarding for key rstan/tidybayes generics
- Option C: on-demand adapter object for external tooling
- Recommendation: Option B
- Final choice: [ ] A  [x] B  [ ] C
- Policy: align with rstan/tidybayes where relevant, but avoid adding extra dependencies unless clearly justified.
- Tracking: maintain an interop bridge watchlist and finalize exact bridge methods incrementally during Phases 1-3.

#### Decision: Cache Semantics
- Global default: pure-by-default methods
- Persistent cache candidates (inside ModelDistribution objects):
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
- [x] Finalize class family contract: Representation, CognitiveModel, ModelDistribution
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
- [x] Implement abstract S7 base classes (Representation, CognitiveModel, ModelDistribution)
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
- [x] Migrate ModelDistribution classes around stanfit outputs
- [x] Standardize constructors and defaults
- [x] Validate migration pattern can be reused by at least one future-family prototype (MUVG, MNIX, or another non-Gaussian family)

Phase gate:
- [x] Constructor normalization tests pass
- [x] No mixed S4/S7 construction paths for migrated classes

### Phase 3: API Unification and Method Coverage
**Goal:** normalize user-facing behavior and signatures

Checklist:
- [ ] Replace mixed S3/S4/direct dispatch with unified S7 API surface. 
- [ ] Update handling of stanfit-related classes. For details, see [docs/planning/phase3-stanfit-migration.md](docs/planning/phase3-stanfit-migration.md).
- [x] Normalize categorization/prediction signatures across model families
- [x] Support both single and list/batch forms consistently
- [ ] Unify high-level plotting entry points with internal dimension specialization
- [ ] Fill method coverage gaps (summary/print/plot/getters/update parity)
- [ ] Ensure generic signatures remain family-agnostic for Gaussian and non-Gaussian model families
- [ ] Close interop bridge watchlist items based on workflow tests and dependency budget

Phase gate:
- [ ] Cross-family API parity tests pass
- [ ] Signature consistency tests pass

### Phase 4: Compatibility Shell
**Goal:** keep old API available but isolated

Checklist:
- [ ] Create deprecated-2026.R and move all wrappers there
- [ ] Keep deprecated-2025.R untouched except for archival maintenance
- [ ] Add deprecation warnings and migration hints
- [ ] Ensure wrappers are thin adapters only
- [ ] Add migration mapping table old -> new API

Phase gate:
- [ ] Wrapper output-equivalence tests pass
- [ ] Deprecation warnings fire as expected

### Phase 5: Data Model and Print Strategy
**Goal:** remove tibble-as-object while preserving readable printing

Checklist:
- [ ] Store canonical fields once (no constant duplication across rows)
- [ ] Keep human-readable tibble-style print methods
- [ ] Add conversion helpers (to_tibble(), as_draws_df(), as_matrix where needed)

Phase gate:
- [ ] Object-size regression checks pass
- [ ] Print snapshot tests pass

### Phase 6: Performance and Caching Framework
**Goal:** fix high-impact bottlenecks with controlled memory use

Planning note: see [docs/planning/phase6-performance-caching-framework.md](docs/planning/phase6-performance-caching-framework.md) for the detailed design and implementation outline for the posterior-kernel caching strategy.

Checklist:
- [ ] Implement operation-aware caching
- [ ] Add cache keys/invalidation (version, transform settings, draw filters, prediction options)
- [ ] Keep large plot-grid computations in dedicated grid/result objects
- [ ] Add optional eager postprocess profile for ModelDistribution objects with size guards

Priority hotspots to benchmark:
- [ ] density plotting recomputation loops
- [ ] repeated draw extraction paths
- [ ] array-to-tibble loop conversions
- [ ] categorization function generation pipelines

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

### Phase 7B: Test Suite Architecture Cleanup (later in process)
**Goal:** systematically clean and modernize tests folder structure and test code quality

Checklist:
- [ ] Reorganize tests into a clear, family-based structure (representations, models, ModelDistribution objects, interoperability)
- [ ] Introduce a systematic tests/testthat/data/ layout for reusable fixtures and generated test inputs
- [ ] Replace ad hoc or repetitive test setup with shared helpers/fixtures
- [ ] Rewrite fragile or outdated existing tests to current API and naming conventions
- [ ] Remove obsolete test files and dead test helpers after migration coverage is in place
- [ ] Add explicit parity snapshots/checks to ensure equivalent behavior across NIW/MVG/exemplar families

Phase gate:
- [ ] New test folder structure documented and used consistently
- [ ] Legacy test paths removed or justified temporarily
- [ ] Updated tests are readable, deterministic, and pass in CI

### Phase 8: Documentation and Vignettes
**Goal:** explain architecture and workflows clearly

Checklist:
- [ ] Add class-level reference docs for all S7 class families (core classes, representation classes, cognitive model classes, and model distribution classes)
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
- [ ] API parity checks across NIW/MVG/exemplar/ModelDistribution families
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