# Phase 0 Architecture Contract: MVBeliefUpdatr v1.0

## Status
Near-final draft pending final sign-off on 3 decision blocks.
This contract becomes binding for Phases 1+ once the three blocks in Section 14 are resolved.

## 1) Architectural Principles
- All core user-facing entities are model objects.
- Model objects explicitly contain representation objects.
- Inferred-model objects are first-class package objects, not raw stanfit objects.
- Interop with rstan/tidybayes is preserved through explicit bridges and selected method forwarding.
- API consistency is prioritized over backward shape compatibility.
- Old API compatibility for the v1 migration lives in deprecated-2026.R.
- Legacy historical deprecations in deprecated-2025.R are treated as out of scope for current migration design.

## 2) Canonical Class Hierarchy (S7)

### Abstract base classes
- MVBU_Object
- MVBU_Representation (extends MVBU_Object)
- MVBU_CognitiveModel (extends MVBU_Object)
- MVBU_InferredModel (extends MVBU_Object)

### Representation families
- NIW_Representation (extends MVBU_Representation)
- MVG_Representation (extends MVBU_Representation)
- Exemplar_Representation

### Cognitive model families
- NIW_IdealAdaptorModel (extends MVBU_CognitiveModel)
- MVG_IdealObserverModel (extends MVBU_CognitiveModel)
- ExemplarModel 

### Inferred model families
- NIW_IdealAdaptorFit (extends MVBU_InferredModel)
- Future families should generally follow <Family>_IdealObserverModel, <Family>_IdealAdaptorModel, and <Family>_IdealAdaptorFit when those concepts apply; family-specific exceptions are allowed for clearer non-ideal comparator models.

## 3) Composition Contract

Each MVBU_CognitiveModel must contain:
- multiple representations, one for each category: an MVBU_Representation subclass instance
- decision_rule: scalar character
- lapse_rate: scalar [0,1]
- lapse_bias: structured numeric value (vector/matrix as needed). one bias per category, each in [0,1], summing to 1
- perceptual_noise: structured numeric value (covariance/scalar as needed)
- priors: optional prior bundle object, one prior per category, each in [0,1], summing to 1
- metadata: optional list for provenance and labels (cue labels, category labels, etc. for ease of access through other functions)

Each MVBU_InferredModel must contain:
- model_family: identifier
- stanfit_ref: raw stanfit object
- staninput_ref: typed staninput object (or validated list)
- data_ref: source data handle/object
- transform_info: transform descriptors/functions
- cache: optional cache container (nullable)
- metadata: provenance/version/schema tags

## 4) Interop Contract for Inferred Models

Default strategy: Wrapper + forwarding (approved direction).

### Required bridge methods
- get_stanfit(x)
- as_stanfit(x)
- get_draws(x, ...)

### Required compatibility surface (initial)
- summary(x, ...)
- print(x, ...)
- loo(x, ...)
- posterior-related extraction methods currently used in package workflows
- tidybayes draw extraction pathways used by current plotting/posterior code

### Constraint
- Forwarding methods must not duplicate heavy computations unnecessarily.
- If forwarding requires transformations, cache keying/invalidation rules must apply.

## 5) Canonical Generic Names (User-facing)

### Construction and validation
- new_model(...)
- new_representation(...)
- validate_object(x)
- is_valid(x)

### Introspection
- get_model_family(x)
- get_representation(x)
- get_parameters(x)
- get_priors(x)
- get_cue_labels(x)
- get_category_labels(x)
- get_group_labels(x)

### Prediction/categorization
- get_categorization(x, new_data, ...)
- get_category_prediction(x, new_data, ...)
- get_posterior_prediction(x, new_data, ...)

### Updating
- update_model(x, data, ...)
- update_representation(x, data, ...)

### Extraction
- get_draws(x, ...)
- get_posterior(x, ...)
- get_expected_category(x, ...)

### Visualization
- plot_categories(x, ...)
- plot_parameters(x, ...)
- plot_diagnostics(x, ...)

### Conversion
- as_tibble(x, ...)
- as_matrix(x, ...)
- as_draws_df(x, ...)

## 6) Naming and Migration Rules

- New implementation uses canonical names only.
- Legacy names are wrappers in deprecated-2026.R only.
- Legacy wrappers already present in deprecated-2025.R are not part of the new compatibility surface.
- Wrappers emit deprecation warnings with migration hints.
- No legacy business logic in wrappers.
- No mixed naming styles inside new class files.

## 7) Method Parity Requirements

For all cognitive model families (NIW/MVG/Exemplar), Phase 3 minimum parity:
- construct
- validate
- introspection getters (family/labels/parameters)
- categorize/predict
- update
- summary/print
- high-level category plotting

For inferred model families, minimum parity:
- bridge to stanfit
- draws extraction
- posterior summaries
- diagnostics plotting/extraction used by workflows

## 8) Data and Storage Contract

- Core objects store canonical values exactly once.
- Do not duplicate constant model parameters row-wise.
- Human-readable tibble representations are produced by converters/print methods.
- Large intermediate plot grids/surfaces are stored in dedicated result objects, not in core model objects.

## 9) Cache Contract

Default behavior:
- Pure-by-default for core methods.
- Optional opt-in cache writes through arguments and global options.

Cache profiles:
- minimal: little/no persistence
- standard: cache common repeated inferred outputs
- eager: precompute broader inferred outputs with memory guards

Cache validity must depend on:
- model schema/version
- transform settings
- draw selection/filtering
- prediction and categorization options

## 10) Testing and Documentation Contract (Per Phase)

Each phase must include:
- testthat updates for all touched interfaces
- method parity tests for impacted families
- documentation updates for all changed public behavior
- migration notes for renamed/removed interfaces

## 11) Cleanup Contract (Per Phase)

Each phase must remove outdated implementation paths when superseded.
Allowed exception: temporary compatibility wrappers in deprecated-2026.R.
Legacy file deprecated-2025.R is retained for historical compatibility only.
When removing obsolete code:
- remove dead functions
- remove dead files
- remove stale references from docs and tests
- keep deprecated surface minimal and explicit

## 12) Release Hygiene Contract

If a phase introduces meaningful package changes, update:
- NEWS.md with concise user-facing changelog entries
- DESCRIPTION metadata and version field as triggered by change scope

Versioning guidance:
- patch: bugfix only
- minor: additive/deprecations
- major: breaking redesign

## 13) Exit Criteria for Phase 0

All must be approved:
- class hierarchy names
- generic naming set
- interop strategy and initial forwarded surface
- cache defaults/profile policy
- parity requirements
- per-phase tests/docs policy
- cleanup policy
- release hygiene policy

## 14) Open Items for Immediate Resolution

### Decision Block A: Constructor Strategy
Scope:
- Whether constructors are family-specific only (e.g., new_niw_ideal_adaptor_model_from_data())
- Or whether to also provide a unified new_model() entry point.

Current recommendation:
- Keep family-specific constructors as canonical; optionally provide new_model() as a thin router only if it remains unambiguous.

Sign-off choice:
- [x] A1 Family-specific only
- [ ] A2 Family-specific + unified router

### Decision Block B: Interop Forwarded Surface (NIW_IdealAdaptorFit)
Scope:
- Finalize exactly which rstan/tidybayes-facing methods are forwarded or bridged in Phase 1-3.

Policy:
- Prefer alignment with rstan/tidybayes interfaces used in real workflows.
- Avoid introducing new dependencies beyond current accepted standards unless they provide clear net benefit.
- Finalize bridge coverage iteratively during implementation, tracked in an interop watchlist.

Minimum required surface (already agreed in principle):
- summary(x, ...)
- print(x, ...)
- loo(x, ...)
- get_draws(x, ...)
- get_stanfit(x)
- as_stanfit(x)

Sign-off additions (choose any required now):
- [x] B1 posterior::as_draws_df bridge method(s)
- [x] B2 additional tidybayes extraction helpers currently used in plotting pipeline
- [ ] B3 no additional methods for Phase 1-3

Operationalization note:
- Exact tidybayes/rstan bridge set is intentionally finalized incrementally in Phase 1-3, based on failing tests, workflow needs, and dependency cost/benefit.

### Decision Block C: Cache Profile Defaults
Scope:
- Confirm default cache profile and mutability policy for inferred-model workflows.

Already agreed:
- Pure-by-default behavior
- Per-call overrides allowed
- Large grid/surface results are external result objects, not stored in core model objects

Remaining sign-off choice:
- [ ] C1 Default profile = minimal
- [x] C2 Default profile = standard
- [ ] C3 Default profile = eager

Optional policy detail:
- [x] C4 Enable package option for session-wide default override
- [ ] C5 Keep default fixed and require explicit per-call overrides

### Near-Final Freeze Note
All other naming and architectural decisions in this contract are considered frozen for Phase 0 unless a blocker is discovered.

## 15) Non-goals (for early phases)
- Stan model expansion
- Advanced optimization beyond high-impact known bottlenecks
- Full redesign of all legacy helper utilities not on active workflows