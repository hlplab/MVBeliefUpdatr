# Phase 0 API and Class Naming Map (Draft)

Purpose: provide an actionable old-to-new naming proposal for v1 so we can lock naming conventions and generate compatibility wrappers in deprecated-2026.R.

Note: deprecated-2025.R is the legacy historical deprecation file and is out of scope for current migration design.

Additional requirement: naming must scale to additional ideal observer/ideal adaptor families (including MNIX and non-Gaussian representations) without renaming core generics.

## Naming Rules
- New public names favor clear domain semantics over historical abbreviations.
- New internals and public APIs use one naming style consistently.
- Legacy names remain callable only via wrappers in deprecated-2026.R during transition.

## Style Decision (Proposed)
- Public functions: snake_case.
- Classes: Family-prefixed snake/Pascal pattern aligned with Phase 0 contract naming (except that abbreviations like NIW, MVG, MNIX, etc. are visually separated by _).
- Accessors/introspection: keep get_* style uniformly.

## Naming Freeze Status (Phase 0)

Status key:
- Final: approved and ready to enforce in implementation.
- Pending: requires one final decision before freeze.

Decision summary:
- Final: get_* accessor style.
- Final: label terminology (get_*_labels) for current API.
- Final: family template supports future ideal observer/adaptor families.
- Final: family-specific exceptions allowed for clearer non-ideal comparator models.
- Final: plot_categories() as single high-level entry point with type argument.
- Final: add_ibbu_stanfit_draw maps to get_draws (extraction semantics).

Phase 0 naming gate checklist:
- [x] Accessor style selected and documented.
- [x] Family naming template selected and documented.
- [x] Label terminology selected and documented.
- [x] Plot entry-point naming selected and documented.
- [x] Remaining ambiguous legacy mappings resolved.
- [x] Naming map copied into architecture contract as final terms.

## Family Naming Template (Proposed)

Use a reusable template so new families can be added with minimal API churn:

- <Family>_Representation
- <Family>_IdealObserverModel
- <Family>_IdealAdaptorModel
- <Family>_IdealAdaptorFit

Examples:
- NIW_Representation, NIW_IdealAdaptorModel, NIW_IdealAdaptorFit
- MVG_Representation, MVG_IdealObserverModel
- MNIX_Representation, MNIX_IdealObserverModel, MNIX_IdealAdaptorModel, MNIX_IdealAdaptorFit
- Exemplar_Representation, ExemplarModel (or Exemplar_IdealObserverModel if observer/adaptor split is introduced)

## Class Name Map

| Current | Proposed v1 | Notes |
| --- | --- | --- |
| MVG | MVG_Representation | Representation-layer object |
| NIW_belief | NIW_Representation | Representation-layer object |
| exemplars | Exemplar_Representation | Representation-layer object |
| MVG_ideal_observer | MVG_IdealObserverModel | Cognitive model |
| NIW_ideal_adaptor | NIW_IdealAdaptorModel | Cognitive model |
| exemplar_model | ExemplarModel | Cognitive model |
| ideal_adaptor_stanfit | NIW_IdealAdaptorFit | Inferred model |
| ideal_adaptor_staninput | NIW_IdealAdaptorStanInput | Inferred-model support object |
| transform_information | TransformInfo | Shared transform metadata object |

## Constructor Name Map

| Current Pattern | Proposed v1 Pattern | Notes |
| --- | --- | --- |
| make_*_from_data() | new_*_from_data() | Keep from_data suffix for clarity |
| lift_X_to_Y() | as_Y(X) | Prefer as_Y for coercion-like behavior |
| implicit S4 constructor calls | explicit new_*() | Avoid direct slot-level construction |

Examples:
- make_NIW_ideal_adaptor_from_data -> new_niw_ideal_adaptor_model_from_data
- lift_NIW_belief_to_NIW_ideal_adaptor -> as_niw_ideal_adaptor_model
- make_MVG_from_data -> new_mvg_representation_from_data
- future extension example: make_MNIX_ideal_observer_from_data -> new_mnix_ideal_observer_model_from_data

## Core Generic Name Map (High Priority)

| Current | Proposed v1 | Rationale |
| --- | --- | --- |
| get_categorization_from_model | get_categorization | Unified get_* prediction API |
| get_categorization_from_NIW_ideal_adaptor | get_categorization.NIW_IdealAdaptorModel | Method form |
| get_categorization_from_MVG_ideal_observer | get_categorization.MVG_IdealObserverModel | Method form |
| get_categorization_from_exemplar_model | get_categorization.ExemplarModel | Method form |
| get_categorization_function | get_categorization_function | Returns callable function object |
| get_posterior_from_model | get_posterior | Standardized posterior accessor |
| get_draws | get_draws | Keep canonical draw accessor |
| get_params | get_parameters | Standardized introspection |
| get_priors_from_model | get_priors | Standardized introspection |
| get_category_labels_from_model | get_category_labels | Label accessor normalization |
| get_group_levels | get_group_labels | Labels terminology chosen |
| get_cue_levels | get_cue_labels | Labels terminology chosen |
| get_stanfit | get_stanfit | Keep canonical bridge accessor |
| plot_expected_categories_* | plot_categories | One top-level plotting entry point |

## Specialized Extraction and Stats Map

| Current | Proposed v1 | Notes |
| --- | --- | --- |
| get_exposure_category_statistic | get_exposure_statistic | Parameterized statistic extraction |
| get_exposure_category_mean | get_exposure_mean | Thin wrapper to get_exposure_statistic |
| get_exposure_category_cov | get_exposure_covariance | Explicit name |
| get_exposure_category_css | get_exposure_centered_ss | Keep explicit abbreviation expansion |
| get_exposure_category_uss | get_exposure_uncentered_ss | Keep explicit abbreviation expansion |
| get_expected_category_statistic | get_expected_statistic | Symmetric with get_exposure_statistic |

## Update API Map

| Current | Proposed v1 | Notes |
| --- | --- | --- |
| update_NIW_belief_by_one_observation | update_representation | Unified update entry |
| update_NIW_beliefs_incrementally | update_representation_incremental | Explicit batching behavior |
| update_model_decision_bias_by_one_observation | update_model_bias | Model-level update |
| update_model_decision_bias_incrementally | update_model_bias_incremental | Explicit batching behavior |

## Interop and Diagnostics Map

| Current | Proposed v1 | Notes |
| --- | --- | --- |
| summary.ideal_adaptor_stanfit | summary.NIW_IdealAdaptorFit | Same behavior under new class |
| loo.ideal_adaptor_stanfit | loo.NIW_IdealAdaptorFit | Forward or bridge to stanfit |
| add_ibbu_stanfit_draw | get_draws | Deprecated wrapper already forwards to get_draws; extraction semantics only |

## Forward Compatibility Requirements

- Adding a new family (e.g., MNIX or non-Gaussian alternatives) must not require renaming any core generics.
- New families should only need:
	- one representation class
	- one or more cognitive model classes (observer/adaptor as appropriate)
	- inferred fit class (if fitted via Stan)
	- family-specific methods on existing generics
- Family-specific helper functions may exist, but public high-level workflows should remain generic-driven.

## Compatibility Wrapper Plan

All legacy names in this map should receive wrappers in deprecated-2026.R with:
- deprecation warning
- direct forwarding to new v1 name
- argument name translation where needed
- no business logic

## Open Naming Decisions to Resolve

No remaining naming decisions are currently open in this map.

## Immediate Use
- Naming gate checklist is complete for this map.
- Freeze final terms in phase0-architecture-contract.md.
- Use final map to generate deprecated-2026.R wrappers during Phase 4.
