## Current Status

Phase 3 is completed and closed. The S7 foundation, model constructors, family coercions, accessors (`get_model_type()`, `get_representation_type()`, `get_noise_treatment()`, `get_lapse_treatment()`), high-level plotting engine (`plot_categories()`, `plot_categorization_function()`), model aggregation (`aggregate_models()`), and Stanfit/Staninput classes are fully implemented and verified. The test suite has been systematically reorganized into sequentially numbered files (`test-01-` through `test-28-` for modern/S7 tests, and `test-80-` through `test-88-` for deprecated tests) with all 1,494 tests passing cleanly.

## Implementation checklist for the S7 Stanfit migration

- [x] Add `MVBU_StanInput` as the base class
- [x] Add family-specific subclasses such as `NIW_IdealAdaptorStanInput`, `NIX_IdealAdaptorStanInput`, and others only as needed
- [x] Add `MVBU_StanFit` as the base fit-result class
- [x] Add model-specific Stan-fit subclasses following the same naming scheme if they are needed
- [x] Add validators for the new classes and enforce the required slots

### Phase 2: Split the Stan input and fit pipeline
- [x] Rename the current Stan-input preparation file `make-staninput.R` to `R/S7-staninput.R`
- [x] Update the input-preparation code to construct the new S7 Stan-input object instead of returning the old list-based structure
- [x] Update the fit pipeline so it constructs a new Stan-fit object and attaches it to the relevant model-distribution object
- [x] Preserve transformed/untransformed input and transform metadata in the new fit object

### Phase 3: Move Stan accessors and helpers to S7
- [x] Rename the core class file to `R/S7-core-classes.R`
- [x] Rename the generics file to `R/S7-generics.R`
- [x] Rename the methods file to `R/S7-core-methods.R`
- [x] Add `get_model_type()` as a generic and implement it for standard model objects and S7 objects
- [x] Add `get_representation_type()` as a generic and implement it for representation objects
- [x] Implement Stan-facing accessors in the S7 methods layer:
  - [x] `get_staninput()`
  - [x] `get_stanfit()`
  - [x] `set_stanfit()`
  - [x] `get_transform_information()`
  - [x] `get_transform_function()`
  - [x] `get_untransform_function()`
- [x] Add `get_noise_treatment()` and `get_lapse_treatment()` generics and methods across S7 hierarchy
- [x] Remove obsolete `sample_observation` alias, keeping `sample_observations()`
- [x] Retire `wide` option from `get_draws()` generic and methods
- [x] Add forwarding methods for standard rstan-style access patterns where appropriate

### Completed adjacent S7 migration work
- [x] Add S7 constructors from data for representations, templates, and models across UVG, NIX, MUVG, MNIX, MVG, NIW, and EXEMPLAR.
- [x] Add generic `type` dispatchers for the from-data constructors.
- [x] Add family coercions for representations, templates, and models, including multi-cue MUVG/MNIX support.
- [x] Move examples to `R/S7-example-objects.R` and provide representation, template, model, and generic example functions.
- [x] Replace legacy sampling helpers with `sample_observations()` for S7 objects.
- [x] Reduce `make_*`/`lift_*` functions to deprecated wrappers.
- [x] Consolidate modern constructor, coercion, example, and deprecated-wrapper tests.
- [x] Unified plotting engine: `plot_categories()`, `plot_categorization_function()`, `plot_parameters()`, `plot_cue_correlations()`, `plot_cue_densities()` across 1D, 2D, 3D, sliced, and interactive plotly modes.
- [x] Test suite architecture cleanup: renumbered and sequenced tests (active `test-01-` to `test-28-`, deprecated `test-80-` to `test-88-`).

### Phase 4: Consolidate helpers and utilities
- [x] Create `R/S7-stanfit-utils.R` for general stanfit-related helper functions
- [x] Rename helper functions to neutral names where appropriate:
  - [x] `read_stanfit()` instead of `read_ideal_adaptor_stanfit()`
  - [x] `write_stanfit()` instead of `write_ideal_adaptor_stanfit()`
- [x] Move caching, persistence, and refit logic into the new utility file (removed from `R/S7-stanfit.R`, added unit tests in `test-28-S7-stanfit-utils.R`)

### Phase 5: Remove obsolete wrappers and fix broken code
- [x] Remove the old S3/S4 wrapper layer for Stanfit-related class types
- [x] Remove temporary compatibility shims that were only introduced to keep the package loadable during migration
- [x] Delete or simplify any now-unnecessary helper functions such as old `is.ideal_adaptor_stanfit_input`-style shims if they are redundant with S7 validators or new class checks
- [x] Fix any broken or convoluted logic introduced during earlier migration attempts, especially around transform handling, fit-object persistence, Stan input validation, and draws detection

### Testing plan
- [x] Add regression tests for the new Stan-input classes
- [x] Add regression tests for the new Stan-fit classes
- [x] Add regression tests for `get_model_type()` and `get_representation_type()`
- [x] Add regression tests for `get_noise_treatment()`, `get_lapse_treatment()`, and `get_draws(wide = )` deprecation
- [x] Standardize test suite naming and execution order (test-01 through test-28 non-deprecated, test-80 through test-88 deprecated)
- [x] Verify the package still loads and all 1,494 tests pass cleanly

## Completed Work & Phase Gate Exit
- [x] Author S7 architecture & workflow vignette introducing S7 class structure, bare/from_data constructors, print/summarize/as_tibble, Stanfit classes, `example_*` functions, model aggregation, and the legacy coercion bridge.
- [x] Separate legacy tibble bridge workflow into dedicated vignette (`backward-compatibility-working-with-old-code.Rmd`).
- [x] Migrate likelihood, categorization, evaluation, and legacy info consumers.
- [x] Clean up obsolete shims (`R/s7-phase2-migration.R`), deprecate `is.ideal_adaptor_stanfit()`, and verify non-deprecated test suite.

