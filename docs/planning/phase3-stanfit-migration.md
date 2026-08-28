## Current Status

Phase 3 is in progress. The S7 foundation, model constructors, family coercions, and modern example/test infrastructure are substantially migrated. The remaining work is concentrated in NIW updating, higher-level legacy consumers, and the Stanfit/model-distribution integration.

## Implementation checklist for the S7 Stanfit migration

- [x] Add `MVBU_StanInput` as the base class
- [x] Add family-specific subclasses such as `NIW_IdealAdaptorStanInput`, `NIX_IdealAdaptorStanInput`, and others only as needed
- [x] Add `MVBU_StanFit` as the base fit-result class
- [x] Add model-specific Stan-fit subclasses following the same naming scheme if they are needed
- [ ] Extend `MVBU_ModelDistribution` so it can hold a Stan fit object and fit metadata
- [x] Add validators for the new classes and enforce the required slots

### Phase 2: Split the Stan input and fit pipeline
- [x] Rename the current Stan-input preparation file `make-staninput.R` to `R/S7-staninput.R`
- [x] Update the input-preparation code to construct the new S7 Stan-input object instead of returning the old list-based structure
- [ ] Update the fit pipeline so it constructs a new Stan-fit object and attaches it to the relevant model-distribution object
- [x] Preserve transformed/untransformed input and transform metadata in the new fit object

### Phase 3: Move Stan accessors and helpers to S7
- [x] Rename the core class file to `R/S7-core-classes.R`
- [x] Rename the generics file to `R/S7-generics.R`
- [x] Rename the methods file to `R/S7-core-methods.R`
- [ ] Add `get_model_type()` as a generic and implement it for standard model objects and S7 objects
- [ ] Add `get_representation_type()` as a generic and implement it for representation objects
- [x] Implement Stan-facing accessors in the S7 methods layer:
  - [x] `get_staninput()`
  - [x] `get_stanfit()`
  - [x] `set_stanfit()`
  - [x] `get_transform_information()`
  - [x] `get_transform_function()`
  - [x] `get_untransform_function()`
- [ ] Add forwarding methods for standard rstan-style access patterns where appropriate

### Completed adjacent S7 migration work
- [x] Add S7 constructors from data for representations, templates, and models across UVG, NIX, MUVG, MNIX, MVG, NIW, and EXEMPLAR.
- [x] Add generic `type` dispatchers for the from-data constructors.
- [x] Add family coercions for representations, templates, and models, including multi-cue MUVG/MNIX support.
- [x] Move examples to `R/S7-example-objects.R` and provide representation, template, model, and generic example functions.
- [x] Replace legacy sampling helpers with `sample_observations()` for S7 objects.
- [x] Reduce `make_*`/`lift_*` functions to deprecated wrappers.
- [x] Consolidate modern constructor, coercion, example, and deprecated-wrapper tests.

### Phase 4: Consolidate helpers and utilities
- [ ] Create `R/S7-stanfit-utils.R` for general stanfit-related helper functions
- [ ] Rename helper functions to neutral names where appropriate:
  - [ ] `read_stanfit()` instead of `read_ideal_adaptor_stanfit()`
  - [ ] `write_stanfit()` instead of `write_ideal_adaptor_stanfit()`
- [ ] Move caching, persistence, and refit logic into the new utility file

### Phase 5: Remove obsolete wrappers and fix broken code
- [ ] Remove the old S3/S4 wrapper layer for Stanfit-related class types
- [ ] Remove temporary compatibility shims that were only introduced to keep the package loadable during migration
- [ ] Delete or simplify any now-unnecessary helper functions such as old `is.ideal_adaptor_stanfit_input`-style shims if they are redundant with S7 validators or new class checks
- [ ] Fix any broken or convoluted logic introduced during earlier migration attempts, especially around transform handling, fit-object persistence, Stan input validation, and draws detection

### Testing plan
- [x] Add regression tests for the new Stan-input classes
- [x] Add regression tests for the new Stan-fit classes
- [ ] Add regression tests for `get_model_type()` and `get_representation_type()`
- [ ] Add regression tests for Stan accessors such as `get_stanfit()`, `set_stanfit()`, and `get_staninput()`
- [x] Verify the package still loads and the S7-focused tests pass after each phase

## Next Work

- [ ] Migrate NIW updating to operate on S7 models.
- [ ] Migrate likelihood, categorization, evaluation, and legacy info consumers that still require tibble models.
- [ ] Finish `tests/functions-to-make-or-load-models.R` S7 migration. `sample_data_from_model()` is gone; one `lift_MVG_ideal_observer_to_NIW_ideal_adaptor()` call remains and depends on the not-yet-migrated NIW workflow.
- [ ] Integrate fitted Stanfit objects with model-distribution objects.
- [ ] Clean up `get-info-from-NIW-IA-stanfit.R`, persistence helpers, compatibility shims, and fully deprecated files.
- [ ] Run the complete Stan-dependent test suite after the NIW workflow migration.

If you want, I can also turn this into a slightly more polished project-plan version with sections like “Goal”, “Scope”, and “Done criteria.”