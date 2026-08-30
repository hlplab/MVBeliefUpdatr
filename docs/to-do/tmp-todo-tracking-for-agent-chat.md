# Done

+ DONE: change all the S3-based classes in S7-stanfit-diagnostics.R to S7 generics (`log_posterior`, `nuts_params`, `rhat`, `neff_ratio`, `control_params`).
+ DONE: rename get_params to get_parameter_names to avoid confusion with get_parameters.
+ DONE: change all deprecation warnings (which currently either refer to 1.0.0 or 0.0.4) to refer to version 0.0.9

+ DONE: `likelihood()` generic added (S7-generics.R) with methods for representation, template, and cognitive model; `get_category_likelihood_function` is now exported. For NIX/MNIX/NIW the likelihood *is* the posterior predictive, so `get_NIW_posterior_predictive()`, `get_posterior_predictive_from_NIW_belief(s)()`, and `get_MVG_likelihood()` can now be reduced to thin wrappers over `likelihood()`. Still to do: plotting code in plot-expected-categories.R continues to call the family-specific helpers directly.

+ DONE: `example_*_from_data()` renamed to `example_*()` (the `_from_data` suffix was unintended). `new_*_from_data()` constructors keep the suffix.

+ DONE: `get_posterior_from_MVG_ideal_observer`, `get_categorization_from_MVG_ideal_observer`, `get_likelihood_from_MVG`, and `get_categorization_from_NIW_ideal_adaptor` are now thin wrappers over the S7 methods (via `.legacy_long_posterior()` / `.legacy_apply_decision_rule()` in get-info-from-model.R). `is.MVG_ideal_observer()` recognizes S7 objects.

+ DONE: it seems documentation from many of the original S3 methods was lost when they were transferred to S7. make sure add that documentation to the S7 methods e.g., in S7-stanfit-methods.R and S7-stanfit-diagnostics.R.

+ DONE:  in get_test_data, perhaps .from_staninput = FALSE should be renamed to .recover_from_staninput. if so, only if the data property is missing and .recover_from_staninput = TRUE should the method try to recover the data from the staninput (adjust documentation)

+ DONE: several functions, in e.g., basics and elsewhere, still seem to maintain checks and conditional code path using legacy is.* functions. but we don't need to maintain this type of backward compatability.

+ DONE: `get_D()` deprecated in favor of dimensionality extraction on representations/models; moment conversions (`get_expected_mu_from_m`, `get_m_from_expected_mu`, `get_expected_Sigma_from_S`, `get_S_from_expected_Sigma`) and posterior predictives (`get_NIW_posterior_predictive`, `get_NIX_posterior_predictive`) moved to `NIX-NIW-basics.R`.
+ DONE: `evaluate_model` turned into an S7 generic with methods for `MVBU_CognitiveModel` in `S7-core-methods.R` and `MVBU_Stanfit` in `S7-stanfit-methods.R`.
+ DONE: all legacy `get-info-from-*.R` files converted to deprecated wrappers in `deprecated-get-info-from-*.R` (including `deprecated-get-info-from-model.R` which wraps S7 generics/accessors and warns with `lifecycle`), and internal `.get_D` and `format_input_for_likelihood_calculation` removed.
+ DONE: verified full test suite passes cleanly with 0 failures across 1,354 tests.

+ DONE: `.infer_noise_treatment` internal helper created in `R/basics.R` and wired across all callers.
+ DONE: verified NIX example models (1-cue shifted prior) are fitted and tested in `test-05-fit-ideal-adaptor-minimal-examples.R`.
+ DONE: `get_sufficient_category_statistics` cleaned up (univariate legacy comments removed, `verbose = FALSE` added to formal args), `make_named_vector` and `make_named_square_matrix` deprecated in `deprecated-nest-model.R`, `make_vector_column` retained for data wrangling.

# To do in Phase 3

## Next steps
+ fit example stanfits for NIX model

+ there is quite a bit of overlap between the different plotting methods for different object types. if it does not make the code to opaque, try to streamline this through shared helper functions. consider whether some compute-intensive tasks during plotting could take advantage of parallelization (e.g. when working with samples from stanfit models; or when calculating grids for densities). keep in mind that we eventually also want to plot 3d plots. 


+ make 3d plot functions, too

WAIT, do not go beyond this point: 

+ Known remaining failures (pre-existing, not caused by the above): MNIX representation validator (test-01, test-07), NIX stanfit validator (test-05), rstan/TBB toolchain dlopen (test-04-stanfit-input-compatibility), and `tests/functions-to-make-or-load-models.R` still calling `mutate()` on S7 models (test-14-get-info-stanfit, test-19-plot-stanfit).

+ for the aggregate_models function, provide more documentation how the aggregation takes place, i.e., how exemplars, means, sigmas, m, S, etc. aggregated? for each model parameter, list the function used to create the aggregate. this might become clearner if we defined an aggregate S7 method that can aggregate category representations, templates, and models (Described in detail on a shared help page)

# General cleanup at end of Phase 3
+ check which utils are needed. 
++ if almost all checks of scalars are actually for non-NA scalars change the .is_X to include requirement for non-NA, remove .is_non_NA_x, and also adjust .assert functions accordingly.
++ for overridden functions check whether they are still necessary.

+ see whether the helper functions for testing can be simplified/reduced. e.g., make_vowel_test_data might be replaced by other test data by adjusting the tests, while yielding the same coverage?

+ check how deprecated functions are marked in terms of their roxygen documentation. is it consistent? ideally, they should not be listed in the table of content of help files, but should have help files. the structure of those help files and the way that deprecation warnings are given should be consistent across deprecated functions.

also check whether we can switch to one warning per session and ensure that warnings are only given "always" during testing?

## Testing

# To do after Phase 3
+ check whether we can get rid of the functions in override.R 
+ check whether there are repeated code chunks that should be consolidated into internal helper functions.

## Validity checks
+  for categorize(), demonstrate the output of the three different decision rules to me with an example


### Stan-related

+ temporarily change stan programs to echo inputs so that the correct structure of inputs can get verified. at that point revisit for broad range of inputs (1d-3d, nix/mnix/niw, w/ or w/o zero exposure) whether a simpler alternative to to_array can be found.

+ Ultimately, we need tests that generate test responses based on posterior (post-exposure) NIX/MNIX/NIW beliefs that were obtained by updating prior (pre-exposure) NIX/MNIX/NIW beliefs with the exposure data. That will let us check whether the 'forward updating' and the inferences of the prior beliefs are aligned---i.e., whether we can recover the prior beliefs provided there are sufficiently information exposure-test combinations in the data used to fit the stanfit model.


