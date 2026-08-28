# Principles for migration of old functionality:

Here's an update list of steps for the conversion of legacy functionality into the new S7 approach:

1) expand the new methods to cover the functionality of the legacy functions. E.g. when get_cue_labels_from_model allows specification of indices, so that only the specified cue labels are returned, integrate that idea into the new get_cue_labels.

2) make sure that similar logic is implemented consistently. e.g., if get_category_labels_from_model doesn't happen to have the indeces functionality, nevertheless add it to the new get_category_label

3) document the new functionality by copying (standardized) descriptions from the legacy roxygen documentation to the documentation of the (revised) new methods.

4) replace all calls  to the old legacy functions anywhere in the legacy code with calls to the new methods

5) mark the legacy functions that have been integrated into the new method as deprecated.

6) port tests of the legacy function to the new method (organizing test files in the same way that they were organized in the legacy code).

7) run ported tests for new methods and run tests for higher-level legacy code that now uses the new methods instead of the old legacy lower-level functions.

8) repeat step 1-7 (each time checking with me), starting from the lowest-level functions to increasingly more high-level functions that call the low-level functions.

9) whenever all code in an R file is marked as deprecated, change the filename by prefixing "deprecated-post-s7-".


# To do in Phase 3

## Next steps


- Migrate NIW update workflows to S7 objects, then migrate the higher-level likelihood, categorization, evaluation, and legacy information consumers that still assume tibble models. 
- then legacy tibble-based workflows can be removed though deprecated functions should be retained as thin wrappers. keep evaluate_model but completely convert it to S7-based workflow (it was previously flagged as depcrecated but I decided to keep it); if necessary adjust its argument naming and workflow, check its logic, remove/reduce dependencies on tidyverse or other packages, and update the documentation accordingly. 

+ Known, expected fallout (out of scope for this task): test-10-evaluate-model, test-12-get-info-MVG, test-13-get-info-NIW, and test-20-update-NIW now fail because they feed example_MVG_ideal_observer()/example_NIW_ideal_adaptor() output into legacy tibble-only consumers (evaluate_model, get_likelihood_from_MVG, get_categorization_from_MVG_ideal_observer, get_NIW_posterior_predictive, update_NIW_ideal_adaptor_incrementally, is.MVG_ideal_observer, and $ tibble-column access). Fixing these requires migrating those legacy consumers to accept S7 model objects — a separate, larger follow-up task. 

+ DONE: `likelihood()` generic added (S7-generics.R) with methods for representation, template, and cognitive model; `get_category_likelihood_function` is now exported. For NIX/MNIX/NIW the likelihood *is* the posterior predictive, so `get_NIW_posterior_predictive()`, `get_posterior_predictive_from_NIW_belief(s)()`, and `get_MVG_likelihood()` can now be reduced to thin wrappers over `likelihood()`. Still to do: plotting code in plot-expected-categories.R continues to call the family-specific helpers directly.

+ DONE: `example_*_from_data()` renamed to `example_*()` (the `_from_data` suffix was unintended). `new_*_from_data()` constructors keep the suffix.

+ DONE: `get_posterior_from_MVG_ideal_observer`, `get_categorization_from_MVG_ideal_observer`, `get_likelihood_from_MVG`, and `get_categorization_from_NIW_ideal_adaptor` are now thin wrappers over the S7 methods (via `.legacy_long_posterior()` / `.legacy_apply_decision_rule()` in get-info-from-model.R). `is.MVG_ideal_observer()` recognizes S7 objects.

+ Known remaining failures (pre-existing, not caused by the above): MNIX representation validator (test-01, test-07), NIX stanfit validator (test-05), rstan/TBB toolchain dlopen (test-04-stanfit-input-compatibility), and `tests/functions-to-make-or-load-models.R` still calling `mutate()` on S7 models (test-14-get-info-stanfit, test-19-plot-stanfit).


+ deprecate function to aggregate models with a indicator for an alternative workflow: use map to create lists of models, aggregate list of model into a model (write a new function for that).


+ now that the issue with corrupt data is fixed, I'm wondering whether we can simplify something that was introduced much earlier in order to avoid a similar issue. We are current re-evoking a new constructor each time we're setting a stanfit and perhaps in some other places:

```
S7::method(set_stanfit, list(MVBU_Stanfit, S7::class_any)) <- function(x, stanfit) {
  if (!is.null(stanfit)) {
    .assert_stanfit(stanfit)
    .assert_that(
      stanfit@model_name %in% names(MVBeliefUpdatr:::stanmodels),
      msg = paste0(
        "stanfit object was not created by one of the accepted stancodes:\n\t",
        paste(names(MVBeliefUpdatr:::stanmodels), collapse = "\n\t"),
        "\n(you can get the name of your model from your_stanfit@model_name)."
      )
    )
  }

  constructor <- get_ideal_adaptor_stanfit_constructor(x@staninput)

  constructor(
    data = x@data,
    staninput = x@staninput,
    stanvars = x@stanvars,
    backend = x@backend,
    save_pars = x@save_pars,
    stan_args = x@stan_args,
    stanfit = stanfit,
    basis = x@basis,
    transform_information = x@transform_information,
    criteria = x@criteria,
    file = x@file,
    version = x@version,
    labels = x@labels
  )
}
```

could it be that simpler code like `x@stanfit <- stanfit; x` instead of:

```
constructor <- get_ideal_adaptor_stanfit_constructor(x@staninput)

  constructor(
    data = x@data,
    staninput = x@staninput,
    stanvars = x@stanvars,
    backend = x@backend,
    save_pars = x@save_pars,
    stan_args = x@stan_args,
    stanfit = stanfit,
    basis = x@basis,
    transform_information = x@transform_information,
    criteria = x@criteria,
    file = x@file,
    version = x@version,
    labels = x@labels
  )
  ```

  would suffice?

+ clean up get-info-from-NIW-IA-stanfit.R, making those functions methods and removing functions no longer needed.
+ change handling/storage of label information in stanfit objects to follow the same structure used in the core model/representation/etc. objects
+ make print methods for the new S7 classes. summary should yield same as print, except for stanfit.
+ check which utils are needed. 
++ if almost all checks of scalars are actually for non-NA scalars change the .is_X to include requirement for non-NA, remove .is_non_NA_x, and also adjust .assert functions accordingly.
++ for overridden functions check whether they are still necessary.

+ integrate the fitted stanfit objects with the model distribution object as part of fit_ideal_adaptor

+ see whether the helper functions for testing can be simplified/reduced. e.g., make_vowel_test_data might be replaced by other test data by adjusting the tests, while yielding the same coverage?

## Code clarity 
+ check how deprecated functions are marked in terms of their roxygen documentation. is it consistent? ideally, they should not be listed in the table of content of help files, but should have help files. the structure of those help files and the way that deprecation warnings are given should be consistent across deprecated functions.

also check whether we can switch to one warning per session and ensure that warnings are only given "always" during testing?

## Testing
+ Minimal example of fitted models that can be evaluated. And where should they be stored?

# To do after Phase 3
+ check whether we can get rid of the functions in override.R 
+ check whether there are repeated code chunks that should be consolidated into internal helper functions.
+ Check MNIX cue weighting: there are at least two possible implementations:
  + adapt after integration (single nix sitting on top of integrated cue dimension)
  + adapt before integraton (one nix for each cue; then integration)
  + make sure that whatever is implemented for Stan is also what is implemented in R
  + make sure that these choices are documented verbosely in roxygen
+ MNIX fitting is current commented out. AI comments for MNIX failure:

" The root cause is visible now: the MNIX Stan program expects a covariance-style summary array, but the current builder is providing a sum-of-squares vector array instead. I’m aligning that data structure with the model’s declared interface before I verify again."

## Extensions
+ write forward-updating functions for NIX, MNIX, NIW, revising existing NIW forward-updating.
+ add specification of category prior and lapse bias via fixed_parameters and extend Stan program to include category prior

## Validity checks
+  for categorize(), demonstrate the output of the three different decision rules to me with an example


### Stan-related

+ temporarily change stan programs to echo inputs so that the correct structure of inputs can get verified. at that point revisit for broad range of inputs (1d-3d, nix/mnix/niw, w/ or w/o zero exposure) whether a simpler alternative to to_array can be found.

+ Ultimately, we need tests that generate test responses based on posterior (post-exposure) NIX/MNIX/NIW beliefs that were obtained by updating prior (pre-exposure) NIX/MNIX/NIW beliefs with the exposure data. That will let us check whether the 'forward updating' and the inferences of the prior beliefs are aligned---i.e., whether we can recover the prior beliefs provided there are sufficiently information exposure-test combinations in the data used to fit the stanfit model.

## Efficiency conisderations
+ Consider making separate versions of stan models for 0, 1, or more observations. functions could be shared between them to ease maintenance
