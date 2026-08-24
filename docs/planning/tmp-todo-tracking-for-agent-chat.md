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
+ create likelihood functions and deprecate old ones. deprecate fully depcrecated R files.


+ currently ongoing: 

I don't like that we have to explicit use the constructor in set_stanfit:
if (S7::S7_inherits(x, NIX_IdealAdaptorStanfit)) {
    constructor <- NIX_IdealAdaptorStanfit
  } else if (S7::S7_inherits(x, MNIX_IdealAdaptorStanfit)) {
    constructor <- MNIX_IdealAdaptorStanfit
  } else if (S7::S7_inherits(x, NIW_IdealAdaptorStanfit)) {
    constructor <- NIW_IdealAdaptorStanfit
  } else if (S7::S7_inherits(x, IdealAdaptorStanfit)) {
    constructor <- IdealAdaptorStanfit
  } else {
    constructor <- MVBU_Stanfit
  }

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

+ clean up get-info-from-NIW-IA-stanfit.R, making those functions methods and removing functions no longer needed.
+ change handling/storage of label information in stanfit objects to follow the same structure used in the core model/representation/etc. objects
+ make new constructors to replace make_*_from data, turn those legacy functions into wrappers and deprecate them.
+ change current wrappers as_{.class}* to replace the lift-* methods
+ make print methods for the new S7 classes. summary should yield same as print, except for stanfit.
+ check which utils are needed. 
++ if almost all checks of scalars are actually for non-NA scalars change the .is_X to include requirement for non-NA, remove .is_non_NA_x, and also adjust .assert functions accordingly.
++ for overridden functions check whether they are still necessary.

+ integrate the fitted stanfit objects with the model distribution object as part of fit_ideal_adaptor

## Code clarity 

## Testing
+ Minimal example of fitted models that can be evaluated. And where should they be stored?

# To do after Phase 3
+ check whether we can get rid of the functions in override.R 
+ check whether there are repeated code chunks that should be consolidated into internal helper functions.

## Extensions
+ write forward-updating functions for NIX, MNIX, NIW, revising existing NIW forward-updating.
+ add specification of category prior and lapse bias via fixed_parameters and extend Stan program to include category prior

## Validity checks
+  for categorize(), demonstrate the output of the three different decision rules to me with an example


### Stan-related
+ Check MNIX cue weighting
+ MNIX fitting is current commented out. AI comments for MNIX failure:

" The root cause is visible now: the MNIX Stan program expects a covariance-style summary array, but the current builder is providing a sum-of-squares vector array instead. I’m aligning that data structure with the model’s declared interface before I verify again."



+ temporarily change stan programs to echo inputs so that the correct structure of inputs can get verified. at that point revisit for broad range of inputs (1d-3d, nix/mnix/niw, w/ or w/o zero exposure) whether a simpler alternative to to_array can be found.

+ Ultimately, we need tests that generate test responses based on posterior (post-exposure) NIX/MNIX/NIW beliefs that were obtained by updating prior (pre-exposure) NIX/MNIX/NIW beliefs with the exposure data. That will let us check whether the 'forward updating' and the inferences of the prior beliefs are aligned---i.e., whether we can recover the prior beliefs provided there are sufficiently information exposure-test combinations in the data used to fit the stanfit model.

## Efficiency conisderations
+ Consider making separate versions of stan models for 0, 1, or more observations. functions could be shared between them to ease maintenance
