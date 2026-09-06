
## Next steps (Phase 4 / Extensions)

+ develop plotting methods that dynamically show model updates, both based on stanfits (from prior to posterior) or from update_model() outputs. 

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

### Extensions
+ implement mixture inference starting with predefined models that are handed to the stan code.

### Stan-related

+ temporarily change stan programs to echo inputs so that the correct structure of inputs can get verified. at that point revisit for broad range of inputs (1d-3d, nix/mnix/niw, w/ or w/o zero exposure) whether a simpler alternative to to_array can be found.

+ Ultimately, we need tests that generate test responses based on posterior (post-exposure) NIX/MNIX/NIW beliefs that were obtained by updating prior (pre-exposure) NIX/MNIX/NIW beliefs with the exposure data. That will let us check whether the 'forward updating' and the inferences of the prior beliefs are aligned---i.e., whether we can recover the prior beliefs provided there are sufficiently information exposure-test combinations in the data used to fit the stanfit model.
