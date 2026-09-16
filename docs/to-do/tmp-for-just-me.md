
## Next steps 

+ rename S7-core-*-classes.R to S7-class-*.R  and rename S7-core-methods to S7-methods

+ re: Resolved rhat method masking by removing redundant library(bayesplot) from the setup chunk and qualifying bayesplot::pp_check(). --> pp_check exists in this library. any reason not to use that? and wouldn't it be more elegant to define the S7 rhat method specifically for MVBU_stanfit objects, so that it only fires then?

+ check get_sufficient_category_statistics, which currently is a weird hybrid. even the method for data frames. standardize and join documentation with the same method for models. also only allow EITHER a model or data. Throw an error if both are provided. also think about in what format the suff. stats for different categories, cues should be returned (list, matrix, array, etc.)

+ for get_sufficient_category_statistics and ALL other functions that prepare, extract, or operate over the sufficient category statistics of stan models, check whether recent changes to what parameters get handed to stan models are correctly considered. ideally, we can standardize the name and the content of the variables we hand to stanfit models. it might be best to store the sufficient statistics in the meta data of the stanfit input (check structure of that) or the stanfit object. so that one can just grab those stats whole and return them to users. that seems to be the case for MVBU_StanfitPosterior already, but not sure whether it's the case for MVBU_Stanfit.

+ for MVBU_StanfitPosterior, check what information is stored and how that is structured: how does it compare to MVBU_Stanfit? make it parallel, but think about which structure is best.

+ in .compute_default_cue_limits, marginal moments are set via get_expectec_category_statistics. But i thought marginal moments marginalize over researcher uncertainty (posterior samples) *and* over the uncertainty of the NIX/NIW/MNIX, whereas get_expected_category_statistics, just marginalizes over the uncertainty of the NIX/NIW/MNIX model.

+ check why we need code like

```
 if (!"Sigma" %in% names(d_sum) && "S" %in% names(d_sum) && "nu" %in% names(d_sum)) {
    d_sum <- dplyr::mutate(
      d_sum,
      Sigma = get_expected_Sigma_from_S(.data$S, .data$nu)
    )
  }
```

in plotting functions. first of all, should expected sigma be cached in objects already? (if not, add it to the cache draws). second, doesn't the posterior density depend on S rather than sigma? it also seems like these plotting functions are still calculating likelihoods and posterior from scratch, rather than using pre-compiled likelihood and categorization functions.








# Updated error log

Codoc mismatches from Rd file 'update_category_representation.Rd': update_category_representation Code: function(x, x_N, x_mean, x_SS, ...) Docs: function(x, x_N, x_mean, x_SS) Argument names in code not in docs: ... update_category_representation Code: function(x, x_N, x_mean, x_SS, ...) Docs: function(x, x_N, x_mean, x_SS) Argument names in code not in docs: ... update_category_representation Code: function(x, x_N, x_mean, x_SS, ...) Docs: function(x, x_N, x_mean, x_SS) Argument names in code not in docs: ...

Codoc mismatches from Rd file 'update_template.Rd': update_template Code: function(x, observations, ...) Docs: function(x, observations, updating = c("batch", "incremental"), keep_history = FALSE, lapse_treatment = "no_lapses", noise_treatment = "no_noise", update_method = "label-certain") Argument names in code not in docs: ... Argument names in docs not in code: updating keep_history lapse_treatment noise_treatment update_method Mismatches in argument names: Position: 3 Code: ... Docs: updating update_template Code: function(x, observations, ...) Docs: function(x, observations, updating = c("batch", "incremental"), keep_history = FALSE, lapse_treatment = "no_lapses", noise_treatment = "no_noise", update_method = "label-certain") Argument names in code not in docs: ... Argument names in docs not in code: updating keep_history lapse_treatment noise_treatment update_method Mismatches in argument names: Position: 3 Code: ... Docs: updating update_template Code: function(x, observations, ...) Docs: function(x, observations, updating = c("batch", "incremental"), keep_history = FALSE, lapse_treatment = "no_lapses", noise_treatment = "no_noise", update_method = "label-certain") Argument names in code not in docs: ... Argument names in docs not in code: updating keep_history lapse_treatment noise_treatment update_method Mismatches in argument names: Position: 3 Code: ... Docs: updating

❯ checking Rd \usage sections ... WARNING Undocumented arguments in Rd file 'MVBU-TransformInformation-class.Rd' ‘transform.parameters’ ‘transform.function’ ‘untransform.function’

Undocumented arguments in Rd file 'MVBU_ModelList.Rd' ‘models’ ‘model_labels’ ‘metadata’

Undocumented arguments in Rd file 'MVBU_Stanfit.Rd' ‘data’ ‘staninput’ ‘stanvars’ ‘backend’ ‘save_pars’ ‘stan_args’ ‘stanfit’ ‘basis’ ‘transform_information’ ‘criteria’ ‘file’ ‘version’ ‘metadata’ ‘cache’

Undocumented arguments in Rd file 'MVBU_StanfitInput.Rd' ‘data’ ‘staninput’ ‘transform_information’ ‘metadata’

Undocumented arguments in Rd file 'MVBU_StanfitPosterior.Rd' ‘cues’ ‘categories’ ‘group’ ‘model_name’ ‘draws’ ‘metadata’

Undocumented arguments in Rd file 'MVBU_Staninput.Rd' ‘values’

Undocumented arguments in Rd file 'as_model_coercion.Rd' ‘category_template’

Undocumented arguments in Rd file 'categorize.Rd' ‘...’

Undocumented arguments in Rd file 'example-s7-objects.Rd' ‘n_cues’ ‘category’ ‘kappa’ ‘nu’ ‘type’ ‘...’ ‘categories’

Undocumented arguments in Rd file 'family-exemplar.Rd' ‘category_likelihood_function’ ‘category_posterior_functions’ ‘lapse_behavior’ ‘noise_behavior’ ‘cache’

Undocumented arguments in Rd file 'family-mnix.Rd' ‘category_likelihood_function’ ‘kappa’ ‘nu’ ‘category_posterior_functions’ ‘lapse_behavior’ ‘noise_behavior’ ‘cache’

Undocumented arguments in Rd file 'family-muvg.Rd' ‘category_likelihood_function’ ‘category_posterior_functions’ ‘lapse_behavior’ ‘noise_behavior’ ‘cache’

Undocumented arguments in Rd file 'family-mvg.Rd' ‘category_likelihood_function’ ‘category_posterior_functions’ ‘lapse_behavior’ ‘noise_behavior’ ‘cache’

Undocumented arguments in Rd file 'family-niw.Rd' ‘category_likelihood_function’ ‘category_posterior_functions’ ‘lapse_behavior’ ‘noise_behavior’ ‘cache’

Undocumented arguments in Rd file 'family-nix.Rd' ‘category_likelihood_function’ ‘category_posterior_functions’ ‘lapse_behavior’ ‘noise_behavior’ ‘cache’

Undocumented arguments in Rd file 'family-uvg.Rd' ‘category_likelihood_function’ ‘category_posterior_functions’ ‘lapse_behavior’ ‘noise_behavior’ ‘cache’

Undocumented arguments in Rd file 'fit_ideal_adaptor.Rd' ‘stanvars’ ‘save_pars’ ‘basis’ ‘chains’ ‘iter’ ‘warmup’ ‘init’ ‘control’ ‘verbose’

Undocumented arguments in Rd file 'get_category_likelihood_function.Rd' ‘...’

Undocumented arguments in Rd file 'get_category_posterior_function.Rd' ‘...’

Undocumented arguments in Rd file 'get_expected_category_statistic.Rd' ‘categories’ ‘groups’

Undocumented arguments in Rd file 'get_exposure_category_statistic.Rd' ‘categories’ ‘groups’ ‘statistic’ ‘untransform_cues’

Undocumented arguments in Rd file 'get_group_labels.Rd' ‘include_prior’

Undocumented arguments in Rd file 'get_metadata.Rd' ‘...’

Undocumented arguments in Rd file 'get_model_family.Rd' ‘...’

Undocumented arguments in Rd file 'get_model_type.Rd' ‘...’

Undocumented arguments in Rd file 'get_number_of_draws.Rd' ‘...’

Undocumented arguments in Rd file 'get_parameters.Rd' ‘original_pars’

Undocumented arguments in Rd file 'get_representation_type.Rd' ‘...’

Undocumented arguments in Rd file 'get_stanfit.Rd' ‘...’

Undocumented arguments in Rd file 'get_staninput.Rd' ‘...’

Undocumented arguments in Rd file 'get_sufficient_category_statistics.Rd' ‘untransform_cues’

Undocumented arguments in Rd file 'get_transform_information.Rd' ‘...’

Undocumented arguments in Rd file 'likelihood.Rd' ‘...’

Undocumented arguments in Rd file 'new_category_representation.Rd' ‘category_labels’ ‘cue_labels’ ‘category_likelihood_function’ ‘metadata’

Undocumented arguments in Rd file 'plot_categorization_functions.Rd' ‘lapse_treatment’

Undocumented arguments in Rd file 'posterior.Rd' ‘...’

Undocumented arguments in Rd file 'reconstruct_update_history.Rd' ‘full_posterior’ ‘step_size’

Undocumented arguments in Rd file 'set_staninput.Rd' ‘...’

Undocumented arguments in Rd file 'update_category_representation.Rd' ‘...’

Undocumented arguments in Rd file 'update_model_decision_bias_incrementally.Rd' ‘beta’

Functions with \usage entries need to have the appropriate \alias entries, and all their arguments documented. The \usage entries must correspond to syntactically valid R code. See chapter ‘Writing R documentation files’ in the ‘Writing R Extensions’ manual.

❯ checking R code for possible problems ... NOTE assert_MVBU_StanfitInput: no visible binding for global variable ‘MVBU_StanfitInput’ Undefined global functions or variables: MVBU_StanfitInput

# Current error log

Ran command: `Rscript -e "devtools::check(cran = FALSE, error_on = 'never')"`

I have launched `devtools::check(cran = FALSE, error_on = 'never')` in the background and will record all issues (errors, warnings, and notes) as soon as the check completes without taking any action.

Ran command: `Rscript -e "devtools::check(cran = FALSE, error_on = 'never', vignettes = FALSE)"`

I have recorded the first issue from the vignette build and launched `devtools::check(vignettes = FALSE)` to capture any additional issues across R CMD check (code, tests, dependencies, and documentation). I will present the complete diagnostic summary as soon as it finishes.

Viewed S7-stanfit-input.R:352-395
Viewed task-1900.log:1-200
Viewed task-1900.log:201-350
Viewed helper-s7-foundation.R:1-40

Here is the complete inventory of all issues identified by `R CMD check` and vignette compilation. No actions have been taken.

---

### 1. Errors (2)

1. **Vignette Build Failure (`visualizing-models-and-categories.Rmd`)**:
   * In chunk `stanfit_load_and_plot_categories` (lines 497–515):
     `plot_categories(fit, sample = TRUE, ndraws = 30)`
   * Fails at `R/S7-plot-methods.R:1839:9`: `as.matrix(d_raw$Sigma[[i]])` with error:
     `'data' must be of a vector type, was 'NULL'`.
   * Occurs because `d_raw` contains draws of scatter `S` rather than pre-computed `Sigma`, leaving `d_raw$Sigma` as `NULL`.

2. **Test Runner Failure under `R CMD check`**:
   * `helper-s7-foundation.R:26` aborts with:
     `stop("Could not locate the package R directory for test helper sourcing.")`.
   * Under `R CMD check`, tests execute against the installed package in a temporary check directory where source `R/` does not exist. The helper's directory detection assumption fails in this environment.

---

### 2. Warnings (8)

1. **Non-Portable File Names & Long Paths**:
   * File names with spaces:
     * `docs/planning/summary of phase 3 completions.md`
     * `docs/to-do/tmp-for just me.md`
   * Rd filenames exceeding the 100-byte tar limit:
     * `man/update_category_representation-MVBeliefUpdatr-MNIX_CategoryRepresentation-any-any-any-method.Rd`
     * `man/update_category_representation-MVBeliefUpdatr-NIX_CategoryRepresentation-any-any-any-method.Rd`

2. **Package Installation Warnings (Unused Arguments / Slot Mismatches in Legacy Helpers)**:
   * `new_muvg_category_representation`: unused arguments `mu`, `sigma2` in `.aggregate_category_representations`.
   * `new_mnix_category_representation`:
     * unused arguments `S`, `kappa`, `nu` in `.aggregate_category_representations`.
     * unused arguments `sigma2`, `kappa`, `nu` in `.create_ideal_adaptor_representation`.
     * unused arguments `kappa`, `nu`, `sigma2` in `.marginalize_representation_to_cues`.
     * unused arguments `m`, `sigma2`, `kappa`, `nu` in `.update_MNIX_category_representation_by_sufficient_statistics`.
     * partial argument matching of `m` to `metadata`.

3. **Non-ASCII Characters**:
   * `R/S7-plot-engine.R` contains non-ASCII characters (e.g. Unicode math/symbols or punctuation).

4. **Dependencies in R Code**:
   * Packages used via `::` or `requireNamespace` but undeclared in `DESCRIPTION` (Imports/Suggests):
     * `bridgesampling`
     * `metR`
   * Declared `Imports` in `DESCRIPTION` with no detected usage in the R code:
     * `RcppParallel`, `av`, `colorspace`, `cowplot`, `forcats`, `ggnewscale`, `ggridges`, `reshape2`, `stringr`, `transformr`, `viridis`.
   * Missing/unexported function call:
     * `plotly::add_isosurface` (not exported by `plotly`).
   * Internal package namespace access:
     * `MVBeliefUpdatr:::stanmodels` (should use `stanmodels` directly or export accessors).

5. **Missing Rd Cross-References**:
   * `stanvars` in `MVBU_Stanfit.Rd`
   * `[loo]{compare.loo}` in `loo_compare.MVBU_Stanfit.Rd`
   * `ideal_adaptor_stanfit` in `reshape_incremental_design_into_unique_exposure_test_combinations.Rd`
   * `update_model` in `update_category_representation-*.Rd` and `update_template-*.Rd`

6. **Missing Documentation Entries**:
   * Undocumented code object: `ideal_adaptor_stanfit`.

7. **Code / Documentation Default Mismatch**:
   * `MVBU-core-classes.Rd`: `MVBU_CognitiveModel` constructor signature default values differ between roxygen and code definition.

8. **Rd Objects in `\usage` without `\alias` / Undocumented Arguments**:
   * `update_template.Rd` and method files have undocumented arguments (`...`) or missing aliases for S7 method dispatch.

---

### 3. Notes (3)

2. **Top-Level Files**:
   * `LICENSE` is present in the repository root but not declared in `DESCRIPTION` (`License:` field).
3. **Undefined Global Functions or Variables in R Code**:
   * Undefined globals: `.eval_exemplar_density_grid`, `MVBU_StanfitInput`, `Phase`, `all_vars`, `bind_rows`, `check_exposure_test_data`, `collapse`, `column_to_rownames`, `everything`, `get_category_statistics_as_list_of_arrays`, `get_test_counts`, `i`, `map_int`, `mixer6`, `palette.colors`, `x_mean_exposure`.


# General cleanup as part of phase 7

# To do after Phase 9
+ think about whether evaluate_model should be renamed. it's job is similar to the add information criteria (IC) functions in brms, just that it's frequentist log-likelihood or accuracy vs. Bayesian ICs. perhaps this could be unified by changing the evaluate_models function into several separate add_* functions that add the IC to the model (similar to brms), allowing also a common compare models function that would extend that compare_models function from non-stanfit to stanfit model (allowing comparison on any IC)?

+ also check whether evaluate_model applied to stanfit or stanfit posterior objects , which currently uses method = "loglik", makes sense 


### Extensions
+ implement mixture inference starting with predefined models that are handed to the stan code.

### Stan-related

+ temporarily change stan programs to echo inputs so that the correct structure of inputs can get verified. at that point revisit for broad range of inputs (1d-3d, nix/mnix/niw, w/ or w/o zero exposure) whether a simpler alternative to to_array can be found.

+ Ultimately, we need tests that generate test responses based on posterior (post-exposure) NIX/MNIX/NIW beliefs that were obtained by updating prior (pre-exposure) NIX/MNIX/NIW beliefs with the exposure data. That will let us check whether the 'forward updating' and the inferences of the prior beliefs are aligned---i.e., whether we can recover the prior beliefs provided there are sufficiently information exposure-test combinations in the data used to fit the stanfit model.
