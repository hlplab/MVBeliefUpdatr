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

+ DONE: comprehensive plotting vignette `vignettes/visualizing-models-and-categories.Rmd` created with clickable table of contents (TOC), interactive and static 1D/2D/3D visualizations, marginal projections, and Stanfit diagnostic plots; renders cleanly with 0 errors to `vignettes/visualizing-models-and-categories.html`.
+ DONE: `plot_categories` aesthetic default aligned to line-only (`aes = "contour"`) for both single representations and multi-category templates.
+ DONE: exemplar plotting rug alpha set to `1 / length(reps)`, rug tick height halved to `0.015`, default exemplar points plotted set to 0 (`n_exemplars = 0L`), and subtitle dynamically indicating `(sampling XXX exemplars)`.
+ DONE: `sample_observations` (and alias `sample_observation`) turned into S7 generic and methods in `R/S7-core-methods.R` with uniform sampling across template categories, prior-weighted sampling across model categories, `with_replacement = TRUE/FALSE` support on exemplar objects with informative error messages, and full roxygen documentation.
+ DONE: category representation template plot subtitle simplified to state representation type without category count; categorization function subtitle concisely lists `all = {1 / K}` when all category priors or lapse biases are identical.
+ DONE: vectorized exemplar likelihood evaluation via Cholesky factor precomputation and `tcrossprod`, accelerating KDE evaluation across all plotting functions.
+ DONE: comprehensive plotting vignette `vignettes/visualizing-models-and-categories.Rmd` updated with single automatic top TOC, embedded interactive 3D WebGL Plotly scene with custom matte lighting/camera controls, mathematical exposition of `levels` and $\chi^2_2$ central density regions, and comprehensive exploration of arguments (`levels`, `limits`, `resolution`, `categories`, `decision_rule`, `combine_plots`, `parameters`, `ndraws`, `aes`). Renders cleanly with 0 errors to `vignettes/visualizing-models-and-categories.html`.
+ DONE: unit tests added for `sample_observations` and plotting defaults; verified full test suite passes cleanly with 0 failures across 1,392 tests.

+ DONE: linked Xie et al. (2023) reference to journal article URL in vignette.
+ DONE: simplified category plot titles to remove redundant "representation" ("1D multivariate Gaussian categories").
+ DONE: clarified argument defaults across documentation and vignette (single unambiguous default per parameter).
+ DONE: added 2D interactive `plot_categories()` with `aes = "fill-discrete"`, `"fill-gradient"`, and `"contour"`, removing mode marker and boosting surface opacity.
+ DONE: added 3D interactive `plot_categories()` with `aes = "fill-discrete"`, `"fill-gradient"`, and `"contour"`, setting default opacity to 0.65.
+ DONE: standardized `levels` default across all `plot_categories` (1:4 sigmas) except 3D ellipsoids (2 sigma); default slices at pooled mean $\pm \sigma$ intervals.
+ DONE: added text labels on 3D slice contour lines displaying enclosed probability mass percentages, with line widths scaling inversely with enclosed mass.
+ DONE: added `aes = "scatter"` for 3D interactive exemplar category plots, deduplicated category legends, and documented category mean diamond markers.
+ DONE: fixed vignette function calls to avoid passing explicit default arguments, stating defaults in the text descriptions.
+ DONE: added dedicated parallelization and performance optimization section in vignette.
+ DONE: unified 1D and 2D category legends using `key_glyph = draw_key_rect` and matching `override.aes` so `patchwork::plot_layout(guides = "collect")` produces a single collected legend.
+ DONE: resolved empty 3D sliced category plots by defaulting `aes` to `"contour"`.
+ DONE: updated `plot_categorization_function()` to apply `decision_rule` ("criterion" vs "proportional"), `noise_treatment`, and `lapse_treatment` through cached S7 posterior closures.
+ DONE: verified all markdown bulleted lists across `vignettes/visualizing-models-and-categories.Rmd` are preceded by a blank line for clean HTML rendering.
+ DONE: removed alias `plot_categorization` across code, documentation, and tests in favor of `plot_categorization_function()`.
+ DONE: updated default `levels` from 1:4 sigmas to 1:3 sigmas (`2 * stats::pnorm(1:3) - 1`) across all category and sliced plots.
+ DONE: set default aesthetic for 2D interactive `plot_categorization_function()` to `"fill-discrete"`, rendering a solid 3D response surface with opacity 0.85.
+ DONE: corrected 3D exemplar sliced contour plot in vignette by setting slices along `vowel_duration` (50, 100, 150 ms) to resolve empty panels.
+ DONE: defined `get_noise_treatment` and `get_lapse_treatment` S7 generics and methods on `MVBU_CognitiveModel`, `MVBU_Object`, and `S7::class_any`, replacing direct slot accesses across `R/S7-core-methods.R`, `R/S7-plot-engine.R`, `R/S7-plot-methods.R`, `R/S7-family-coercion.R`, `R/S7-make-objects.R`, and `R/S7-update-model.R`.
+ DONE: removed `sample_observation` alias, retaining only `sample_observations`.
+ DONE: retired `wide` argument from `get_draws()` generic and method (issuing `lifecycle::deprecate_warn("0.0.9", "get_draws(wide = )")` if supplied), and removed redundant `wide = FALSE` calls across plotting and stanfit methods.
+ DONE: reorganized and standardized test suite naming and execution order: all 26 non-deprecated tests are sequentially numbered `test-01-` through `test-26-` (with `S7-` prefix for S7 tests), and all 9 deprecated test files are numbered `test-80-` through `test-88-` (`test-80-deprecated-...` through `test-88-deprecated-...`) to run strictly after all non-deprecated tests; verified all 1,414 tests pass cleanly with 0 failures.

+ DONE: Phase 3 is COMPLETED and CLOSED. S7 class structure, bare/from_data constructors, print/summary/as_tibble, Stanfit classes, model aggregation (`aggregate_models()`), and vignettes (`s7-class-structure-and-workflows.Rmd`, `backward-compatibility-working-with-old-code.Rmd`) are fully implemented and verified with 1,494 passing tests.

# To do in Phase 3
Phase 3 is closed.


