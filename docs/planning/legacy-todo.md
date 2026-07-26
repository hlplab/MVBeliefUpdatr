# Legacy To-Do (Migrated from to-do.txt)

This file preserves the historical task list, converted from plain text into Markdown.

## Very Important
- Correct MNIX formulation in stan code to match the ideal cue integration MNIX in the R code implemented as 
  part of the new S7 implementation
- Change handling of cue naming in mu, Sigma, Sigma_noise, m, S. Use attributes instead of naming each column.
- Decide whether univariate vs multivariate differences should be handled by umbrella sampling/density functions; simplify model print methods rather than storing mu, Sigma, Sigma_noise, m, S differently by dimensionality.
- Fix transform_cues when center = FALSE and return untransform function = TRUE.

## Important
- Change fit argument to model argument in get-info-from-NIW-IA-stanfit.R functions.
- Consider class definitions for MVBU representations and MVBU models.
- Check why get_IBBU_predicted_response reports many-to-many join warning.

## General
- Consider switching from assertthat to assertive.
- Consider making NIW_MCMC and NIW_belief one class, or NIW_MCMC as subclass of NIW_belief.
- Decide relationship between NIW_belief and NIW_beliefs classes.
- Add ability to add noise to data (cues_noisy) independently from update function.
- Update get_category_levels(), get_group_level(), and get_original_levels().

## NIW Belief Related Functions
- Expand update_beliefs to handle multiple exposure groups.

## IBBU Stanfit Related Functions
- Make add_ibbu_plot more efficient by filtering draws once.
- Add function to extract prior from stanfit in NIW_belief format.

## New IO Functions (Generalizable)
- Consider Dave's cross-validation functions; discuss with Xin.
- Build plot_category() for category plotting; keep expected-category extraction separate.
- Add options to plot mode (not only mean) and allow sampling from NIW.
- Add function to categorize one group's data using another group's beliefs.
- Add accessible plots for category means, log SDs, and correlations; update get_exposure_mean accordingly.
- Add accessible exposure covariance output (not uncentered squares); update getter accordingly.

## NIW Belief Related Plotting
- Add plot to check assumptions of NIW updating (distribution and covariance across talkers).
- Add plot of m and S for univariate/bivariate/trivariate cases (Murphy 2012, p. 127).
- Add standard plotting function for NIW_belief parameters.

## NIW IBBU Stanfit Related Plotting
- Add more X-like and more Y-like x-axis labels in test_categorization plots.

## Efficiency Improvements
- Speed up marginalization when prior is uniform; possibly derive analytic form instead of row-wise summation.
- Consider reworking add_ibbu_draws non-nested sample extraction, then nesting; could reduce code and improve transparency.

## Catch Up
- Implement tests for NIW functions.

## Done (Historical)
- Tested whether plot_expected_categories_2D works with plot.exposure = TRUE after preceding steps.
- Added information about groups and categories (number and order of levels).
- Added information about input to new class.
- Updated add_ibbu_draws() to return M and S style values while noting .stan variable rename still pending.
- Changed plot_expected_categorization to show category label in legend instead of category id.
- Explored class addition so new class is added to stanfit and inherits from it.
- Implemented class for mv_ibbu stanfits in class.R.
- Tested that add_ibbu_draws works on new object class.
- Tested that plotting functions work on new object class.