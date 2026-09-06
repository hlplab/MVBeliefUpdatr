# Legacy To-Do (Migrated from to-do.txt)

This file preserves the historical task list, converted from plain text into Markdown.

## Very Important
- Correct MNIX formulation in stan code to match the ideal cue integration MNIX in the R code implemented as 
  part of the new S7 implementation
- Decide whether univariate vs multivariate differences should be handled by umbrella sampling/density functions; simplify model print methods rather than storing mu, Sigma, Sigma_noise, m, S differently by dimensionality.
- Fix transform_cues when center = FALSE and return untransform function = TRUE.

## NIW Belief Related Functions
- Expand update_beliefs to handle multiple exposure groups.

## New IO Functions (Generalizable)
- Consider Dave's cross-validation functions; discuss with Xin.
- Add options to plot mode (not only mean) and allow sampling from NIW.
- Add function to categorize one group's data using another group's beliefs.

## NIW Belief Related Plotting
- Add plot to check assumptions of NIW updating (distribution and covariance across talkers).
- Add plot of m and S for univariate/bivariate/trivariate cases (Murphy 2012, p. 127).

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