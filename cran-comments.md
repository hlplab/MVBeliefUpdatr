## CRAN Submission Comments: MVBeliefUpdatr 0.1.0

Initial submission of MVBeliefUpdatr to CRAN.

## Test environments
* local macOS Sonoma 14.5 (Apple Silicon), R 4.5.2
* win-builder (x86_64-w64-mingw32), R Under development (unpatched)
* win-builder (x86_64-w64-mingw32), R 4.4.2 (release)

## R CMD check results

Status: 0 errors | 0 warnings | 0 notes

## Additional Repositories

* `cmdstanr` is in `Suggests` and its repository is declared in the
  `Additional_repositories` field in `DESCRIPTION`
  (https://stan-dev.r-universe.dev) in accordance with CRAN repository policy.

## Domain-Specific Terminology and Abbreviations

The following capitalized words or abbreviations in documentation are domain
terms in cognitive science, speech perception, and Bayesian statistics:
* `Adaptor`, `Adaptors`: Ideal adaptors in normative speech adaptation models
  (Kleinschmidt & Jaeger, 2015, Psychological Review).
* `NIW`: Normal-Inverse-Wishart conjugate distribution.
* `NIX`: Normal-Inverse-Chi-squared conjugate distribution.
* `MNIX`: Multivariate Normal-Inverse-Chi-squared distribution.

## Package Size and Pre-compiled Stan Models

* Built package archive size: ~4.8 MB (clearing the strict < 5.0 MB limit).
  Datasets in `data/` are compressed using xz (3.27 MB total).
* Installed library size is attributable to `libs/` (compiled C++ shared
  library). Like other Stan-based packages on CRAN (such as `rstanarm` and
  `brms`), MVBeliefUpdatr includes pre-compiled Stan models to provide instant
  MCMC sampling functionality without requiring users to configure local C++
  toolchains or recompile models during runtime.

## Companion Documentation Website

* Interactive 3D WebGL scenes and animated belief-updating figures that exceed
  CRAN package size guidelines are hosted online via pkgdown at:
  https://hlplab.github.io/MVBeliefUpdatr/
