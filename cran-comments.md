## CRAN Submission Comments: MVBeliefUpdatr 0.1.0

Initial submission of MVBeliefUpdatr to CRAN.

## Test environments
* local macOS Sonoma 14.5 (Apple M-series), R 4.5.2
* win-builder (x86_64-w64-mingw32), R Under development (unpatched)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Package Size and Compiled Shared Library Notes

* Built package archive size: 4.8 MB (clearing the strict < 5.0 MB limit).
  Datasets in `data/` are compressed using xz (3.27 MB total).
* Installed library size is primarily attributable to `libs/` (compiled C++
  shared library). Like other Stan-based packages on CRAN (such as `rstanarm`
  and `brms`), MVBeliefUpdatr includes pre-compiled Stan models to provide
  instant MCMC sampling functionality without requiring users to configure
  local C++ compilers or compile models during runtime.
