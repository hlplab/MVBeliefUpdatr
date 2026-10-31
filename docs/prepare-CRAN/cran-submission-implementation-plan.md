# Implementation Plan: CRAN Preparation and Submission Readiness

This implementation plan outlines the required steps, gates, and technical
criteria to prepare the `MVBeliefUpdatr` package for its initial submission
to CRAN.

## Overview & Current Status

- **Architecture:** Complete v0.1.0 S7 hierarchy (Phases 0 through 9).
- **Test Suite:** 45 test files, 100% passing (3,042 assertions, 0 failures).
- **R CMD Check (macOS):** 0 errors, 0 warnings, 0 notes.
- **Repository State:** Branch `S7-classes` up to date with remote.

---

## CRAN Submission Requirements & Action Items

```
+-----------------------------------------------------------------------------+
|                      CRAN Submission Preparation Path                       |
+-----------------------------------------------------------------------------+
| 1. Archive Size Optimization (< 5 MB tarball limit) [COMPLETED]             |
|    - Remove deprecated ChodroffWilson2018 (Done, saved 1.1 MB)              |
|    - Compress/resave data objects with xz (Done, data/ = 3.27 MB)           |
|    - Optimize vignette assets and figures (Done, archive = 4.8 MB)          |
|                                                                             |
| 2. Documentation & Help Page Examples                                       |
|    - Add compact, executable @examples (< 5 sec) to exported functions      |
|    - Protect slow MCMC fitting examples with \donttest{}                    |
|                                                                             |
| 3. Dependency Footprint Review                                              |
|    - Audit 30 non-default Imports                                           |
|    - Move optional/heavy libraries (plotly, gganimate) to Suggests          |
|                                                                             |
| 4. CRAN Policy Formatting in DESCRIPTION                                    |
|    - Expand acronyms on first use (MCMC, NIW, NIX, MNIX)                    |
|    - Enclose software/package names in single quotes ('Stan', 'S7', etc.)   |
|                                                                             |
| 5. Multi-Platform Validation (R-devel, Windows, Linux)                      |
|    - Run devtools::check_win_devel() and devtools::check_win_release()      |
|    - Verify GCC and Clang builds with rhub / GitHub Actions                 |
|                                                                             |
| 6. Authoring cran-comments.md                                               |
|    - Document test environments, check results, and explain compiled libs   |
+-----------------------------------------------------------------------------+
```

---

## Detailed Task Breakdown

### Phase 1: Archive and Package Size Optimization [COMPLETED]

CRAN strictly limits the submitted `.tar.gz` archive to under 5.0 MB.
The initial build measured 13.0 MB.

1. **Dataset Optimization (`data/`):** [DONE]
   - Removed `data/ChodroffWilson2018.RData` (-1.1 MB). `mixer6` now serves
     as the single canonical phonetic stop corpus.
   - Preserved active binding deprecation warning in `R/zzz.R` returning
     `mixer6`.
   - Executed `tools::resaveRdaFiles("data", compress = "xz")` to ensure
     maximum compression ratio across `h95.rda`, `pb52.rda`, `swehvd.rda`,
     and `mixer6.rda`.
   - Result: `data/` directory total size = 3.27 MB (cleared < 3.5 MB gate).

2. **Vignette Asset Optimization (`vignettes/`):** [DONE]
   - In `vignettes/visualizing-models-and-categories.Rmd`, set chunk
     evaluation to `eval = isTRUE(getOption("pkgdown"))`, allowing the
     online `pkgdown` website to render full 3D WebGL scenes and animations
     while compiling a lean 112 KB HTML vignette for the CRAN submission
     (down from 14.5 MB).
   - Added a prominent top banner directing readers to the rendered online
     documentation website.
   - Rescaled and quantized 8 factor graph diagrams in `vignettes/figures/`
     using ImageMagick (from 1.2 MB to 241 KB), shrinking the S7 vignette
     HTML from 1.98 MB to 693 KB.
   - Added `!vignettes/figures/*.png` to `.gitignore` and `^_pkgdown\.yml$`
     to `.Rbuildignore`.
   - Gate: Built `.tar.gz` package archive size = 4.8 MB (5,032,868 bytes),
     clearing the strict < 5.0 MB CRAN limit (down from 13.0 MB).

---

### Phase 2: Executable Examples (`@examples`)

CRAN incoming policy requires that all exported user-facing functions and
datasets include executable examples in their Rd help pages.

1. **Core Constructor Examples:**
   - Add minimal (< 1 second) examples to:
     - `cognitive_model()`
     - `prior()`
     - `fix_parameter()`
     - `as_model_list()`
     - `update_category_representation()`

2. **Evaluation & Plotting Examples:**
   - Add quick examples using built-in datasets (`pb52` or `h95`):
     - `categorize()`
     - `plot_categories()`
     - `plot_categorization_functions()`

3. **Computationally Intensive Stan Functions:**
   - For functions involving MCMC sampling or Stan compilation
     (`fit_ideal_adaptor()`, `add_criterion()`), wrap examples in
     `\donttest{...}` to comply with CRAN's 10-minute total check limit.
   - Gate: `checking examples ... OK` with 0 failures and execution < 5 min.

---

### Phase 3: Dependency Pruning & DESCRIPTION Polishing

1. **Dependency Audit (`Imports:` vs. `Suggests:`):**
   - Current `Imports` contains 30 non-default packages, triggering an INFO
     notice during CRAN incoming checks.
   - Move non-essential visualization engines to `Suggests:`:
     - `plotly` (used exclusively for interactive 3D mesh rendering)
     - `gganimate` (used exclusively for sequential update animations)
     - `tidybayes` (used for post-processing draws)
   - Protect their internal calls with `.assert_installed()` or
     `requireNamespace()`.

2. **DESCRIPTION Field Refinements:**
   - Expand all acronyms in `Description:` field:
     - MCMC -> Markov Chain Monte Carlo (MCMC)
     - NIW -> Normal-Inverse-Wishart (NIW)
     - NIX -> Normal-Inverse-Chi-squared (NIX)
     - MNIX -> Mixture of Normal-Inverse-Chi-squared (MNIX)
   - Wrap package and system names in single quotes:
     - `'Stan'`, `'S7'`, `'brms'`, `'ggplot2'`.
   - Gate: 0 incoming check notes on DESCRIPTION or dependencies.

---

### Phase 4: Cross-Platform & Compiler Verification

Stan packages compile C++ headers via `StanHeaders` and `RcppEigen`.
CRAN runs checks on diverse toolchains:

1. **Windows Checks:**
   - Submit to Win-Builder via `devtools::check_win_devel()` and
     `devtools::check_win_release()`.
   - Verify `configure.win` and compilation under Rtools44 UCRT.

2. **Linux & R-devel Checks:**
   - Check against Debian/Ubuntu R-devel via R-hub or containerized runner.
   - Ensure GCC-14 and Clang-19 do not raise fatal template warnings or
     strict aliasing errors.

3. **Gate:**
   - Zero errors, zero warnings across all three major platforms.

---

### Phase 5: Submission Preparation & Documentation

1. **`cran-comments.md` Creation:**
   - Draft standard CRAN submission comments:
     - Declaration of initial package submission.
     - Description of test environments (macOS Sonoma, Win-builder, R-devel).
     - Explanation of package size and compiled shared library (`libs/`)
       resulting from pre-compiled Stan models (matching `rstanarm`/`brms`).
     - Summary of 0 errors, 0 warnings, 0 notes.

2. **Final Verification Gate:**
   - Run `rcmdcheck::rcmdcheck(args = c("--as-cran"))`.
   - Confirm archive size <= 5.0 MB.
   - Proceed with `devtools::submit_cran()`.
