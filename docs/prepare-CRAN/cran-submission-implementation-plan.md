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
| 2. Documentation & Help Page Examples [COMPLETED]                           |
|    - Compact executable @examples on all 199 help pages (< 1 sec default)   |
|    - Protect MCMC sampling with \donttest{} (total examples = 53.6s)        |
|                                                                             |
| 3. Dependency Footprint Review [COMPLETED]                                  |
|    - Audit non-default Imports (pruned plotly, gganimate, tidybayes)        |
|    - Move optional/heavy libraries to Suggests, patchwork to Imports        |
|    - Protect optional package calls with .assert_installed()                |
|                                                                             |
| 4. CRAN Policy Formatting in DESCRIPTION [COMPLETED]                         |
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

### Phase 2: Executable Examples (`@examples`) [COMPLETED]

CRAN incoming policy requires that all exported user-facing functions and
datasets include executable examples in their Rd help pages.

1. **Core Constructor Examples:**
   - Added minimal (< 1 second) examples to all constructors, generics,
     methods, and dataset documentation across all 199 Rd help pages.
   - Verified 0 missing examples across exported functions.

2. **Evaluation & Plotting Examples:**
   - Added quick examples using built-in datasets (`pb52`, `h95`):
     - `categorize()`
     - `plot_categories()`
     - `plot_categorization_functions()`
     - `plot_parameters()` (with native plotmath expressions for μ, τ, ρ,
       κ, ν, λ)

3. **Computationally Intensive Stan Functions & Optimization:**
   - Optimized `example_ideal_adaptor_stanfit()` default parameters
     (`chains = 1L, iter = 50L, warmup = 25L`, 95-row dataset) reducing
     sampling overhead from 85s down to 2.69s.
   - Wrapped long-running MCMC sampling examples in `\donttest{...}`.
   - Re-exported `loo::loo_compare` generic to resolve example dispatch.
   - Gate: `checking examples ... OK (21.1s)` and
     `checking examples with --run-donttest ... OK (32.5s)`. Total runtime
     53.6s with 0 errors across all 199 help pages.

---

### Phase 3: Dependency Pruning & DESCRIPTION Polishing [COMPLETED]

1. **Dependency Audit (`Imports:` vs. `Suggests:`):** [DONE]
   - Moved heavy/optional packages (`plotly`, `gganimate`, `tidybayes`) to
     `Suggests:`.
   - Moved `patchwork` from `Suggests:` to `Imports:` (used directly for
     multi-panel `plot_parameters` layouts) and added `@importFrom`.
   - Added `.assert_installed()` helper to `R/internal-asserts.R` to guard
     all calls to optional suggested packages with informative installation
     guidance.

2. **DESCRIPTION Field Refinements:** [DONE]
   - Expanded all acronyms in `Description:` field:
     - MCMC -> Markov Chain Monte Carlo (MCMC)
     - NIW -> Normal-Inverse-Wishart (NIW)
     - NIX -> Normal-Inverse-Chi-squared (NIX)
     - MNIX -> Multivariate Normal-Inverse-Chi-squared (MNIX)
   - Wrapped package and system names in single quotes:
     - `'Stan'`, `'S7'`, `'loo'`, `'ggplot2'`, `'patchwork'`.
   - Outcome: 0 errors, 0 warnings, 0 notes.

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
