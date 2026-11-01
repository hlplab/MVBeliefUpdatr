# Refactor Plotting Architecture, Add Parallelization & 3D Visualization

## Overview

This document specifies the architecture and implementation plan to refactor
the `MVBeliefUpdatr` plotting system. The goals are:
1. Eliminate redundant evaluation and layering code across object types
   (representations, templates, cognitive models, Stanfit objects).
2. Standardize on a single parallelization mechanism (`parallel` package) for
   compute-intensive plotting steps (MCMC posterior predictive draws,
   high-resolution 2D grids, and 3D voxel spaces).
3. Support 3-cue spaces ($D = 3$) seamlessly:
   - Default for `plot_categories()` / `plot_categorization()`: multi-panel 2D
     slices across quantiles of the 3rd cue.
   - Interactive 3D visualizations via `plot_categories_interactive()` and
     `plot_categorization_interactive()` for both 2-cue 3D surfaces ($z =$
     density or posterior) and 3-cue 3D ellipsoids / isosurfaces.
4. Maintain strict line length compliance ($\le 80$ characters per line) and
   full backwards compatibility.

---

## Architectural Schema

```
                  ┌────────────────────────────────────────┐
                  │          S7 Plotting Generics          │
                  │  plot_categories() plot_categorize()   │
                  │  plot_parameters() plot_correlate()    │
                  │  plot_*_interactive()                  │
                  └───────────────────┬────────────────────┘
                                      │
           ┌──────────────────────────┴──────────────────────────┐
           ▼                                                     ▼
┌───────────────────────┐                             ┌───────────────────────┐
│ Core S7 Objects       │                             │ Stanfit S7 Objects    │
│ Reps, Templates,      │                             │ MVBU_Stanfit,         │
│ Cognitive Models      │                             │ MCMC Posterior Draws  │
└──────────┬────────────┘                             └──────────┬────────────┘
           │                                                     │
           └──────────────────────────┬──────────────────────────┘
                                      │
                                      ▼
                  ┌────────────────────────────────────────┐
                  │    Shared Grid & Evaluation Engine     │
                  │ • Grid generation (1D, 2D, 3D)         │
                  │ • Density & posterior calculations     │
                  │ • Parallel execution (parallel pkg)    │
                  └───────────────────┬────────────────────┘
                                      │
          ┌───────────────────────────┼───────────────────────────┐
          ▼                           ▼                           ▼
┌───────────────────┐       ┌───────────────────┐       ┌───────────────────┐
│ 1D Visualizations │       │ 2D Visualizations │       │ 3D Visualizations │
│ • Density curves  │       │ • Contours/fills  │       │ • 3D 2-cue surface│
│ • Decision bounds │       │ • Gradient probs  │       │ • Multi-2D slices │
│ • Parameter HDIs  │       │ • Pairwise matrix │       │ • 3D 3-cue spheres│
└───────────────────┘       └───────────────────┘       └───────────────────┘
```

---

## Modular File Structure

The plotting system will be organized into four focused files:

### 1. `R/S7-plot-engine.R` (Shared Internal Evaluation Engine)
Contains internal computation and layer construction functions:
- `.mvbu_make_grid(limits, n_points, dims)`: Generates regular 1D, 2D, and
  3D coordinate grids.
- `.mvbu_eval_category_densities(object, grid, parallel, n_cores)`: Computes
  category likelihoods $p(x \mid c)$ across grid points using `parallel`.
- `.mvbu_eval_categorization_posteriors(object, grid, parallel, n_cores)`:
  Computes posterior probabilities $p(c \mid x)$ and decision choices across
  grid points using `parallel`.
- `.mvbu_build_1d_density_layers()` & `.mvbu_build_2d_density_layers()`:
  Standardized ggplot layer constructors for contour lines, quantile fill
  bands, gradient fills, and ellipses.
- `.mvbu_build_sliced_3d_layers()`: Standardized 2D sliced facet layers for
  3-cue spaces ($D = 3$).
- `.mvbu_build_parameter_intervals()`: Parameter mean and highest density
  interval (HDI) extraction and visualization.

### 2. `R/S7-plot-methods.R` (Core S7 Plotting Methods)
High-level methods for `MVBU_CategoryRepresentation`,
`MVBU_CategoryRepresentationTemplate`, and `MVBU_CognitiveModel`:
- `plot_categories(x, ...)`: 1D density curves, 2D contours / filled bands,
  and 3-cue sliced facet panels (default for $D = 3$).
- `plot_categorization(x, ...)`: 1D probability curves, 2D decision boundary
  maps, and 3-cue sliced facet panels (default for $D = 3$).
- `plot_parameters(x, ...)`: Category means, covariances, prior belief
  intervals, and exemplar distributions.

### 3. `R/S7-plot-stanfit-methods.R` (Stanfit S7 Plotting Methods)
Methods for `MVBU_Stanfit` and `IdealAdaptorStanfit`:
- `plot_categories(x, ...)`: Posterior predictive category densities across
  MCMC draws (with parallelization support).
- `plot_categorization(x, ...)`: Posterior predictive categorization
  functions across draws (with parallelization support).
- `plot_parameters(x, ...)`: Marginal and joint posterior parameter
  distributions.
- `plot_correlations(x, ...)`: Correlation matrices and pairwise parameter
  relationships.

### 4. `R/S7-plot-3d.R` (Interactive 3D Visualization Capabilities)
Interactive 3D visualization generics and S7 methods:
- `plot_categories_interactive(x, ...)`:
  - 2 cues: 3D density surface $z = p(x_1, x_2 \mid c)$ via `plotly`.
  - 3 cues: 3D ellipsoids and isosurfaces in 3D cue space via `plotly`.
- `plot_categorization_interactive(x, ...)`:
  - 2 cues: 3D posterior probability surface $z = p(c \mid x_1, x_2)$ via
    `plotly`.
  - 3 cues: 3D decision boundaries and probability volumes via `plotly`.

---

## Parallelization Strategy

- Standardized exclusively on the `parallel` package (built-in base R).
- Evaluates grid chunks in parallel using `parallel::mclapply` (Unix/macOS)
  or `parallel::parLapply` (Windows fallback) when `parallel = TRUE` or
  `n_cores > 1`.
- Defaults to `parallel = FALSE` for interactive single-core execution and
  deterministic reproducibility.

---

## Verification & Testing Plan

1. **Unit Testing:**
   - Existing plot test suites (`test-18-plot-categories.R`,
     `test-19-plot-categorization.R`, `test-20-plot-parameters.R`,
     `test-21-plot-stanfit.R`).
   - New test suite `test-22-plot-parallel-and-3d.R` testing:
     - Parallel vs sequential equivalence (`parallel = TRUE` vs
       `parallel = FALSE`).
     - 3-cue 2D slicing defaults for `plot_categories()` and
       `plot_categorization()`.
     - Interactive 3D plotting functions (`plot_categories_interactive()`,
       `plot_categorization_interactive()`).
2. **Quality & Standards:**
   - Strict $\le 80$ characters per line limit across all files.
   - Re-generate documentation with `devtools::document()`.
   - Validate with `tools::checkRd` (0 issues).
