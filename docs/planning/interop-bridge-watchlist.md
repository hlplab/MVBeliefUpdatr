# Interop Bridge Watchlist

Purpose: track which rstan/tidybayes-facing methods should be bridged for ModelDistribution classes, and why.

Policy:
- Prioritize methods used in real package workflows.
- Align syntax/behavior with rstan and tidybayes where relevant.
- Avoid adding new dependencies unless there is a clear maintenance and user-value win.
- Resolve items incrementally during Phases 1-3.

## Baseline (already agreed)
- [x] get_stanfit(x) - Adopted (R/S7-stanfit-methods.R)
- [x] as_stanfit(x) - Rejected in favor of canonical `get_stanfit(x)` to avoid redundant accessors
- [x] get_draws(x, ...) - Adopted (R/S7-stanfit-methods.R)
- [x] summary(x, ...) - Adopted (R/S7-stanfit-methods.R)
- [x] print(x, ...) - Adopted (R/S7-stanfit-methods.R)
- [x] loo(x, ...) - Adopted (R/S7-stanfit-diagnostics.R)
- [x] posterior::as_draws_df bridge - Adopted (as_draws_df, as_draws_matrix, as_draws_list, as_draws_rvars in R/S7-stanfit-methods.R)

## Tidybayes and tidyverse bridge candidates
- [x] Draw extraction helpers used in current plotting pipeline - Adopted via `get_draws()` and `posterior::as_draws_*`
- [x] `recover_types()` integration during model fitting - Adopted (R/S7-stanfit-fitting.R)
- [x] Legacy tibble bridge: `as_tibble()` methods and dplyr verb forwarding - Adopted (R/S7-to-legacy-tibble-compatibility.R)

## Decision Log

1. `get_stanfit(x)`: Adopted in Phase 3. Canonical S7 generic and method returning underlying `stanfit` object.
2. `as_stanfit(x)`: Rejected. Redundant with `get_stanfit()`.
3. `get_draws(x, ...)`: Adopted in Phase 3. Canonical S7 generic returning posterior draws in long format.
4. `summary(x, ...)` & `print(x, ...)`: Adopted in Phase 3. Returning `Summary_MVBU_Stanfit` with parameter summaries.
5. `loo(x, ...)`: Adopted in Phase 3. S7 generic and method wrapping `loo::loo`.
6. `posterior::as_draws_*`: Adopted in Phase 3. Full method suite for seamless interop with `posterior` and `bayesplot`.
7. `as_tibble(x)` & `dplyr::*`: Adopted in Phase 4. Backward-compatibility bridge for legacy scripts with deprecation warnings.

## Exit Criteria
- [x] All methods required by package tests and vignettes are bridged or deliberately deferred.
- [x] No bridge added without recorded rationale.
- [x] Dependency additions (if any) documented in DESCRIPTION rationale and NEWS entry when applicable.
