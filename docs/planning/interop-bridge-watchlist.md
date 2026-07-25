# Interop Bridge Watchlist

Purpose: track which rstan/tidybayes-facing methods should be bridged for ModelDistribution classes, and why.

Policy:
- Prioritize methods used in real package workflows.
- Align syntax/behavior with rstan and tidybayes where relevant.
- Avoid adding new dependencies unless there is a clear maintenance and user-value win.
- Resolve items incrementally during Phases 1-3.

## Baseline (already agreed)
- [ ] get_stanfit(x)
- [ ] as_stanfit(x)
- [ ] get_draws(x, ...)
- [ ] summary(x, ...)
- [ ] print(x, ...)
- [ ] loo(x, ...)
- [ ] posterior::as_draws_df bridge

## Tidybayes bridge candidates
- [ ] Draw extraction helpers used in current plotting pipeline
- [ ] Additional helpers discovered by failing tests during migration

## Decision Log Template

For each candidate method:
- Method/signature:
- Used by (file/function/workflow):
- User-facing value:
- Dependency impact:
- Performance impact:
- Decision: adopt / defer / reject
- Phase decided:

## Exit Criteria
- [ ] All methods required by package tests and vignettes are bridged or deliberately deferred.
- [ ] No bridge added without recorded rationale.
- [ ] Dependency additions (if any) documented in DESCRIPTION rationale and NEWS entry when applicable.
