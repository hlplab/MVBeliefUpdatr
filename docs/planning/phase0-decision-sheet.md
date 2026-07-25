# Phase 0 Decision Sheet: MVBeliefUpdatr v1.0 Overhaul

## Purpose
Lock governance and architecture decisions before implementation.
This sheet is the contract for Phases 1+.

## Locked Decisions
- [x] Migration strategy: big-bang v1.0 redesign
- [x] Naming policy: aggressive renaming in new code
- [x] Backward compatibility: wrappers only, isolated in deprecated-2026.R (legacy historical wrappers remain in deprecated-2025.R and are out of current scope)
- [x] Wrapper horizon: maintain for 2-3 minor releases
- [x] Core modeling design: explicit composition (model contains representation)
- [x] Stan expansion timeline: after core class/API stabilization

## New Governance Rules (added)
- [x] Documentation rewrite occurs in each phase, not deferred to the end
- [x] All new code and all code integrated into the S7 scaffold must be roxygen documented within the phase that introduces/integrates it
- [x] End-of-phase cleanup must include documentation QA (roxygen generation, unresolved-link checks, malformed Rd/macros checks)
- [x] testthat rewrite/expansion occurs in each phase, tied to phase scope
- [x] Obsolete code is removed promptly when superseded
- [x] Obsolete files are deleted or moved to deprecated only when transitional value exists
- [x] Keep repository tidy continuously (no long-lived dead paths)
- [x] NEWS.md and DESCRIPTION are updated whenever a phase lands meaningful API/behavior/package changes

## Open Decisions

### 1) ModelDistribution interoperability with rstan/tidybayes
Goal: ModelDistribution objects remain first-class package objects while staying easy to use with rstan/tidybayes.

Options:
- [ ] A. Strict wrapper only (always call get_stanfit()/as_stanfit())
- [x] B. Wrapper + S3 method forwarding for key generics
- [ ] C. On-demand adapter object for external tooling

Current recommendation: B
Reason: best ergonomics with controlled internal structure.

### 2) Cache behavior defaults
Goal: speed where useful, avoid unexpected object mutation and memory bloat.

Decisions:
- [x] Base policy: pure-by-default methods
- [x] Allow opt-in persistent caches for ModelDistribution objects
- [x] Do not store large plot-grid/surface results inside core model objects
- [x] Add dedicated grid/result objects for large plotting/intermediate computations
- [x] Provide profile options: minimal/standard/eager
- [x] Allow per-call overrides and global options

Default profile proposal:
- [ ] minimal
- [x] standard
- [ ] eager

## Release Hygiene Policy (NEWS, DESCRIPTION, Version)
Apply this at the end of every phase that changes behavior.

### Trigger to update NEWS.md
Update NEWS.md when any of the following happen:
- API added/removed/renamed
- behavior changed (even if signature is unchanged)
- performance changed materially
- deprecations introduced/removed
- tests/docs overhauled in user-visible ways

### Trigger to update DESCRIPTION
Update DESCRIPTION when any of the following happen:
- Version bump
- Dependency changes (Imports/Suggests/LinkingTo)
- Metadata changes (Title/Description/URLs/Authors)
- Minimum R version changes

### Version bump convention for this overhaul
- Patch (x.y.Z): bug fixes only, no API change
- Minor (x.Y.z): additive API changes + deprecations
- Major (X.y.z): breaking changes, class redesign, incompatible semantics

Current intent:
- continue on development version line until phase gate for v1.0 readiness is met
- switch to v1.0.0 when compatibility layer + migration docs + parity tests pass

## Phase 0 Acceptance Criteria (Gate)
Phase 0 is complete only when all are checked:
- [x] Architecture contract is written and approved
- [x] Canonical class names and hierarchy are fixed
- [x] Generic/method naming conventions are fixed
- [x] Wrapper/deprecation strategy text is finalized
- [x] Cache policy and defaults are finalized
- [x] Per-phase tests/docs/cleanup policy is finalized
- [x] NEWS/DESCRIPTION/versioning policy is finalized
- [x] Definition of done for each later phase is documented

## Definition of Done Template (for every phase)
A phase is done only if all apply:
- [ ] Implementation completed for phase scope
- [ ] testthat tests added/updated for phase scope
- [ ] Documentation updated for phase scope (Rd + vignette/README as relevant)
- [ ] Roxygen documentation exists for all newly added functions/classes and all code integrated into the S7 scaffold in phase scope
- [ ] Roxygen generation and documentation QA checks pass for phase-touched files (no new broken links; no malformed Rd issues)
- [ ] Obsolete code removed or moved to deprecated-2026.R (do not expand deprecated-2025.R)
- [ ] Obsolete files deleted or archived with justification
- [ ] NEWS.md updated if triggered
- [ ] DESCRIPTION + version updated if triggered
- [ ] Migration notes updated when user-facing changes occur

## Initial Risk Controls
- Keep wrappers thin (no business logic in deprecated-2026.R)
- Treat deprecated-2025.R as legacy-only and avoid adding new migration logic there
- Add parity tests across NIW/MVG/exemplar/ModelDistribution object families
- Add memory guards for cache-heavy pathways
- Add benchmark checks for known hotspots before and after refactors

## Immediate Next Actions
1. Kick off Phase 1 S7 foundations (base classes, generics, validator stubs, and baseline tests).
2. Initialize interop bridge watchlist items in implementation and resolve incrementally via workflow tests.
3. Draft first versioned migration table old API -> new API namespace.
4. Define first phase-scoped cleanup candidates (files/functions likely to become obsolete).