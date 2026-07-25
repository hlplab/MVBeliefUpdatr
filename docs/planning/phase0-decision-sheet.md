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
- [x] testthat rewrite/expansion occurs in each phase, tied to phase scope
- [x] Obsolete code is removed promptly when superseded
- [x] Obsolete files are deleted or moved to deprecated only when transitional value exists
- [x] Keep repository tidy continuously (no long-lived dead paths)
- [x] NEWS.md and DESCRIPTION are updated whenever a phase lands meaningful API/behavior/package changes

## Open Decisions

### 1) Inferred-model interoperability with rstan/tidybayes
Goal: inferred objects remain first-class package objects while staying easy to use with rstan/tidybayes.

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
- [x] Allow opt-in persistent caches for inferred models
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
- [ ] Architecture contract is written and approved
- [ ] Canonical class names and hierarchy are fixed
- [ ] Generic/method naming conventions are fixed
- [ ] Wrapper/deprecation strategy text is finalized
- [ ] Cache policy and defaults are finalized
- [ ] Per-phase tests/docs/cleanup policy is finalized
- [ ] NEWS/DESCRIPTION/versioning policy is finalized
- [ ] Definition of done for each later phase is documented

## Definition of Done Template (for every phase)
A phase is done only if all apply:
- [ ] Implementation completed for phase scope
- [ ] testthat tests added/updated for phase scope
- [ ] Documentation updated for phase scope (Rd + vignette/README as relevant)
- [ ] Obsolete code removed or moved to deprecated-2026.R (do not expand deprecated-2025.R)
- [ ] Obsolete files deleted or archived with justification
- [ ] NEWS.md updated if triggered
- [ ] DESCRIPTION + version updated if triggered
- [ ] Migration notes updated when user-facing changes occur

## Initial Risk Controls
- Keep wrappers thin (no business logic in deprecated-2026.R)
- Treat deprecated-2025.R as legacy-only and avoid adding new migration logic there
- Add parity tests across NIW/MVG/exemplar/inferred object families
- Add memory guards for cache-heavy pathways
- Add benchmark checks for known hotspots before and after refactors

## Immediate Next Actions
1. Finalize canonical class names and generic names (Phase 0 architecture contract).
2. Finalize inferred-model interop list: exactly which rstan/tidybayes generics are forwarded.
3. Draft first versioned migration table old API -> new API namespace.
4. Define first phase-scoped cleanup candidates (files/functions likely to become obsolete).