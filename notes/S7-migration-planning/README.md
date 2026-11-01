# Planning Index

This folder contains the overhaul planning artifacts for the MVBeliefUpdatr v1.0 architecture migration.

## Suggested Reading and Working Order

1. [Phase 0 Decision Sheet](phase0-decision-sheet.md)
- Governance rules, release hygiene, deprecation policy, and phase-level definition of done.
- Use this to lock process decisions first.

2. [Phase 0 Architecture Contract](phase0-architecture-contract.md)
- Canonical class hierarchy, composition rules, generic naming direction, interop and cache contracts.
- Use this to lock technical architecture before implementation.

3. [Phase 0 API and Class Naming Map](phase0-api-name-map.md)
- Old API to proposed v1 API mapping for classes, constructors, generics, and wrappers.
- Use this to quickly approve naming and unblock implementation.

4. [Editable Overhaul Plan](plan-v1-overhaul-editable.md)
- Phase-by-phase execution plan with checklists and phase gates.
- Use this as the active implementation tracker.

5. [Legacy To-Do (Migrated)](legacy-todo.md)
- Historical backlog migrated from the old plain-text list.
- Use this as an idea bank and triage source, not as the primary execution tracker.

6. [Interop Bridge Watchlist](interop-bridge-watchlist.md)
- Tracks rstan/tidybayes bridge-method decisions and dependency rationale during Phases 1-3.
- Use this to finalize bridge coverage incrementally without over-committing dependencies early.

## Maintenance Rules for This Folder

- Keep planning files current as phase decisions change.
- If a plan item is superseded, either update in place or archive with explicit rationale.
- Keep naming and phase numbering synchronized across all documents.
- Reflect major planning changes in [NEWS.md](../../NEWS.md) when they result in user-visible package changes.

## Current Source of Truth

- Process and governance: [phase0-decision-sheet.md](phase0-decision-sheet.md)
- Technical contract: [phase0-architecture-contract.md](phase0-architecture-contract.md)
- Naming decisions: [phase0-api-name-map.md](phase0-api-name-map.md)
- Interop bridge decisions: [interop-bridge-watchlist.md](interop-bridge-watchlist.md)
- Execution tracking: [plan-v1-overhaul-editable.md](plan-v1-overhaul-editable.md)