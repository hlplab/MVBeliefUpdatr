# Phase 6: Performance and Caching Framework Plan

## Purpose

This document captures the proposed Phase 6 strategy for making posterior and categorization paths faster and more stable within the new S7 architecture. The goal is not merely to cache results ad hoc, but to make the cached posterior objects genuinely reusable computation kernels that stay aligned with the current model state.

## Guiding Principle

The long-term target is to treat each treatment-specific posterior configuration as a compiled kernel that is built once from the model's immutable or slowly changing properties and then reused across many calls. The kernel should not repeatedly rebuild the posterior pipeline, re-parse the category template, re-derive priors, or re-prepare family-specific likelihood machinery on every invocation.

## Proposed Architecture

### 1. Split model state from compiled kernels

The cognitive model should continue to own the mutable model state:

- category prior
- lapse rate and lapse bias
- noise behavior and noise treatment
- category template and representation ordering
- decision rule
- metadata and label information

The cached posterior object should instead own a compiled kernel that is built from the current model state and the selected treatment configuration.

This creates a clean separation between:

- the stateful model object, and
- the stateless or mostly stateless execution kernel used to compute posteriors.

### 2. Use versioned invalidation

Every mutation to model properties should update a model state version. The compiled posterior kernel should record the state version it was built from.

Suggested design:

- add a lightweight internal version counter to the model object
- each mutation path increments the version
- each cached kernel stores the version it was compiled with
- before reuse, the kernel checks whether the current model version matches its stored version
- if not, it is invalidated and rebuilt

This keeps the cache correct without forcing eager recompilation after every minor update.

### 3. Compile kernels around immutable or precomputed pieces

For each treatment combination such as:

- `(noise_treatment = "no_noise", lapse_treatment = "no_lapses")`
- `(noise_treatment = "sample", lapse_treatment = "sample")`
- `(noise_treatment = "marginalize", lapse_treatment = "marginalize")`

build a kernel that pre-binds:

- the ordered category names
- the prior vector aligned to the current category order
- the relevant lapse parameters
- the relevant noise parameters
- the family-specific likelihood adapters or closures
- any expensive precomputed quantities such as:
  - Cholesky factors or covariance decompositions
  - log-determinants
  - precomputed log-prior vectors
  - precomputed category index maps for subsetting

The purpose is for the kernel to be able to evaluate the posterior matrix for new observations with minimal setup overhead.

### 4. Keep the public API unchanged

The public-facing methods should remain intuitive:

- `posterior()`
- `categorize()`
- `get_category_posterior_function()`

Internally, these should route through the compiled kernel cache rather than rebuilding the entire posterior path each call.

The change should be mostly internal and should not force callers to manage the cache explicitly.

### 5. Keep subsetting and batching lightweight

Category subsetting should not force recompilation of the whole kernel. Instead, the kernel can either:

- store an index map for the current category ordering and apply it cheaply, or
- support a `categories` argument that performs a light column selection after the main posterior computation.

List or batch inputs should still work through the existing outer methods, but the inner per-batch evaluation should reuse the same compiled kernel object.

### 6. Treat mutation paths as cache invalidators

Any method that changes model properties should be treated as a cache invalidation boundary. This includes:

- updating the category prior
- changing lapse rate or lapse bias
- changing noise treatment or noise variance
- replacing the category template or category representations
- switching the decision rule if it changes the computation pipeline

The implementation should make this explicit so that stale kernels are never silently reused.

## Suggested Implementation Shape

### Option A: lightweight kernel cache inside the model object

A simple first version would store a named list of compiled kernels inside the cognitive model object, similar to the existing `category_posterior_functions` field, but with a more explicit structure.

Each entry would contain:

- the treatment key
- the state version used for compilation
- the compiled function
- optionally, small metadata such as category ordering and prior vector

This is the lowest-risk implementation path and fits the current S7 architecture well.

### Option B: dedicated kernel class

A more explicit future version would introduce a dedicated S7 class such as:

- `MVBU_PosteriorKernel`

This class would hold:

- the treatment key
- the state version
- the compiled function
- the category index ordering
- the precomputed prior vector
- any family-specific precomputed parameters

This is more explicit and easier to evolve, especially if other expensive operations also need similar caching semantics later.

For the current Phase 6 scope, Option A is probably the best starting point, with Option B as a future refactor if the kernel interface grows.

## Initial Scope for Phase 6

Phase 6 should focus first on the posterior path and the categorization path, because these are the most likely to benefit immediately from reduced setup overhead.

Priority targets:

1. posterior computation for single observations and batches
2. categorization helpers that repeatedly call posterior internally
3. treatment-specific posterior kernels for noise and lapse settings
4. invalidation and version checks for model mutation paths

## Later Extensions

Once the core kernel cache works, the same pattern can be extended to:

- prediction and posterior-predictive helpers
- plotting and evaluation grids
- draw extraction and summary caches
- larger ModelDistribution-level caches

## Acceptance Criteria

Phase 6 can be considered successful when:

- posterior calls reuse precompiled kernels instead of rebuilding the full workflow each time
- model mutations invalidate stale kernels correctly
- the behavior remains identical to the current S7 posterior semantics
- the implementation does not introduce significant memory growth for ordinary use
- the caching strategy is explicit, testable, and easy to extend

## Implementation Order

1. Add a model state version counter and cache invalidation hooks.
2. Refactor the existing posterior-function cache into a more explicit compiled-kernel cache.
3. Pre-bind category order, priors, and treatment-specific constants.
4. Keep the public API unchanged while routing through the compiled kernel path.
5. Add regression tests for cache invalidation and treatment-specific reuse.
6. Benchmark the posterior path before and after the refactor.
