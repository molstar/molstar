# Drift: who coarsens when what is in view passes the store's bound

**What**: protocol.md (Lifetime and eviction) says that when what is in view alone exceeds the byte bound, "the store
asks for the next coarser level for the widest layers" and reports it once. The code's store reports it
(`overBudget`, `seq.store.over_budget`) and never asks a coarser level itself; nothing reads `overBudget` to coarsen.
What bounds a layer's reads today is before the store asks: `maxTilesPerLayer`, `maxBytesPerLayer` and `maxTileBytes`,
with a region layer reading the first coarser level that fits (`fittedLevel`) or saying "zoom in".

**Where**: `packages/seq/core/src/data/store.ts` (`overBudget`); the core spec, Store, "Eviction" ("a coarser level is
the planner's") and "A layer's wants together".

**Options**: (a) the design follows the code: the store evicts and reports, the per-layer bounds keep the view's tiles
within the cache before it is asked, and a coarser level is the planner's; (b) build the design: on `overBudget`, the
render planner reads the widest layers' next coarser level.

**Recommendation**: (a). The per-layer bounds (tiles, read cost, tile bytes) bound a view's reads before they reach the
cache; a planner reacting to `overBudget` would be a second, later path to the same end. Revisit if a view is ever
seen over budget with its layers each within their bounds.

**Resolved** (the author, 2026-11-23, as recommended): option (a): protocol.md's Lifetime and eviction says the store evicts and reports and a coarser level is the planner's (per-layer bounds), region-reads.md's second-grid text likewise; the rows point at core › Eviction and core › A layer's wants together.
