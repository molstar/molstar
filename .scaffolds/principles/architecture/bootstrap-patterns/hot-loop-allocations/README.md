# Zero-allocation hot loops {#hot-loop-allocations}

- **Implementations:** [TypeScript](typescript.md)

Rules §no-accidental-quadratics to §zero-alloc-hot-loops: no allocations in steady-state render and hot loops; minimize them elsewhere, with scratch arenas for temporary memory; no accidental quadratics; structures chosen for the access pattern. Models may allocate when state *changes*; the loop only reads preallocated storage.

## Problem

A render or simulation loop that allocates per frame (a filtered copy, a closure, one small record per item) produces
garbage at frame rate, and collection pauses show up as dropped frames. A lookup inside a loop over another collection
turns linear work quadratic, invisible at test sizes and fatal at real ones. Neither is caught by good intentions; only
budgeted tests catch them (§performance-tested).

## Shape

- **Allocate in setup, reuse in the loop.** Numeric data iterated in bulk is stored as struct-of-arrays (one flat
  numeric buffer per field, a `count`, a fixed capacity); a step updates in place and compacts survivors in one pass.
- **Out-parameters**: a function writes into a caller-supplied buffer and returns it (`normalizeInto(out, x, y, z)`).
- **Pool**: `acquire() -> item` reuses a released item or makes a new one; `release(item)` resets it and keeps it.
- **Scratch arena**: one preallocated block; temporary memory is carved by advancing an offset, and freed by moving
  the offset back. Allocation is O(1), there is no per-object garbage, and a whole scope of temporaries goes at once.
  Use it for per-call or per-frame working memory (intermediate vertex data, sort keys, neighbor lists, parse buffers).
  - `alloc(n) -> offset` into the arena's one buffer (offsets, not new views, where creating a view allocates);
  - `mark() -> mark`, `release(mark)`, `reset()` (release to 0);
  - `peak`: the highest use since creation, the number a guard test budgets;
  - a debug mode poisons released memory so a use after release shows in tests.
- **Per-frame arena**: reset once per frame; passes allocate from it and free nothing individually. The arena is
  created by the runtime or composition and passed in, never a module-level singleton (§dependencies-passed-in).
- **Structures for the access pattern** (§efficient-structures):

| Access pattern | Structure |
|---|---|
| membership or lookup by key, many queries | hash set or map |
| ordered iteration, few lookups, small N | array (N ≤ ~16: a scan beats hashing) |
| ordered with range queries, frequent inserts | balanced tree, or sorted array with binary search |
| numeric data iterated in bulk | flat numeric arrays, struct-of-arrays |
| FIFO with many pops | ring buffer (removing from an array's front is O(n)) |
| repeated "next min/max" | binary heap, not a re-sort each time |
| incremental recompute | dirty flags or a keyed diff, not a full rebuild |

A knowingly simple choice says so where it matters: `N <= 8 (max open panes); linear scan is intentional.`

## Invariants

- A steady-state frame or hot-loop iteration allocates nothing; models may allocate when state changes.
- Nothing from an arena escapes its scope (stored in a model, returned, captured by a closure); results go into
  caller-supplied buffers or are copied to long-lived storage; an offset is never returned past its mark.
- Nested arena use is a stack: release to the mark you took, in reverse order, also on error.
- One arena per thread or lane; an arena is never held across a suspension point (an await).
- Arena capacity is a budgeted, tested number; it never grows on demand in a hot path; exhaustion fails loudly and says
  to raise the capacity.
- `alloc` is O(1) and allocates nothing.
- For long-lived or irregular lifetimes, use pools or normal allocation, not an arena.
- No lookup over a collection inside a loop over another unless N is provably small and a comment says so.
- Allocation budgets and scaling are checked by guard-tier tests (§test-tiers).
- A per-frame path at rest does nothing: its guard compares against the rest case, where a check that is always true
  still costs the frame. A baseline that moves is bisected to the commit that moved it before it is raised; a raise
  names that commit and its reason.

## Conformance cases

- Arena of capacity 10: `alloc(4)` → 0; `alloc(4)` → 4; `alloc(4)` → fails (exhausted); `peak` → 8.
- Arena: `alloc(4)`; `m = mark()` → 4; `alloc(3)` → 4; `release(m)`; `alloc(2)` → 4.
- `reset()` after use → next `alloc` returns 0; `peak` is unchanged.
- Debug arena: `release(m)` → every slot from `m` to the previous top reads as NaN.
- `normalizeInto(out, 0, 0, 0)` → `(0, 0, 0)`; `normalizeInto(out, 3, 0, 4)` → `(0.6, 0, 0.8)`, written into `out`.
- Trail points with lives `[1, 0.05, 2]`, `step(0.1)` → `count` 2, lives `[0.9, 1.9]`, order kept.
- Pool: `release(a)`, then `acquire()` → `a`, reset.
- Guard: after warm-up, 1,000 steady-state frames grow the heap by less than 16 KiB.
- Guard: 100 frames on a per-frame arena keep `peak` below half its capacity.
- Guard: 4× the input takes less than 8× the time (best of 5 runs); quadratic would be ~16×.

## Open

- `release(mark)` does not check that the mark is at or below the current top, so a stale mark can move the top up;
  the code does not say whether that is a caller's bug to assert on.
