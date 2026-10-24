# Slice 03: core state and task (v2)

Scope: `packages/core/src/state` (tree, transient tree, builder, reconciler, diff, selections) and
`packages/core/src/task` (observable runtime, scheduler, `chunkedSubtask`), plus `core/src/util/async-queue.ts`. The
slice is not on the profiled load path: cost scales with the number of state nodes n, not atoms. Measured at
`0944823d4` with Node 25.6 (M4), trivial transformers, 10 `stateUpdated` subscribers. Nothing here touches frozen gl or
Canvas3D code. MVS frame building is covered by [.v6/mvs-optimizations/research](../mvs-optimizations/research/README.md);
P15 is the core side of applying a frame, which that research lists as not measured.

Since v1: #2013 (14/F4) diffs a builder's edits instead of both trees (`diff.ts`); #1989 fixed the flag loss (B1);
#2014 fixed the Immutable shared-entry write (14/N7); #2022 (S9) caches the default param schema per transformer (P7)
and gave `setSnapshot` update options. The task runtime is unchanged.

## Hot paths now

| Path | Cost (10k nodes unless noted) | Notes |
| --- | --- | --- |
| Builder commit, edits path (≤1 removal, no move, not rebased, cells in sync) | 0 edits 0.002 ms, 1 update 0.04 ms, 2 updates 0.10 ms, add 2 nodes 0.03 ms, 1 removal 0.02 ms | was 5.2-7.5 ms in v1 |
| Builder commit with 2+ removals, a move, or a rebase | 0.55 ms at 1k, 5.9 ms at 10k, 12.4 ms at 20k | `findDeletes` + `findUpdateRoots`, about 0.6 µs/node (P1) |
| A `StateTree` passed in (undo, transaction revert) | 8-9 ms | whole-tree diff + `syncNewStates` |
| `setSnapshot`, one frame changing one node | 1.4-1.9 ms at 1k, 20-23 ms at 10k, 46 ms at 20k | `fromJSON` (10-12 ms) + whole-tree diff + `syncNewStates` (P15) |
| Whole-tree walks (`doPreOrder`, `ofTransformer`, `findTagInSubtree` miss, `getSnapshot`) | 2.4 / 2.9 / 2.4 / 2.8 ms | 2 Immutable lookups per node (P2) |
| Builder commit with a flag written while it was outstanding | 0.31 ms at 1k, 2.6 ms at 10k flagged refs | scans the whole flag log (P17) |

Measured and cheap: `stateUpdated` (3 per evaluated node, about 0.2 µs per event with 10 listeners);
`getEffectiveDependsOn`'s structural scan for the 7 plugin transformers without a default schema (about 1 µs per call);
retained param definitions (about 4 KB per cell, 1.9 MB for 491 cells); a transaction around one fast-track edit
(0.02 ms).

## v1 → now

Counts: 21 v1 items (P1-P9, B1-B10, and the core-state items 14/D2, 14/N7 of slice 14): fixed 3, partly fixed 3, still
open 15, regressed 0, dropped 0. The B1 fix introduced a new bug (B11).

| v1 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| P1 | Every commit except a one-edit builder diffs both trees | partly fixed | #2013 | Edits path for ≤1 removal; 2+ removals, moves, rebases, trees passed in still whole-tree (P1 below) |
| P2 | Whole-tree walks: 2 Immutable lookups per node, about 10 per change | partly fixed | #2013, #2008 | Fewer plugin walks per change; core walk unchanged, no pre-order index |
| P3 | 22-char refs thrash Immutable's string-hash cache | still open | | Also 86-95 % of counter refs collide in the 31-hash |
| P4 | `collectScheduled` repeats the `initCellStatus` walk | still open | | [state.ts:1285](../../packages/core/src/state/state.ts#L1285) |
| P5 | No batched `updateCellState` | still open | | per-ref, one `syncVisibility` per representation |
| P6 | `AsyncQueue` broadcasts every dequeue: O(k²) | still open | | [async-queue.ts:24](../../packages/core/src/util/async-queue.ts#L24) |
| P7 | `tryGetDefaultParamDefinition` throws once per transform | fixed | #2022 | `WeakMap` per transformer ([transform.ts:188](../../packages/core/src/state/transform.ts#L188)); dependencies.md §3 records the no-registry rule |
| P8 | `setImmediate` polyfill prefers `window.postMessage('*')` | still open | | [scheduler.ts:160](../../packages/core/src/task/util/scheduler.ts#L160) |
| P9 | Builder `updateState` scans actions with `find` | still open | | latent |
| B1 | State flags set while a builder is outstanding never reach `state.tree` | fixed | #1989 | `asImmutable` adopts its result; `flagLog` + `carryFlagWrites`; re-verified (`verify/v-b1.ts`). Fix introduced B11 |
| B2 | Builder `updateState` never reaches `cell.state`, re-evaluates | still open | | reproduced |
| B3 | Rebase replays `tag`/`dependsOn`/`updateState` as `setParams` | still open | | code unchanged |
| B4 | `setDependsOn` skips `tree.dependencies` | still open | | code unchanged |
| B5 | Hash versions cover params only | still open | | reproduced |
| B6 | `doError` clears objects without `dispose` | still open | | reproduced |
| B7 | Late abort discards a completed result | still open | | code unchanged |
| B8 | `chunkedSubtask` sizing compounds | still open | | code unchanged |
| B9 | `subtree()`/`children()` return undefined cells | still open | | reproduced |
| B10 | `findTagInSubtree` never stops, returns the last match | still open | | reproduced |
| 14/D2 | `findDeletes` + `findUpdateRoots` per multi-edit commit (many structures) | partly fixed | #2013 | same as P1 |
| 14/N7 | A builder's removal writes into `state.tree` (Immutable HashCollisionNode) | fixed | #2014 | `change*`/`delete*` gates, `renewOwner`; canary + random tests in `transient.test.ts` |

## Open performance findings (ranked)

| Rank | Id | vs v1 | Title | File | Kind | Impact | Estimate | Effort | Verdict | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | P1 | partly fixed; **fixed** #2025 | Commits removing 2+ nodes (or moving, or rebased) diff both trees: O(N²) per-structure | [state.ts:646](../../packages/core/src/state/state.ts#L646) | time | high | preset re-apply over N=500 structures 1.4 s → a few ms of diff; N=1000 6.1 s → ~10 ms; **measured**: `many --count 500` `represet` (the preset menu over all) 4.1–15.1 → 0.5–0.7 s, alloc 3.08 → 0.37 GB; `remove:all` unchanged | M | confirmed | no |
| 2 | P15 | new | Every snapshot frame rebuilds and diffs the whole tree | [state.ts:203](../../packages/core/src/state/state.ts#L203) | time | medium | 20-23 ms → ~0.1-0.3 ms per frame at 10k nodes (1.4-1.9 ms at 1k) | M | confirmed | no |
| 3 | P2 | partly fixed | Whole-tree walks: 2 Immutable lookups per node | [immutable.ts:134](../../packages/core/src/state/tree/immutable.ts#L134) | time | medium | 2.4-2.9 ms per walk at 10k; 20-50x with a memoized pre-order index | L | confirmed | no |
| 4 | P17 | new; **fixed** #2046 | `carryFlagWrites` scans every ref ever flagged | [state.ts:117](../../packages/core/src/state/state.ts#L117) | time | low | 2.6 ms per queued commit at 10k flagged refs, 0.31 ms at 1k | S | confirmed | no |
| 5 | P16 | new; partly fixed #2034 (focus callers, 12/P7) | `update(recipe)` returning a new object: mutative walks the returned params | [builder.ts:341](../../packages/core/src/state/state/builder.ts#L341) | time | low | 0.07 / 0.63 / 2.25 ms for 1k / 100k / 1M-int bundles vs 0.002-0.01 ms mutating | S | confirmed | no |
| 6 | P3 | still open | Counter refs: 22 chars (string-hash cache) and 31-hash collisions | [transform.ts:114](../../packages/core/src/state/transform.ts#L114) | time | low | Immutable get 153-161 vs 125 ns (random 22) vs 50 ns (16 chars); end to end within noise | S | confirmed | no |
| 7 | P4 | still open; **fixed** #2046 | `collectScheduled` repeats the `initCellStatus` walk | [state.ts:1285](../../packages/core/src/state/state.ts#L1285) | time | low | one walk of the evaluated subtrees, ~2.5 ms at 10k evaluated nodes | S | confirmed | no |
| 8 | P5 | still open | No batched `updateCellState` | [state.ts:212](../../packages/core/src/state/state.ts#L212) | time | low | k `syncVisibility` → 1 for hide-all | M | plausible | no |
| 9 | P6 | still open; **fixed** #2046 | `AsyncQueue` broadcasts every dequeue | [async-queue.ts:24](../../packages/core/src/util/async-queue.ts#L24) | time | low | O(k²) → O(k) callbacks | S | confirmed | no |
| 10 | P8 | still open | `setImmediate` polyfill uses `postMessage('*')` first | [scheduler.ts:160](../../packages/core/src/task/util/scheduler.ts#L160) | time | low | one stray `message` event per yield | S | confirmed | no |
| 11 | P9 | still open; **fixed** #2046 | Builder `updateState` uses `actions.find` | [builder.ts:363](../../packages/core/src/state/state/builder.ts#L363) | time | low | latent | S | confirmed | no |

### P1: multi-removal commits take the whole-tree diff

- Mechanism: `editsOnly` needs `edits.removed.length <= 1 && !edits.moved` and an un-rebased builder
  ([state.ts:641-647](../../packages/core/src/state/state.ts#L641)). `removed` also holds dependents a removal cascaded
  to, so one cascading delete is already 2+. Everything else runs `findDeletes` (post-order of the old tree) and
  `findUpdateRoots` (pre-order of the new tree). Per-structure multi-removal commits in the plugin:
  `component.ts` `syncPreset` (each old component or representation not kept) and `clearComponents`, and the hierarchy
  manager's `remove`, `updateStructure`, `clearTrajectory`. `applyPreset` runs `syncPreset` for each of N structures,
  each commit O(total nodes).
- Evidence (`bench-commit.ts`, `multi-remove.ts`, re-run): 2 removals 0.55 / 5.9 / 12.4 ms at 1k / 10k / 20k; 1 removal
  0.02 ms. N structures of 13 nodes, one commit each: 2 removals 54 ms (N=100), 492 ms (300), 1.41 s (500), 6.06 s
  (1000); 1 removal 2-6 ms total. Pure state diff, no transformer work.
- Fix: generalize `findEditedDeletes` to k removals: dedup `removed`, drop refs absent from the old tree or with an
  ancestor in `removed`, sort the rest in old-tree pre-order (`pathOf`/`childOrder`/`comparePreOrder` on the old tree),
  concatenate each root's post-order subtree (equal to the whole tree's post-order restricted to disjoint subtrees).
  Relax `editsOnly` to any `removed.length`; keep the whole tree for moves and rebases. Extend `diff.test.ts`'s
  differential generator with multi-removal builders; update reconcile.md §2/§8.
- Risk: the order must stay the whole old tree's post-order (dispose and event order, reconcile.md §3/§7).

### P15: snapshot frames rebuild and diff the whole tree

- Mechanism: `setSnapshot` runs `StateTree.fromJSON` (a `StateTransform.fromJSON` and 3 Immutable sets per node;
  each child set is set twice, mutable then immutable) and `updateTree` with a plain tree: `findDeletes`,
  `findUpdateRoots`, `syncNewStates`, three O(n) walks. Callers per frame: `PluginState.setAnimationSnapshot`
  (scrubbing, [plugin state.ts:346](../../packages/plugin/core/src/state.ts#L346)), the built-in snapshot animation
  (`state-snapshots.ts:36`), transitions via `setSnapshotAnimationFrame`. Since #2022 once per applied frame, not per
  tick. Frame transforms are new objects, so unchanged nodes still pay hash and lookup.
- Evidence (`snap-frames.ts`, `bench-commit.ts`, re-run): one leaf changed per frame: 1.39 ms at 1k, 20.4-22.8 ms at
  10k, 46 ms at 20k; `fromJSON` alone 10.2-11.8 ms at 10k. Profile at 10k (`prof-snap/`): `fromJSON` ~45 %,
  `findDeletes` 19 %, `findUpdateRoots` + `syncNewStates` ~26 %; self time Immutable `HashArrayMapNode.get` 17 %, `is`
  9 %, `cachedHashString` 8 %, GC 11 %.
- Fix: a frame-delta path in `setSnapshot`. Compare the frame's serialized transforms with the current tree (or the
  last applied frame, kept per `State`) by ref, parent, version and flags in a native loop; when the shape is the same,
  edit a transient of the current tree in place (`fromJSON` only for changed transforms, `assignState` for flags),
  commit through the edits path (k roots) and sync cell flags for those refs only. Fall back to `fromJSON` + whole-tree
  diff when the shape differs. Deltas per transition can be cached on the frames array (as S4's `WeakMap`). Smaller:
  `fromJSON` builds child sets in a native `Map` and sets each once.
- Risk: keep tree-passed-in semantics (the tree's flags win, B11's log reset, deleted cells, current-object move); add
  a differential test against today's `setSnapshot`.

### P2: whole-tree walks

- Mechanism: `_doPreOrder`/`_doPostOrder` do `children.get` + `transforms.get` per node with 22-char keys. Still per
  change outside transactions: the structure hierarchy rebuild, measurement sync (1 walk), the focus representation's
  `byRef(root).subtree().withTag` selects, selection generators, `findTagInSubtree`, `findWithAllTags`; and inside the
  reconciler on every whole-tree pass (P1, P15).
- Evidence (`bench-commit.ts`, 10k / 20k): `doPreOrder` 2.45 / 4.96 ms; `ofTransformer` 2.9 / 5.7 ms;
  `findTagInSubtree` miss 2.4 / 4.6 ms; `getSnapshot` 2.8 / 5.3 ms. Slice 14: the hierarchy rebuild is 12 ms of each
  single action at N=500.
- Fix: (a) order-free searches by `transforms.forEach` or a scan of `state.cells`; (b) lazily memoize per immutable
  `StateTree` a pre-order ref array + subtree-end indices; a visitor returning `false` jumps to the end. Move the plugin
  hierarchy managers onto it.
- Risk: O(n) per retained tree (history keeps 5): build lazily; `syncNewStates` order stays pre-order.

### Low

- **P17**: `flagLog` (ref → flags) only shrinks on deletion; any builder commit with a flag write after its `build()`
  mark iterates the whole map with an Immutable `has` + `get` per entry (`flaglog.ts`: 2.59 ms at 10k flagged refs vs
  0.006 ms without a write; 0.31 ms at 1k). The scan runs when the queued commit executes, so it hits commits waiting
  behind a long update during UI toggles. Fix: an append-only write list ordered by `n`, scanned from the first
  `n > mark` (binary search), trimmed on B11's reset.
- **P16**: `To.update(transformer, recipe)` runs `produce(old.params, recipe)`; a recipe returning `{ ...old, bundle }`
  spreads the draft (creating child drafts), so mutative traverses the returned value, including the new bundle, to
  finalize them (`produce-bench.ts`: 0.07 / 0.63 / 2.25 ms for 1k / 100k / 1M ints, vs 0.002-0.010 ms mutating the
  draft). Callers: `structure-focus-representation.ts:212, 217, 271, 272` on every focus change. Fix: pass the params
  (`update({ ...cell.transform.params, bundle })`) or mutate the draft; document in `builder.update`.
- **P3**: `UUID.create22` differs only in its last chars; 31-hash pairs (`Aa`/`BB`) collide (`hash-coll.ts`: 86 % of
  interleaved refs, 95 % of consecutive). Immutable `get` 153-161 ns vs 125 ns random 22-char vs 50 ns 16-char; a 10k
  `setSnapshot` with 15-char random refs was within noise. Fix: ≤16-char refs with the counter through a bijective mixer.
  Mostly subsumed by P2/P15.
- **P4**: fill `ctx.scheduled` in `initCellStatusVisitor` and drop `collectScheduled`.
- **P5**: `setSubtreeVisibility` calls `updateCellState` per node; `UpdateRepresentationVisibility` calls
  `canvas3d.syncVisibility()` per representation. Add `State.updateCellStates` with a batch end, or coalesce in the
  behavior (no Canvas3D change). Not measured in a plugin.
- **P6**: `Map<T, resolve>`; `handled` resolves the new head, `remove` the removed entry.
- **P8**: prefer `MessageChannel` (or `scheduler.postTask`), keep `postMessage` as fallback.
- **P9**: a `Set` of refs with an update action, or a separate `state` action (with B2, B3).

## Bugs

| Id | vs v1 | Severity | Title | File | Verdict |
| --- | --- | --- | --- | --- | --- |
| B11 | new | medium | Carried flags override a tree passed in after them (setSnapshot, undo, revert): tree and cell disagree | [state.ts:467](../../packages/core/src/state/state.ts#L467) | confirmed |
| B6 | still open | medium | `doError` clears objects and caches without `dispose` | [state.ts:1257](../../packages/core/src/state/state.ts#L1257) | confirmed |
| B9 | still open; **fixed** #2034 | low | `subtree()`/`children()` return undefined cells; chained filters throw | [selection.ts:278](../../packages/core/src/state/state/selection.ts#L278) | confirmed |
| B12 | new; **fixed** #2033 | low | `applyOrUpdateTagged` drops extra tags when the params changed | [builder.ts:318](../../packages/core/src/state/state/builder.ts#L318) | confirmed |
| B5 | still open | low | Hash versions cover params only | [transform.ts:281](../../packages/core/src/state/transform.ts#L281) | confirmed |
| B10 | still open; **fixed** #2034 | low | `findTagInSubtree` never stops, returns the last match | [selection.ts:448](../../packages/core/src/state/state/selection.ts#L448) | confirmed |
| B7 | still open | low | A late abort discards a completed result, undisposed | [observable.ts:106](../../packages/core/src/task/execution/observable.ts#L106) | confirmed |
| B8 | still open; **fixed** #2045 | low | `chunkedSubtask` sizing compounds | [chunked.ts:35](../../packages/core/src/task/util/chunked.ts#L35) | confirmed |
| B3 | still open; **fixed** #2046 | low | Rebase replays `tag`/`dependsOn`/`updateState` as `setParams` | [builder.ts:89](../../packages/core/src/state/state/builder.ts#L89) | confirmed (mostly latent) |
| B4 | still open; **fixed** #2046 | low | `setDependsOn` skips `tree.dependencies` | [transient.ts:259](../../packages/core/src/state/tree/transient.ts#L259) | confirmed (latent) |
| B2 | still open; partly fixed #2046 (no evaluation, replayed on a rebase; the flags still reach the tree only) | low | Builder `updateState` never reaches `cell.state`, re-evaluates | [builder.ts:357](../../packages/core/src/state/state/builder.ts#L357) | confirmed (latent) |

- **B11** (new, from #1989): `updateCellState` logs each flag write with a sequence number; `carryFlagWrites` writes
  every logged flag newer than the builder's `build()` mark into the builder's tree. A plain tree passed to `updateTree`
  resets cell flags (`syncNewStates`) but not `flagLog`. A builder made before the write and committed after that tree
  is rebased onto it and gets the stale flag back; its pass has `editInfo`, so `syncNewStates` does not run. Contradicts
  reconcile.md §6 ("a tree passed in resets flags written after it was taken"). Failure (`bugs.ts`): `build()`, hide X,
  `setSnapshot(snapshot with X visible)`, commit the builder → `cell.isHidden` undefined (shown), tree and
  `getSnapshot()` hide X. Same with `undo()` in place of `setSnapshot`. In the plugin: a builder outstanding across an
  `await` while the user hides something and a session load, snapshot frame or undo is queued ahead of its commit; the
  saved session and later undo hide what is shown. Fix: on a non-builder tree in `updateTreeAndCreateCtx`, record
  `flagResetAt = flagWrites` (or clear `flagLog`), and skip writes with `n <= max(mark, flagResetAt)`. Add both
  repros to `cell-state.test.ts`; a line in tree.md §4.
- **B6**: `doError` sets `obj`/`cache` undefined for the failed node and recursively its children (also below a node
  turned Null, via `updateSubtree`), with no `dispose`; a later delete calls `dispose` with `b` and `cache` undefined.
  Failure (`bugs.ts`): parent update throws, parent deleted: 2 `dispose` calls, 0 with an object. Skipped bodies:
  StructureRepresentation3D theme/custom-property release, StructureComponent descriptors, Download asset; refcounts
  grow per failed or Null re-apply. Fix: dispose before clearing when `obj` is set and not Null, top first; idempotent
  with the plugin's destroy on `object.removed`.
- **B9**: `subtree` maps refs to `s.cells.get(x)!`, `children` pushes `s.cells.get(c)!`; nodes below a Null or error
  have no cell. Failure (`bugs.ts`): two leaves added below a Null node, `byRef(null).subtree().withTag('x')` throws
  `Cannot read properties of undefined (reading 'transform')`. Callers: focus representation (`:205`), superposition,
  MVS util. Fix: skip refs with no cell; replace `filter`'s `flatMap(p ? [n] : [])` with a loop.
- **B12**: `updateTagged` runs `setParams(...) || setTags(...)`; when params differ, `setTags` never runs. Failure
  (`bugs.ts`): `applyOrUpdateTagged('a', Leaf, {value:1})`, then `applyOrUpdateTagged('a', Leaf, {value:2}, {tags:['extra']})`
  → value 2, tags `['a']`. Callers: `builder/structure.ts:181` `tryCreateComponent`, `representation.ts:213`,
  volume-streaming, anvil, assembly-symmetry. Fix: evaluate both; record a tags action for rebase (B3).
- **B5**: `setParamsHashVersion` overwrites the version `withTags`/`withDependsOn` made with a params hash equal to the
  cell's (`bugs.ts`: tree tags `['new']`, cell and obj tags undefined). Fix: hash `{ params, tags, dependsOn }`, mirrored
  in the MVS virtual builder.
- **B10**: the visitor's `false` only skips the subtree (`bugs.ts`: returns `second`, expected `first`); every miss
  walks the whole tree. Fix: a stoppable `doPreOrder`, return the first match.
- **B7**: `execute` throws `Task.Aborted` after `task.f` returned when an abort was requested; the reconciler's
  `doError` then drops the new object (GPU textures, a Download asset) unstored and undisposed. Fix: return the value;
  the abort lands at the next checkpoint (`revertIfAborted` then sees none on the last node).
- **B8**: `start` resets only after an update, `currentTime`/`globalTime` sum cumulative deltas, so chunks grow to a
  whole window (v1 `verify/v-b8.ts`: 527-578 ms gaps at a 250 ms update rate vs 270-280 ms fixed). Fix: per-chunk
  timing, size by throughput.
- **B3**: `buildTree`'s replay knows add/update/delete/insert; `tag`, `dependsOn`, `updateState` record
  `{kind:'update', params: unchanged}`, a no-op `setParams` on replay. In-repo exposure: B12's path on a rebased
  builder. Fix: record and replay `tags`/`dependsOn`/`state` actions.
- **B4**: `setDependsOn` replaces the transform without `mutateDependency` (which `add` and `fromJSON` maintain), so
  deleting the dependency errors the dependent (`Unresolved dependency.`) instead of cascading. Fix: diff old/new
  `dependsOn`; `add` replacing a ref should also drop old entries.
- **B2**: `updateState` counts as an edit with `lastUpdate`, so it takes the fast track and runs `update()`; cells are
  synced only by `syncNewStates`, skipped for builders (`bugs.ts`: cell undefined, tree true, 1 `update()`). No in-repo
  caller. Fix: a `state` action that schedules nothing and syncs `cell.state` on commit.

## Checked and dropped

- P7's per-transformer cache as a stale-registry regression: not a bug by spec. dependencies.md §3 requires
  `params(undefined, undefined)` to read nothing that changes; transformers whose schema reads the plugin context throw
  and take the structural scan.
- P15 as an MVS duplicate: kept. The MVS research covers building frames; it lists applying them (`setSnapshot`) as
  not measured, and the built-in snapshot animation and scrubbing use the same path.
- B1 re-check: fixed. Hide A, abandoned `build()`, hide B → tree and snapshot hide B; a builder made before a toggle
  commits it; `apply`, `getTree()`, `apply` no longer throws (`verify/v-b1.ts`).
- A transaction around a fast-track edit as a whole-tree path: refuted, 0.02 ms at 10k.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/1f4e2fbb-da6c-4948-aacb-385a4dab465c/scratchpad/perf-audit-2/core-state-task/`.
Bundle with `node build.mjs x.ts x.mjs && node x.mjs [args]` (from `verify/`: `../build.mjs`).

| Script | Settles |
| --- | --- |
| `bench-commit.ts [N]` | P1, P2, P15 (commit paths, walks, `fromJSON`, `setSnapshot`, undo path) |
| `multi-remove.ts` | P1 (N per-structure commits with 1 vs 2 removals) |
| `snap-frames.ts N F [uuid\|rnd15\|rnd22]` (+ `prof-self.mjs` on `prof-snap/*.cpuprofile`) | P15, P3 end to end |
| `hash-coll.ts [N]` | P3 (31-hash collisions, Immutable `get` by key kind) |
| `flaglog.ts [N]` | P17 |
| `produce-bench.ts` | P16 |
| `schemas.ts`, `defmem.ts` | P7 follow-up (structural scan cost), retained param definitions |
| `bugs.ts` | B11 (setSnapshot, undo), B12, B2, B6, B10, B5, B9 |
| `verify/v-b1.ts` | B1 fixed |
