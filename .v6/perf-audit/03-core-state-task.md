# Slice 03: core state and task

Scope: `packages/core/src/state` (the tree, the transient tree, the builder, the reconciler, selections) and
`packages/core/src/task` (the observable runtime, the scheduler, `chunkedSubtask`). The slice is not in the profiled
load pipeline ([00-profile.md](00-profile.md)). The reconciler is linear, and its per-node overhead is small next to
transformer work. The costs only add up on large trees (mesoscale sessions, thousands of entities) and in the plugin's
whole-tree walks on every `changed`. Most of the value is in bugs: B1 loses visibility flags from the tree that gets
saved. Measured on `90d717f49` with Node 25.6 on an M4, using trivial transformers. Nothing here touches frozen gl or
Canvas3D code. The MVS virtual builder is not audited; see [.v6/mvs-optimizations/research](../mvs-optimizations/research/README.md).

## Hot paths

| Path | Cost | Notes |
| --- | --- | --- |
| `State.update`, one-edit fast track | 0.02-0.08 ms at any size | roots = `[lastUpdate]`, no walk |
| `State.update`, anything else (0 edits, 2+ edits, a delete, a tree passed in) | 0.5-0.6 µs/node: 1k 0.45 ms, 10k 5.2-5.8 ms, 20k 12-14 ms, 50k 30 ms | `findDeletes` + `findUpdateRoots`, plus `syncNewStates` for a tree |
| Full evaluation | 3.4-4.6 µs/node of reconciler overhead (20k about 75 ms), plus 1.3 µs/node in the builder | three status events, `params(a)`, log entries, Immutable lookups |
| `fromJSON` / `setSnapshot` with an unchanged tree | 15 / 16 ms at 10k; 36 / 28 ms at 20k | |
| `StateSelection` generators and `findTagInSubtree` | 0.25-0.34 µs/node: 2.5-3.4 ms at 10k | a scan of `state.cells` takes 0.09 ms |
| Immutable `Map.get` with 22-char refs | 140-200 ns, against 60-110 ns with 16-char refs and 6-15 ns for a native `Map` | refs are longer than `STRING_HASH_CACHE_MIN_STRLEN` (16), so every lookup goes through a 255-entry cache that thrashes |
| Retained memory | about 1.5 KB per node, without the `StateObject` | `state-mem.ts` |

Measured and cheap, so not findings: per-node parameter resolution (`params(a)`, `normalizeParams`, `resolveRefs`,
under 0.1 ms even for StructureRepresentation3D on 186k atoms), dependency linking, `DependencyWaiters`, deferral, and
the task runtime's progress overhead (callers check `shouldUpdate` first).

## Performance findings (confirmed and plausible)

| Id | Title | File | Kind | Impact | Estimate | Effort | Verdict | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| P2 | Whole-tree walks pay about 2 Immutable lookups per node; the plugin runs about 10 per `changed` | [immutable.ts:134](../../packages/core/src/state/tree/immutable.ts#L134) | time | medium | 2.7 ms per walk at 10k nodes; about 25 ms per change at 10k (measurement 6, hierarchies 3); 20-50x less with an index or `forEach` | L | confirmed | no |
| P1 | Every commit except a one-edit builder diffs both trees, even with 0 edits | [state.ts:587](../../packages/core/src/state/state.ts#L587) | time | medium | 5.2-5.8 ms per 0- or 2-edit commit at 10k, 12-14 ms at 20k; under 0.5 ms below 1k nodes | M | confirmed | no |
| P3 | 22-char refs thrash Immutable's string-hash cache | [transform.ts:114](../../packages/core/src/state/transform.ts#L114) | time | low | 20k nodes: 0-edit commit 12.7 → 9 ms (−30 %), `fromJSON` 36 → 22 ms (−40 %), `setSnapshot` −10 % | S | confirmed | no |
| P4 | `collectScheduled` repeats the `initCellStatus` walk | [state.ts:1235](../../packages/core/src/state/state.ts#L1235) | time | low | one walk of the evaluated subtrees, about 0.25 µs/node (5-7 % of reconciler overhead) | S | confirmed | no |
| P7 | `tryGetDefaultParamDefinition` throws once per transform object for plugin transformers | [transform.ts:189](../../packages/core/src/state/transform.ts#L189) | time | low | 3.3-4.5 µs per caught throw: about 40 ms per snapshot load at 10k representations, 0.4 ms at 100 | S | confirmed | no |
| P5 | No batched `updateCellState`: hiding a subtree runs one `syncVisibility` per representation | [state.ts:178](../../packages/core/src/state/state.ts#L178) | time | low | k reprs × O(R renderables) hash; about 5-10 ms for k=1000, R=5000 | M | plausible | no |
| P6 | `AsyncQueue` broadcasts every dequeue to all waiters: O(k²) | [async-queue.ts:24](../../packages/core/src/util/async-queue.ts#L24) | time | low | about 500k callbacks for 1000 unawaited commits | S | confirmed | no |
| P9 | Builder `updateState` scans all actions with `find`: O(n²) | [builder.ts:358](../../packages/core/src/state/state/builder.ts#L358) | time | low | latent: no in-repo caller of the real builder's `updateState` | S | confirmed | no |
| P8 | The browser `setImmediate` polyfill uses `window.postMessage('*')` before `MessageChannel` | [scheduler.ts:160](../../packages/core/src/task/util/scheduler.ts#L160) | time | low | one stray `message` event per yield or progress notify, seen by every listener on the page | S | confirmed | no |

### P2: whole-tree walks

- Mechanism: `_doPreOrder` and `_doPostOrder` do `children.get(ref)` and then `transforms.get(child)` for each node.
  Callers that walk the whole tree: `findDeletes`, `findUpdateRoots`, `syncNewStates`, the `StateSelection` generators
  (`ofType`, `ofTransformer`, `rootsOfType`, `subtree`), `findTagInSubtree` and `findWithAllTags`. On each
  `state.data` `changed`, the plugin runs: `StructureMeasurementManager.sync`, whose 6 `getTransforms` calls each run
  `findTagInSubtree` from the root ([measurement.ts:429](../../packages/plugin/core/src/state/manager/structure/measurement.ts#L429));
  and the structure, volume and particle hierarchy managers, each a full pre-order walk with its own copy of the
  visitor (`manager/*/hierarchy-state.ts`).
- Evidence (`walk-bench.ts`, 10k nodes): `doPreOrder` 2.7 ms; building a pre-order array 3.15 ms; scanning it
  0.05 ms; Immutable `transforms.forEach` 0.13 ms (unordered). `sel-bench.ts`: `subtree().withTag()` 3.4 ms,
  `ofTransformer` 2.7 ms, `findTagInSubtree` 2.6 ms, a native scan of `cells` 0.09 ms.
- Fix: (a) Cheap: run order-free searches with `transforms.forEach` (`findTagInSubtree` from the root when the tag is
  unique, `ofType`/`ofTransformer` from the root, then sort the few hits by pre-order if needed). (b) Full: lazily
  memoize on each immutable `StateTree` a pre-order array plus a subtree-end index (an Euler tour). `doPreOrder`
  iterates the array, and a visitor that returns `false` jumps to the subtree end. Building the index costs about one
  walk, so it pays off from the second walk of the same tree, and the plugin does about 10 per change. The live
  `TransientTree` delegates to its source tree when only state flags changed. Move the hierarchy managers onto
  `StateTree.doPreOrder` so they use the index.
- Risk: memory is O(n) per retained tree (history keeps 5), so build the index lazily. If `syncNewStates` switches to
  `forEach`, its event order changes (reconcile.md §7 does not fix that order).

### P1: commits that are not one edit diff both trees

- Mechanism: the fast track needs `count === 1 && lastUpdate && sourceTree === oldTree`. Everything else runs
  `findDeletes` (post-order over the old tree) and `findUpdateRoots` (pre-order over the new tree): a builder with 0
  edits or with 2+ edits, a single delete (`lastUpdate` undefined), or two updates of one node. Only one caller
  guards against a 0-edit commit (`measurement.ts:119`).
- Evidence (`state-bench.ts`): at 10k nodes a 0-edit commit takes 5.2-5.8 ms and a 2-edit commit 5.7-7.5 ms; at 20k,
  12-14 ms; at 1k, 0.45 ms. The fast track takes 0.03-0.08 ms. In the profile, `findDeletes` is 15 % of total and
  `findUpdateRoots` 12.6 %.
- Fix: (1) when `count === 0` and `sourceTree === state.tree`, return from `updateTree` before
  `updateTreeAndCreateCtx`. This also stops a 0-edit commit from overwriting state flags written after `build()` (B1).
  (2) When the builder recorded no delete and did not rebase, use `deletes = []`. (3) A k-edit fast track for an
  un-rebased builder: roots are the topmost touched refs in action order, and deletes are the refs `remove` returned
  that have cells. Keep the whole-tree diff for trees passed in and for rebased builders.
- Risk: the fast track evaluates below Null and error cells and never repairs cells left behind by an update that
  threw (reconcile.md §2, Decisions). A k-edit fast track inherits that, so update the spec and extend the differential
  test (testing.md).

### P3-P9 (low)

- **P3**: Immutable 5.1.9 hashes strings longer than 16 characters through `cachedHashString`, a 255-entry object that
  is reset when full. Generate refs of at most 16 characters, for example a `UUID.create16`, in `StateTransform.create`.
  Versions are not map keys and can stay 22 characters. Refs also appear in snapshots and URLs: grep for code that
  assumes 22 characters. P2 removes most lookups anyway.
- **P4**: `initCellStatus` and `collectScheduled` walk the same `roots` (dependents included). Fill `ctx.scheduled`
  in `initCellStatusVisitor`. Optionally format the `Created X in Y ms` log entry lazily.
- **P7**: `StructureRepresentation3D.params` reads `ctx.representation`, so `params(undefined, undefined)` throws. The
  cache lives on the transform object, so every `withParams` and every `fromJSON` throws again, and the fallback
  `collectStructuralRefs` scans the params. Cache per transformer (a `WeakMap`), or pass the plugin as `globalContext`.
- **P5**: `setSubtreeVisibility` calls `updateCellState` per node, and `UpdateRepresentationVisibility` calls
  `canvas3d.syncVisibility()` per representation. That call is an O(renderables) hash plus dirty flags, not a scene
  rebuild, and Prism's per-event handler is O(1). Add `State.updateCellStates(refs, update)` with a batch boundary, and
  coalesce `syncVisibility` in the behavior. Not measured in a plugin.
- **P6**: keep a `Map<T, resolve>`. `handled` resolves the new head and `remove` resolves the removed entry.
- **P9**: a `Set` of refs that already have an update action, or the separate `state` action from B3.
- **P8**: prefer `MessageChannel` (or `scheduler.postTask`) when it exists, and keep `postMessage` as the fallback.

## Bugs

| Id | Title | File | Severity | Verdict |
| --- | --- | --- | --- | --- |
| B1 | State flags set while a builder is outstanding never reach `state.tree` (frozen maps; the stale builder tree) | [transient.ts:294](../../packages/core/src/state/tree/transient.ts#L294) | high | confirmed |
| B6 | `doError` clears objects and caches without calling transformer `dispose` | [state.ts:1207](../../packages/core/src/state/state.ts#L1207) | medium | confirmed |
| B9 | `subtree()` and `children()` return undefined cells, so the chained filters throw `TypeError` | [selection.ts:272](../../packages/core/src/state/state/selection.ts#L272) | low | confirmed |
| B7 | A late abort discards a completed task's result; nobody disposes it | [observable.ts:106](../../packages/core/src/task/execution/observable.ts#L106) | low | confirmed |
| B8 | `chunkedSubtask` sizing compounds: gaps between yields grow to about 2x the update rate | [chunked.ts:35](../../packages/core/src/task/util/chunked.ts#L35) | low | confirmed |
| B5 | Hash versions cover params only: tag and `dependsOn` edits never reach cells | [transform.ts:277](../../packages/core/src/state/transform.ts#L277) | low | confirmed |
| B10 | `findTagInSubtree` never stops early and returns the last match | [selection.ts:440](../../packages/core/src/state/state/selection.ts#L440) | low | confirmed |
| B3 | Builder rebase replays `tag`, `dependsOn` and `updateState` as `setParams`, which drops them | [builder.ts:86](../../packages/core/src/state/state/builder.ts#L86) | low | confirmed (mostly latent) |
| B2 | Builder `updateState` never reaches `cell.state`, and as a single edit it re-evaluates the subtree | [builder.ts:352](../../packages/core/src/state/state/builder.ts#L352) | low | confirmed (latent) |
| B4 | `setDependsOn` does not update `tree.dependencies`, so deleting the dependency does not cascade | [transient.ts:208](../../packages/core/src/state/tree/transient.ts#L208) | low | confirmed (latent) |

- **B1** (high). There are two ways to lose a flag. (a) `asImmutable()` calls `transforms.asImmutable()`, which in
  Immutable 5 clears the owner of the same map object. The `TransientTree` keeps `changedNodes = true`, so the next
  `assignState` writes into a frozen map and throws the result away. `state.build()`, `transaction()` and a `canUndo`
  update all call it on the live `_tree`. (b) Even without the freeze, a builder's tree is a snapshot taken when
  `build()` ran. `state.tree` keeps its identity, so the builder does not rebase, and the commit installs a tree without
  the flags written in between. Failure (`bug-transient.ts`, `verify/v-b1-b9.ts`): hide A; `build()` (never
  committed, for example the `count === 0` return at `measurement.ts:119`, or a commit waiting in the queue during a
  long load); hide B. `cell(B).state.isHidden` is true, but in the tree and in `getSnapshot()` it is undefined. The
  saved session shows B, and a later `setSnapshot` or undo (`syncNewStates`) shows it again. Also
  (`bug-builder.ts`): `apply`, `getTree()`, `apply` throws `Could not find node`. Fix: in `asImmutable`, adopt the
  result as the new base and reset `changed*`, `_childMutations`, `_dependencyMutations` and `_stateUpdates`. The
  last also ends the in-place `assignState` on transforms that snapshots share (tree.md §3, an open item). For (b),
  carry `cell.state` into the committed tree for builder commits, or give `updateCellState` a tree identity the
  builder detects. Update tree.md §2 and §4.
- **B6** (medium): a node that throws, or the children of a node that turned Null (`updateSubtree` calls
  `doError(child, undefined, true)`), loses `obj` and `cache` without `dispose`. Later deletes call `dispose` with
  `b` and `cache` undefined, so nothing is released. Skipped dispose bodies: StructureRepresentation3D
  (`ensureCustomProperties.detach`, `Theme.releaseDependencies`), StructureComponent
  (`customPropertyDescriptors.dispose()`), Download (`asset.dispose()`). Failure: a structure update returns Null or
  throws. Every component and representation below it keeps its custom-property and theme refcounts, and they grow on
  every re-apply. Fix: call `dispose` before clearing when `obj` is set and not Null, top of the subtree first, as the
  deletes path does. It must be idempotent with the plugin's destroy on `object.removed`.
- **B9**: invariant 1 allows tree nodes without cells (below a Null object). Failure (`verify/v-b9.ts`): add a child
  below a Null node in a 2-edit builder; `byRef(null).subtree().withTag('x')` throws `Cannot read properties of
  undefined (reading 'transform')`, and `children().filter` throws on `status`. Callers include
  structure-focus-representation, superposition and MVS `util.ts`. Fix: skip refs with no cell. Also replace
  `filter`'s `flatMap(p ? [n] : [])`, which allocates an array per cell, with a loop.
- **B7**: `execute` throws `Task.Aborted` after `task.f` has returned when an abort was requested. In the reconciler
  this becomes `doError`: the new object (for example texture-mesh GPU textures from `createOrUpdate`, or a Download
  asset in the cell cache, which B6 also drops) is never stored, disposed or announced. Fix: return the value and let
  the abort land at the next checkpoint. An abort on the last node is then ignored, so `revertIfAborted` sees no abort.
- **B8**: `start` is reset only after an update, so each `delta` is the time since the last update, and
  `currentTime`/`globalTime` sum cumulative values. After an update, `chunkSize` comes out as the number of items in
  the last window, so one chunk lasts a whole window, and a chunk longer than the window is a fixed point. Yields
  happen only in `update`, which runs only when `shouldUpdate`. Evidence (`verify/v-b8.ts`, fake clock, 250 ms rate):
  the maximum gap between yields is 527 ms (CIF-like 300M tokens) and 578 ms (60 s of 1 ms items). Per-chunk timing
  with throughput sizing keeps it at 270-280 ms. Correction: the fix bounds blocking at about the update rate plus one
  chunk, not 30-50 ms. Fewer blocking milliseconds would need yields between chunks or a lower update rate.
- **B5**: `setParamsHashVersion` overwrites the new UUID version that `withTags`/`withDependsOn` made with a params
  hash equal to the cell's. Failure (`bug-hashver.ts`): after `tag('new-tag')` and `dependsOn('x')` through
  `getTree({ useHashVersion: true })`, the tree has both, but `cell.transform.tags`, `obj.tags` and the dependency
  links do not. The MVS `VirtualStateBuilder` mirrors this hash (`load-generic.ts:757`, load.md §3), so both change
  together: hash `{ params, tags, dependsOn }`.
- **B10**: `false` only skips the subtree, so later matches overwrite `s.ref`. Make `doPreOrder` stoppable and return
  the first match. It does not help the measurement manager (see Checked and dropped).
- **B3**: the replay calls `setParams(ref, sameParams)`, which is a no-op. Failure (`bug-rebase.ts`): `tag`,
  `dependsOn` and `updateState({ isCollapsed })`, another commit lands, then commit: all three are undefined. In-repo
  exposure: only `applyOrUpdateTagged` with extra `options.tags`. `tag`, `dependsOn` and `updateState` on the real
  builder have no callers outside tests (MVS uses its virtual builder). Fix: record `tags`, `dependsOn` and `state`
  actions and replay them with `setTags`, `setDependsOn` and `assignState`.
- **B2**: the edit counts toward the fast track, so the node and its 4-node subtree run `update()` 4 times, and
  `cell.state` is never synced: `syncNewStates` runs only without a builder (`bug-updatestate.ts`). This contradicts
  tree.md Decisions ("must not evaluate it"). Fix: a `state` action that schedules nothing and syncs `cell.state` on
  commit.
- **B4**: `setDependsOn` skips `mutateDependency`, which `add` and `fromJSON` both call (`bug-dependson.ts`: the
  dependent errors with `Unresolved dependency.` instead of being removed; after a snapshot round trip it cascades).
  Fix: diff the old and new `dependsOn`. `add`, when it replaces an existing ref, also never removes old entries.

## Checked and dropped

- B10 as a speed-up for `StructureMeasurementManager.sync`: refuted. The measurement group is a child of the root
  created after the structures, so it comes late in pre-order, and before the first measurement the tag does not exist
  and every call walks the whole tree. Use `forEach`, the root's children, or a tag index instead (P2a).
- B8's "blocks at most about 30-50 ms instead of 270-530 ms": refuted. Yields happen only on update, so the fixed loop
  still blocks for about the update rate (270-280 ms in the sim).
- P5's "Prism's scene source reacts per event" as a cost: refuted. The handler sets one signal. `syncVisibility` is a
  hash over renderables with lazy dirty flags, not O(scene) work.
- B2 and B3 "happen routinely while a long update is running": refuted for this repo. There are no callers of the real
  builder's `updateState`, `tag` or `dependsOn` outside tests. Kept as latent API bugs.
- B1's fix in `asImmutable` alone: incomplete. A builder made before a toggle still commits a tree without the flag
  (`verify/v-b1-b9.ts`, case B1b).
- P2's "the index never costs more than a walk": building it is 3.15 ms against 2.7 ms for a walk, so it pays off
  only from the second walk of the same tree.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/ac6cef1a-5954-43a5-ae83-c22fef35f065/scratchpad/perf-audit/core-state-task/`.
Bundle with `node ../build.mjs x.ts x.mjs && node x.mjs` (use `build2.mjs` for plugin imports; from `verify/` use
`../../build.mjs`).

| Script | Settles |
| --- | --- |
| `state-bench.ts [N]`, `state-bench-short.ts [N] --short` | P1, P3 (commit, snapshot and `fromJSON` costs; 22-char vs 15-char refs) |
| `state-build.ts` (+ `prof-self.mjs` on its `.cpuprofile`) | full evaluation overhead, P4 |
| `walk-bench.ts`, `sel-bench.ts` | P2, B10 |
| `imm-keys.ts` | P3 (Immutable `get` by key length) |
| `state-mem.ts` | retained bytes per node |
| `params-bench.ts` | parameter resolution (dropped as cheap) |
| `throw-cost.mjs` | P7 |
| `chunked-sim.ts`, `verify/v-b8.ts` | B8 (chunk durations; gaps between yields, current vs fixed) |
| `bug-transient.ts`, `bug-builder.ts`, `verify/v-b1-b9.ts` | B1 (a), the builder after `getTree()`, B1 (b) |
| `bug-updatestate.ts`, `bug-rebase.ts`, `bug-dependson.ts`, `bug-hashver.ts` | B2, B3, B4, B5 |
| `verify/v-b9.ts` | B9 |
