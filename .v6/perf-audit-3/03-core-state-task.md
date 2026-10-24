# Slice 03: core state and task (v3)

Scope: `packages/core/src/state` (tree, transient tree, builder, reconciler, diff, selections), `packages/core/src/task`
(observable runtime, scheduler, `chunkedSubtask`) and `core/src/util/async-queue.ts`. HEAD `16277d822`, v2 base
`0944823d4`. Not on the profiled load path: cost scales with state nodes n, not atoms.

Measured with Node 26.12 (M4) on a shared machine, **load 5 → 41 during the runs**: wall times carry far more than the
usual noise, so the A/B uses allocation (`total_allocated_bytes`) and back-to-back HEAD/base runs of the same bundle
(base = `git archive 0944823d4 packages/core`, bundled with the same esbuild setup). Trivial transformers, 10
`stateUpdated` subscribers, trees of groups with 9 leaves.

PRs since v2 that touched the slice: #2025 (03/P1: k removals from the edits), #2033 (03/B12), #2034 (03/B9, 03/B10),
#2042 (`relinkCells` set for > 8 dependencies), #2045 (03/B8 chunk sizing), #2046 (03/P4, P6, P9, P17, B2 partly, B3,
B4).

## Hot paths now

| Path | Cost | Where the time goes |
| --- | --- | --- |
| Builder commit, edits path (any number of removals, no move, not rebased, cells in sync) | 10k: 1 update 0.04 ms, 2 updates 0.07-0.13 ms, 2 removals 0.018 ms (base 6.3 ms), 20 removals 209 KB alloc (base 4.3 MB) | edits, their ancestors, the removed subtrees |
| `setSnapshot`, one frame changing one node (P15) | 1k: 1.6-1.8 ms, 2.6 MB alloc; 10k: 24-26 ms (load 5), 27.6-28.3 MB alloc | `fromJSON` 21.3 MB of the 28 MB (77 %), then whole-tree diff + `syncNewStates` |
| A `StateTree` passed in (undo, transaction revert) | 10k: 9.4-9.7 ms, 6.6 MB alloc; 1k: 0.7 ms, 0.68 MB | `findDeletes` + `findUpdateRoots` + `syncNewStates` |
| Whole-tree walks (P2): `doPreOrder`, `ofTransformer`, `findTagInSubtree` miss, `getSnapshot` | 10k: 2.7 / 3.0 / 2.5 / 3.1 ms (load 5); a bare `doPreOrder` allocates 2.0 MB (~200 B/node) | Immutable `forEach` per child set + 2 lookups per node |
| Builder commit with a flag written while it was outstanding | 10k flagged refs: 0.09 ms (base 16.9 ms under load; v2 2.6 ms) | `flagWriteOrder` from the mark |
| Multi-file parse with several chunked tasks at once (P18) | 4 concurrent 40M-item tasks: 0.7-1.1 M chunks (base 0.3-38 k) | chunk sizing from wall clock between own updates |

Benches: `bench-commit.log`, `bench-commit-base.log`, `alloc.log` in the scratch folder (see the end).

## v2 → now

Counts (22 v2 items: 11 P, 11 B; 14/D2 = P1): **fixed 12** (P1, P4, P6, P9, P17, B3, B4, B8, B9, B10, B12, 14/D2), **partly
fixed 2** (P16, B2), **still open 9** (P2, P3, P5, P8, P15, B5, B6, B7, B11), **regressed 0**, **dropped 0**. One new
cost from a fix: P18 (#2045, the B8 fix).

| v2 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| P1 | Commits removing 2+ nodes diff both trees | fixed | #2025 | `editsOnly` drops the `removed.length <= 1` gate; 2 removals 6.3 → 0.018 ms at 10k; `multi-remove`: N = 1000 per-structure 2-removal commits 51 s (base, loaded) → 50 ms total. Differential test covers k removals and ancestor pairs |
| P2 | Whole-tree walks, 2 Immutable lookups per node | still open | | Unchanged; new datum: ~200 B allocated per visited node |
| P3 | 22-char counter refs, 31-hash collisions | still open | | `transform.ts` unchanged |
| P4 | `collectScheduled` repeats the `initCellStatus` walk | fixed | #2046 | `ctx.scheduled` filled in `initCellStatusVisitor`; same set (dependents pushed before the walk, as before) |
| P5 | No batched `updateCellState` | still open | | unchanged |
| P6 | `AsyncQueue` broadcasts every dequeue | fixed | #2046 | `Map<T, resolve[]>`; same semantics (head resolved by `handled`, removed entry by `remove`) |
| P8 | `setImmediate` polyfill prefers `postMessage('*')` | still open | | [scheduler.ts:160](../../packages/core/src/task/util/scheduler.ts#L160) unchanged |
| P9 | Builder `updateState` scans actions with `find` | fixed | #2046 | a `state` action instead |
| P15 | Every snapshot frame rebuilds and diffs the whole tree | still open | | re-measured, see below |
| P16 | `update(recipe)` returning a new object walks the params | partly fixed | #2034 | callers in the focus representation pass params; `builder.update` unchanged |
| P17 | `carryFlagWrites` scans every ref ever flagged | fixed | #2046 | `flagWriteOrder` + binary search, compacted to each ref's last write; checked correct (a ref's last write > mark iff any write > mark) |
| B2 | Builder `updateState` never reaches `cell.state` | partly fixed | #2046 | no evaluation, replayed on a rebase; cell still stale, see Bugs |
| B3 | Rebase replays `tag`/`dependsOn`/`updateState` as `setParams` | fixed | #2046 | `tags`/`dependsOn`/`state` actions |
| B4 | `setDependsOn` skips `tree.dependencies` | fixed | #2046 | `relinkDependencies` in `setDependsOn` and in `add` replacing a node; `OrderedSet.add` keeps the second add in `add` idempotent |
| B5 | Hash versions cover params only | still open | | probe: tree tags `['new']`, cell and obj tags undefined |
| B6 | `doError` clears objects without `dispose` | still open | | probe: 2 dispose calls, 0 with an object |
| B7 | Late abort discards a completed result | still open | | [observable.ts:106](../../packages/core/src/task/execution/observable.ts#L106) unchanged |
| B8 | `chunkedSubtask` sizing compounds | fixed | #2045 | one task: max gap 316-346 ms at 250 ms rate (base 488 ms; 1044 ms with a stalling thread vs 325 ms). New cost under concurrency: P18 |
| B9 | `subtree()`/`children()` return undefined cells | fixed | #2034 | probe: no throw |
| B10 | `findTagInSubtree` never stops | fixed | #2034 | probe: returns `first`; recursion stops at the first match |
| B11 | Carried flags override a tree passed in after them | still open | | probe unchanged: after `setSnapshot` (or `undo`) then the builder commit, cell shown, tree and snapshot hidden |
| B12 | `applyOrUpdateTagged` drops tags when params changed | fixed | #2033 | probe: tags `a, extra` |

Fix review, nothing wrong found in: `findEditedDeletes` for k roots (dedup, nested-root drop, pre-order sort by
`childOrder`, post-order runs; equal to the whole-tree diff for removed-then-re-added refs and cascaded dependents),
`initCellStatus` collecting `scheduled`, `AsyncQueue` waiters (abort of the running head is resolved by `handled` in the
`finally`, as before), `flagWriteOrder` compaction, `updateTagged` both setters, `_findTagInSubtree` (stops via
Immutable `forEach` returning false), `relinkDependencies`. Rebase replay of a `tags`/`dependsOn` action on a node
deleted in the newer tree throws a `TypeError` (`transforms.get(ref).tags`) where `setParams` threw "not present":
same outcome, worse message.

## Open performance findings (ranked)

| Rank | Id | Status | Title | File | Kind | Estimate | Effort | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | P15 | still open | Every snapshot frame rebuilds (`fromJSON`) and diffs the whole tree | [state.ts:233](../../packages/core/src/state/state.ts#L233) | time, memory | measured: 24-26 ms and 28 MB alloc per frame at 10k, 1.6-1.8 ms / 2.6 MB at 1k → estimated 0.1-0.3 ms with a frame-delta path | M | no |
| 2 | P2 | still open | Whole-tree walks: Immutable `forEach` + 2 lookups per node, ~200 B alloc per node | [immutable.ts:134](../../packages/core/src/state/tree/immutable.ts#L134) | time, memory | measured: 2.5-3.1 ms and 2.0 MB per 10k walk; estimated 20-50x with a memoized pre-order index | L | no |
| 3 | P18 | new (from #2045) | `chunkedSubtask` chunks collapse when several chunked tasks interleave | [chunked.ts:47](../../packages/core/src/task/util/chunked.ts#L47) | time | measured (sim): 4 concurrent tasks 0.7-1.1 M chunks vs 0.3-38 k base, wall 1.5-2.4x base | S | no |
| 4 | P16 | partly | `builder.update(recipe)` returning a new object: mutative walks it | [builder.ts:367](../../packages/core/src/state/state/builder.ts#L367) | time | v2: 0.07 / 0.63 / 2.25 ms for 1k / 100k / 1M-int bundles; callers fixed | S | no |
| 5 | P3 | still open | Counter refs: 22 chars, 31-hash collisions | [transform.ts:114](../../packages/core/src/state/transform.ts#L114) | time | v2: Immutable get 153-161 vs 50 ns (16 chars); subsumed by P2/P15 | S | no |
| 6 | P5 | still open | No batched `updateCellState` | [state.ts:242](../../packages/core/src/state/state.ts#L242) | time | v2: k `syncVisibility` → 1 for hide-all (plausible, not measured) | M | no |
| 7 | P8 | still open | `setImmediate` polyfill uses `postMessage('*')` first | [scheduler.ts:160](../../packages/core/src/task/util/scheduler.ts#L160) | time | one stray `message` event per yield | S | no |

### P15: snapshot frames

- Re-measured (`bench-commit.ts`, `alloc.ts`, same process, HEAD vs base back to back): one leaf changed per frame
  1.77 / 1.65 ms at 1k, 26.2 / 24.2 ms at 10k (load 5); allocation per frame 2.6 MB (1k), 27.6-28.3 MB (10k), of which
  `StateTree.fromJSON` 2.3 / 21.3 MB. Undo path (a tree passed in, no change) 6.6 MB, 9.4-9.7 ms at 10k. Unchanged since
  v2; P1's fix does not touch it (a tree passed in has no `editInfo`).
- Now the top item of the slice. Mechanism and fix as v2: a frame-delta path in `setSnapshot` (compare the serialized
  transforms by ref, parent, version and flags against the current tree; when the shape is equal, edit a transient of
  the current tree for the changed transforms and commit through the edits path; fall back otherwise). Keep B11's
  semantics in mind: the delta path is a tree passed in, so its flags must win (and reset `flagLog`, see B11).

### P2: whole-tree walks

- New datum (`alloc.ts`): a bare `doPreOrder` over 10k nodes allocates 2.0 MB (~200 B/node), 213 KB at 1k:
  `children.forEach(_preOrderFunc, ctx)` on an Immutable `OrderedSet` allocates iteration closures per child set, plus
  the lookups. `getSnapshot` 3.5 MB at 10k. Fix as v2 (memoized pre-order ref array per immutable tree; order-free
  searches over `transforms`/`cells`).

### P18 (new): chunk sizes collapse with interleaved chunked tasks

- Mechanism: after an update, `chunkSize = processedSinceUpdate / 4`. `shouldUpdate` is wall clock since this
  context's last `update`, which includes the time other tasks ran while this one yielded. When the others hold the
  thread for 250 ms or more, a resumed task's first chunk already triggers an update with `processedSinceUpdate` = one
  chunk, so the chunk shrinks 4x per turn; the others see the same, and sizes fall to 1. The pre-#2045 formula
  (`currentTime × globalProcessed / globalTime`) also collapsed in this case, less often.
- Evidence (`chunk-sim.ts`, `chunk-ab.sh`; real `chunkedSubtask` and `ObservableRuntimeContext`, 250 ms rate; 3
  alternating runs): four concurrent tasks of 40M cheap items: HEAD 0.68-1.09 M chunks, 57-96 updates, wall 5.7-11.5 s;
  base 0.3-38 k chunks, 24-54 updates, 2.4-6.5 s. Most HEAD tasks end at chunk size 1. One task alone: HEAD 346-470
  chunks, 1.2-1.7 s; base 41-96 chunks, 1.5-2.8 s (HEAD better). Wall is load-sensitive; the chunk counts are not.
- Exposure: low in the plugin. Transformers run one at a time through the state's `AsyncQueue`, so concurrent chunked
  parses need tasks outside it (several `plugin.runTask` parses at once, a worker-less app that parses in parallel).
- Fix: size from the processing rate over `f`'s own time (one `now()` before and after each `f`; `shouldUpdate` already
  calls `now()`), e.g. `chunk = rate × interval / 4`, so time spent in other tasks does not count; or shrink only when
  the chunk itself took longer than its target. Keep the B8 test (gaps within ~1.25x) and add a two-task case to
  `chunked.test.ts`.

### Low

- **P16**: core `update(recipe)` still runs `produce` over a returned object; document that a recipe should mutate the
  draft, or detect a returned non-draft and use it as is.
- **P3, P5, P8**: code unchanged; v2 evidence and fixes stand.

## Bugs

| Id | Severity | Status | Title | File | Failure scenario |
| --- | --- | --- | --- | --- | --- |
| B11 | medium | still open | Carried flags override a tree passed in after them | [state.ts:492](../../packages/core/src/state/state.ts#L492) | `build()`, hide X, `setSnapshot` (or `undo`) with X shown, commit the builder → cell shown, tree and `getSnapshot()` hide X |
| B6 | medium | still open | `doError` clears objects and caches without `dispose` | [state.ts:1264](../../packages/core/src/state/state.ts#L1264) | parent update throws, then parent deleted → 2 `dispose` calls, 0 with an object; theme/custom-property refcounts and Download assets below leak |
| B5 | low | still open | Hash versions cover params only | [transform.ts:281](../../packages/core/src/state/transform.ts#L281) | hash-versioned tag-only edit → tree tags `['new']`, cell and obj tags undefined |
| B7 | low | still open | A late abort discards a completed result, undisposed | [observable.ts:106](../../packages/core/src/task/execution/observable.ts#L106) | abort requested after `task.f` returned → `Task.Aborted` thrown, the built object leaks |
| B2 | low | partly fixed (#2046) | Builder `updateState` never reaches `cell.state`; the cell can no longer undo it | [builder.ts:383](../../packages/core/src/state/state/builder.ts#L383) | `build().to(x).updateState({isHidden: true})`, commit → tree hidden, cell shown, 0 `changed` and 0 `stateUpdated` events; then `updateCellState(x, {isHidden: false})` is a no-op (the cell already says shown), so the tree, the session and undo keep X hidden |

- **B11**: unchanged since v2; `updateTreeAndCreateCtx` still resets nothing for a non-builder tree. Fix as v2: on a
  plain tree record `flagResetAt = flagWrites` and carry only writes with `n > max(mark, flagResetAt)` (with
  `flagWriteOrder` now: start the binary search at `max(mark, flagResetAt) + 1`). Add both repros to
  `cell-state.test.ts`.
- **B2** (`b2.ts`): #2046 made `updateState` cost nothing (count 0, no fast track, a `state` action for rebases), but
  cells are still synced only by `syncNewStates`, which a builder pass skips. The new detail: the stale cell turns the
  user's later toggle into a no-op (`updateCellState` compares against `cell.state`), so the tree's flag sticks. No
  in-repo caller of the builder's `updateState` (`grep`), so latent. Fix: on commit, sync `cell.state` from the tree for
  the refs of `state` actions (and emit `stateUpdated`), or route the builder's flags through `updateCellState` after
  the pass.
- **B6, B5, B7**: probes (`bugs.ts`) and code as v2; fixes as v2.

## Already good, re-checked

- `0-edit commit` 0.002-0.003 ms; fast-track commit 0.03-0.04 ms at 1k and 10k; transaction around a fast-track edit
  0.02 ms.
- `findEditedDeletes`' sort cost: one `childOrder` pass per parent where paths part, `pathOf` per removal (depth); 20
  removals at 10k 209 KB alloc.
- `logFlagWrite` compaction: amortized O(log F) per flag write (sort of F entries every ~F + 64 writes).
- `pnpm test:unit packages/core/src/state packages/core/src/task`: 11 files, 66 tests pass.

## Not audited

- `detectDependenciesCycle` runs a DFS per added cell over its dependency closure: O(n²) for a long `dependsOn` chain
  added in one commit. No in-repo chain found; not measured.
- Plugin-side callers of `AsyncQueue` and `updateCellState` (slice 12).

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/31ee827f-4715-4967-bec3-085f5013d11f/scratchpad/pa3/03/`.
`node build.mjs x.ts x.mjs [base]` bundles against HEAD or against `base/` (v2 core).

| Script | Settles |
| --- | --- |
| `bench-commit.ts [N]` | P1, P2, P15, commit paths (logs: `bench-commit.log`, `bench-commit-base.log`) |
| `alloc.ts [N]` | allocation and CPU per commit kind, `fromJSON`, frame, undo, walk (`alloc.log`) |
| `multi-remove.ts` | P1 at N = 100-1000 structures |
| `flaglog.ts [N]` | P17 |
| `snap-frames.ts N F` | P15 (too load-sensitive this run; use `alloc.ts`) |
| `chunk-sim.ts mode items work`, `chunk-ab.sh` | B8 fix, P18 |
| `bugs.ts` | B11, B12, B2, B6, B10, B5, B9 |
| `b2.ts` | B2 remainder |
