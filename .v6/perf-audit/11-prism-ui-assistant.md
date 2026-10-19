# Slice 11: `prism/ui` and `prism/assistant`

The Prism DOM runtime and the assistant UI: the reactive core (`reactive/state.ts` `Listened`/`State`/`derive`,
keyed state and selectors, the frame lane, work queues), the DOM runtime (`dom/scope.ts` bindings and flush,
`dom/render.ts` insert/regions/delegated events, `dom/control-flow.ts` `Show`/`For`), virtual lists, the undo
`HistoryManager`, `CommandManager`, and in the assistant the turn card (`viewer/agent-model.ts`, `views/agent-turn.tsx`)
and the Inspector (`inspector/model.ts`, `inspector/trace.tsx`). It is outside the profiled data pipeline: costs scale
with UI items (rows, commands, trace steps, history steps), except the history checkpoint, which scales with the scene.
Measured on `v6-opts-261019` (`90d717f49`), Node 25.6 + jsdom, Apple M4. Each finding was re-read with its callers by a
second pass; disputed impacts were re-measured. Nothing here touches the frozen Canvas3D/gl code or MVS.

## Hot paths

| Path                                                        | Runs                                                       | Cost driver                         |
| ----------------------------------------------------------- | ---------------------------------------------------------- | ----------------------------------- |
| `HistoryManager.run` → `#take` (checkpoint + `size`)        | every non-read API op: each form edit, each agent tool call | snapshot JSON size (~4 ms/MB)       |
| `For` reconcile                                             | every list change                                          | O(n) at ~0.1 µs/row (8000 rows 0.84 ms) |
| `Show`/`Switch`/`Dynamic` regions                           | per mounted instance                                       | rxjs chain memory                   |
| Delegated dispatch (`dispatchDelegated`)                    | every bubbling event of 21 types, per delegation root      | `composedPath()` per root           |
| Record-state sets (`sameShallow` → `shallowEqual`)          | turn card per token, list windows per scroll frame         | 4 arrays per compare                |
| Inspector overview and trace list                           | every 80 ms settle while a turn runs (diagnostic page)     | running turn's step count, session rows |
| `CommandManager.handleKey`                                  | every keydown                                              | O(registered commands)              |

Already fine: listener slots with O(1) removal; counting-sort flush with reused arrays; allocation-free frame lane and
work queues; keyed per-row notification (`KeyedState.replace`, `KeyedSelector`); `For`'s unchanged-list fast path,
LIS placement and time-sliced creation; `RecycledList` slot pooling; the streamed reply is one Text node; the Inspector
keeps rows per turn and per step (`turnRowsOf`); `derive` is lazy while unlistened. No quadratic path in the primitives.

## Performance findings

| Id  | Title                                                                                     | Where                                                                                                         | Kind        | Impact | Estimate                                                                                   | Effort | Verdict   | Frozen |
| --- | ----------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- | ----------- | ------ | ------------------------------------------------------------------------------------------ | ------ | --------- | ------ |
| P1  | History: full snapshot + JSON size on every run, also in a turn and past the hard limit   | [history.ts:608](../../packages/prism/ui/src/managers/history.ts#L608)                                        | time        | medium | −~4 ms/MB per API op; 20–200 ms per agent call at 5–50 MB; ≥ 0.5 s per call past the hard limit | S      | confirmed | no     |
| P3  | `Show`/`Switch`/`Dynamic` built on rxjs chains                                            | [control-flow.ts:38](../../packages/prism/ui/src/dom/control-flow.ts#L38)                                     | memory      | low    | −1.3 KB per instance (4.1 → 2.8 KB); ~0.3–3 MB at realistic mounted counts                 | M      | confirmed | no     |
| P2  | Inspector rebuilds the running turn's column every settle; trace list not virtualized     | [trace.tsx:151](../../packages/prism/assistant/src/inspector/trace.tsx#L151)                                  | time+memory | low    | O(S) blocks per settle → O(new steps); trace DOM bounded to a window                        | M      | confirmed | no     |
| P4  | `handleKey` maps every command to an object per keydown                                   | [manager.ts:740](../../packages/prism/ui/src/commands/manager.ts#L740)                                        | time+memory | low    | 17 µs (300 cmds) / 190 µs (3000) per key → O(commands on that chord), no allocation        | S      | confirmed | no     |
| P5  | `shallowEqual` allocates four arrays per compare                                          | [shallow-equal.ts:9](../../packages/prism/ui/src/reactive/shallow-equal.ts#L9)                                | time+memory | low    | 920 → 46 ns (12 keys), 228 → 15 ns ({first,last}); zero allocation                          | S      | confirmed | no     |
| P6  | Delegated dispatch calls `composedPath()` per root for every pointermove                  | [render.ts:293](../../packages/prism/ui/src/dom/render.ts#L293)                                               | time+memory | low    | no allocation per pointer move while no handler of that type is mounted                    | S      | confirmed | no     |
| P9  | `TreeNavigation.#indexOf` scans all rows per key                                          | [tree-keys.ts:118](../../packages/prism/ui/src/dom/tree-keys.ts#L118)                                         | time        | low    | O(rows) → O(1) per arrow key; µs today (state-tree cells only)                             | S      | confirmed | no     |

### P1: history checkpoint per run

- **Mechanism**: every non-read API op goes through `history.run` ([run.ts:262](../../packages/prism/viewer/src/api/run.ts#L262));
  an undoable run calls `#take(!!group)`: `port.checkpoint()` (`plugin.state.getSnapshot`), then `port.size()`
  = `checkpointSize` ([checkpoint.ts:81](../../packages/prism/viewer/src/sources/checkpoint.ts#L81)), a second full
  `JSON.stringify` with a replacer. `ms` is taken before `size`, so the time budget never sees it. A run that joins an
  open turn (`beginGroup`, [session.ts:197](../../packages/prism/viewer/src/agent/session.ts#L197)) passes
  `group=true`, which bypasses the checkpointing-off flag: every tool call still snapshots and measures. Past the hard
  limit nothing is remembered: each call of each later turn takes the snapshot, measures it, and discards it ("too
  large to undo"). The turn's per-run checkpoint is rollback-only; `#record` appends to the open entry and drops it.
- **Evidence**: `bench-ckpt.mjs` (snapshot-shaped, cartoon params with a 30-color list): 500 transforms 0.35 MB 1.4 ms,
  2000 1.41 MB 5.6 ms, 8000 5.65 MB 22.1 ms (~4 ms/MB, size alone; `getSnapshot` comes on top). Typical scenes
  (tens of transforms, < 0.5 MB) pay < 2 ms; the cost matters on large sessions and in agent turns.
- **Fix**: (1) Skip `size()` for runs that join an open group and for nested runs (budget decided at the scene's last
  measured checkpoint). (2) Reuse the last `bytes` unless the data tree's transform count or version changed since.
  (3) Count `size()` in `ms`. (4) Remember "over the hard limit" until the scene shrinks (`session.clear`, or a
  measured drop) and return `unsafe` at once instead of snapshotting per call. (5) Optional: a total-bytes cap over
  `#entries` + `#redo`, evicting the oldest checkpoints.
- **Risk**: a fast-growing scene stays over budget for a few steps before checkpoints turn off; a run in a turn past
  the hard limit has no rollback either way (same as today). Tests: history turn-off (budget, hard limit,
  `SlowCheckpointsToTurnOff`), rollback of a failing run in an open group, `undoTurn`. `sequence.md §7.4` says "each
  checkpoint is measured": update it with (1)–(2).

### P3: rxjs inside `Show`

- **Mechanism**: `Show` wraps `toObservable(when).pipe(distinctUntilChanged()).subscribe` in a new `Observable`, handed
  to `insertRegion`; each instance keeps the Observable, the operator Observable, Subscriber chains and Subscriptions on
  top of the region's anchors and scope. `Switch` adds `combineLatest`, `Dynamic` a `map`.
- **Evidence**: `bench-show2.mjs`, 4000 instances, half shown: `Show` 4061–4204 B retained per instance, the same
  region from `derive` + `insertRegion` 2795–2809 B; make 42–51 vs 38–41 ms, toggle-all ~22 ms both. Mounted counts
  are bounded by the visible UI (virtual lists keep a window), so realistic totals are hundreds to a few thousand
  instances: 0.3–3 MB, not the 6.5 MB per 5k rows the auditor projected.
- **Fix**: when `when` is a state (`isState`), give `insertRegion` a state built by `derive(when, …)` with the memoized
  render function and an `Object.is` check; keep rxjs only for foreign Observables. Same for `Switch` (a multi-source
  `derive`) and `Dynamic`.
- **Risk**: the pause/resume (`fresh`/`resumed` with `sameShallow`) and the closing animation must be reproduced.
  Tests: control-flow, `hidden-quiet.test.tsx`, `leaks.heap.tsx`, the waste budget.

### Low

- **P2**: `model.turns` changes on each settle while a turn runs (its `TurnView` is new), so `For` (identity keyed,
  see B1) disposes and rebuilds the running turn's `Column` and all its `StepBlock`s (button + watch + selector key)
  every 80 ms; O(S²) block builds per turn. The trace list mounts every row of the session (1.6k rows for 50 × 30) and
  re-reconciles O(rows) per settle (~0.2 ms at 2000, `bench-for.mjs`). At typical S (5–30) it is < 0.3 ms per settle;
  the Inspector is a diagnostic page. Fix: inner `<For by={s => s.id}>` over the turn's steps (StepViews are kept by
  the TraceCache), flex-grow from a derived sum; `RecycledList` for the rows (sticky turn head as a slot kind). Risk:
  `reveal`/`scrollToIndex` and j/k with virtualized rows; tests `views.test.tsx`, `rows.scaling.ts`.
- **P4**: `handleKey` does `map` (one object per command) + `filter` + `sort` before the `typing` check. Shape bench
  (`bench-handlekey.mjs`): 17 µs per keydown at 300 commands, 190 µs at 3000; lazy commands are not made there, so C is
  the eager set. Fix: a `Map<chordId, CommandDefinition[]>` invalidated on `#commands` or `keymap.layers` change,
  looked up for the event's chord (both mac variants and the `code` fallback of `chordMatches`); bail early when typing
  and the chord types text. Risk: chord normalization must equal `chordMatches`.
- **P5**: `shallowEqual` builds `Reflect.ownKeys` × 2 and filters both. `bench-shallow.mjs`: 920 → 46 ns (12 keys),
  228 → 15 ns ({first,last}); `RecordState.set` 0.265 → 0.029 µs. Absolute cost is small (per token < 1 µs), but it is
  the default compare of every record state and breaks the zero-allocation rule. Fix: a `for…in` + `Object.hasOwn`
  count loop; symbol keys only if `getOwnPropertySymbols` is non-empty.
- **P6**: `document` plus every `render()` container and `MspElement` is a root; each calls `composedPath()` and two
  `indexOf` per bubbling event, and writes `Walked` on the event. Pointer moves over the canvas reach 2–3 roots at
  60–240 Hz with no handler anywhere (only the command menu's 4 `onPointerMove` and 2 `onPointerOver` exist, mounted
  while open). Fix: a live per-type handler count (inc in `setHandler`, dec when the element's scope disposes), return
  at once at 0; a "types ever registered" flag is not enough (the command menu sets it on first open). Or attach the
  pointermove/mousemove/over/out handlers natively via `on()`. Keep the composed-path walk otherwise (detached-subtree
  and shadow-root behaviour).
- **P9**: `#indexOf` scans `rows()` per arrow, focus, expand/collapse and leave. Only `state-tree-tab.tsx` uses
  `TreeNavigation`; rows are state cells (hundreds, rarely thousands), so it is µs per key. Fix: check
  `rows[this.#index]` first, scan only on a miss.

## Bugs

| Id  | Title                                                                         | Where                                                                                    | Severity | Verdict   |
| --- | ----------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------- | -------- | --------- |
| B1  | `key=` on `<For>` is dropped by the JSX runtime: ~40 lists keyed by identity  | [control-flow.ts:366](../../packages/prism/ui/src/dom/control-flow.ts#L366)              | low      | confirmed |
| B2  | Past the hard limit, every agent call still snapshots and measures the scene  | [history.ts:608](../../packages/prism/ui/src/managers/history.ts#L608)                   | low      | confirmed |

- **B1**: with the automatic runtime, `<For each={x} key={f}>` compiles to `jsx(For, { each, children }, f)`
  (esbuild, `k.tsx`); `jsx()` ignores `_key`, so `For` falls back to `props.by ?? props.key ?? identity` = identity.
  `For`'s props type declares `key?`, so TypeScript accepts it (`IntrinsicAttributes` is empty). 40 `key=` sites vs 2
  `by=` in 33 files (prism ui, assistant, viewer, stories; e.g. trace.tsx:95/151, tools.tsx:44, canvas-chrome.tsx:137,
  viewer history.tsx:402). Behaviour coincides today because `For` re-creates an entry whose item changed under the same
  key; it differs only for duplicates (distinct objects with one key all render instead of one), and the test
  "renders an item again when its key stays" ([control-flow.test.tsx:108](../../packages/prism/ui/src/dom/_test/control-flow.test.tsx#L108))
  passes without keying. `hot-cases.tsx`/`runtime-bench` measure identity keying. The canvas-chrome duplicate concern
  does not hold: `sectionsOf` makes one ungrouped section (only a group `''` would collide with `undefined`). Fix:
  remove `key?` from `For`'s props (TS then rejects every site), rename to `by=`; any future in-place update by key
  would otherwise silently not apply.
- **B2**: part of P1. `#take(true)` ignores the off flag and has no memory of "too large": a scene at 600 ms or
  > 50 MB per snapshot makes each tool call of each turn block the main thread ≥ 0.5 s for a checkpoint that is then
  thrown away. Fix: P1 (4).

## Checked and dropped

- **P7 ScrollArea per scroll event**: scroll events are dispatched once per frame per target in "run the scroll steps"
  of update-the-rendering (all engines), so there is one read/write per frame already; reads happen on clean layout and
  the thumb write is a `transform`. Moving `update` to `frames()` adds a frame of lag for no gain.
- **P8 `presetIndex` JSON compare**: the always-listened path is the group header's `meta` (partial mode), which exits
  on the first differing key: 1.5 µs per change (`bench-presets.mjs`, 10 presets). The list and line-graph paths cost
  ~26 µs (40 presets × 64 colors) but run only while the presets popover is open (`derive` is lazy, `current` is read
  inside `PresetsList`'s `Show`), not during a drag.
- Corrections: P1's typical-scene cost is < 2 ms (it is medium only for large sessions and turns); P3's totals are
  bounded by mounted UI; P4 is 17–190 µs, not 0.1–1 ms; P9 is µs (state cells, not 1e5 scene rows); B1 is 40 sites,
  not 21.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/ac6cef1a-5954-43a5-ae83-c22fef35f065/scratchpad/perf-audit/prism-ui-assistant/`,
bundled with `node build.mjs X.ts X.mjs` (esbuild, `molstar-src` condition; `dom-env.ts` installs jsdom), run with
`node --expose-gc X.mjs`.

| Script                 | What                                                                       |
| ---------------------- | -------------------------------------------------------------------------- |
| `bench-ckpt.mjs`       | `JSON.stringify` with the `checkpointSize` replacer on 500/2000/8000 transforms |
| `bench-show.ts`        | `Show` vs `bind` text, make time and retained heap                         |
| `bench-show2.ts`       | `Show` vs `derive` + `insertRegion`, retained heap and toggle-all          |
| `bench-for.ts`         | `For` append / replace-last per change at 500/2000/8000 rows               |
| `bench-shallow.ts`     | `shallowEqual` vs `for…in`, `RecordState.set` with each                    |
| `bench-handlekey.mjs`  | `handleKey`'s map/filter/sort shape at 300 and 3000 commands               |
| `bench-presets.mjs`    | `presetIndex` partial group and full list modes                           |
| `k.tsx`                | esbuild transform showing `key=` as `jsx`'s third argument                 |
