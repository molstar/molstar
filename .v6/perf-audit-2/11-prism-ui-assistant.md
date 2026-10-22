# Slice 11: `prism/ui` and `prism/assistant` (v2)

The Prism DOM runtime and the assistant UI: the reactive core (`reactive/state.ts`, keyed state and selectors, the
frame lane), the DOM runtime (`dom/scope.ts`, `dom/render.ts`, `dom/control-flow.ts`), virtual lists, the undo
`HistoryManager`, `CommandManager`, the shell (status line, Log), and in the assistant the turn card, the History tab
and the Inspector. Outside the profiled data pipeline: costs scale with UI items (rows, commands, trace steps, tasks,
log entries); the undo checkpoint scales with the scene. Audited at `0944823d4` (v1 base `60fe4c2f1`), Node 25 +
jsdom, Apple M4. Each finding re-read with its callers; disputed impacts re-measured.

Since v1: #2008 (incremental checkpoints, 10-step cap, 30 s pause over the hard limit, size counted in ms), #2009 and
#2010 (Log: log-only sources with a bound per source, filter icons), #2021 (command menu own forms,
`invocation.setParam`, a pinch guard on the shell root), #2002 (an icon), #2015 (test lanes). `prism/assistant/src`
has no commits since v1. The reactive core, `For`/`Show`, `RecycledList`/`VirtualList` and the frame lane are
unchanged and still well optimized (O(1) listener removal, allocation-free counting-sort flush, `For`'s unchanged-list
fast path and LIS placement, slot pooling).

## Hot paths now

| Path                                                    | Runs                                                   | Cost driver                                                                 |
| ------------------------------------------------------- | ------------------------------------------------------ | --------------------------------------------------------------------------- |
| Wheel over any shell panel (`shell-view.tsx:385`)       | every wheel gesture                                    | a non-passive listener on the shell root: the gesture waits for main thread |
| Status task pills (`StatusLine` → `TaskPill`)           | every Mol\* task progress event (≤ 10/s per task)      | every pill disposed and rebuilt                                             |
| Log tab (`LogView`)                                     | every log flush (a microtask burst) while the Log shows | O(rows ≤ 300) reconcile + one row build per new entry, unvirtualized        |
| `HistoryManager.#take` → checkpoint port                | every undoable run, every agent tool call              | O(cells) cached walk: 0.5 ms at 1401 cells, 1.7 ms at 4201 (was ~4 ms/MB)   |
| `CommandManager.handleKey`                              | every keydown                                          | O(commands): 27 µs at 437                                                   |
| `Show`/`Switch`/`Dynamic`                               | per mounted instance                                   | rxjs chain memory                                                           |
| Inspector overview and trace list                       | every 80 ms settle while a turn runs (diagnostic page) | running turn's steps, session rows                                          |

The command menu's per-keystroke ranking is ~0.14 ms at 437 commands (`bench-rank.mjs`): nothing recorded. Outside
this slice, for slice 10: `viewer/src/sources/target-format.ts:47` `complete()` walks every unit of every structure per
keystroke in a target field (11k units at flat-192) and keeps 40.

## v1 → now

Fixed 2, partly fixed 0, still open 7, regressed 0, dropped 2 (as audited; all 7 fixed since in #PRNUM, see [Fixed in #PRNUM](#fixed-in-prnum)).

| v1 id | Title                                                                            | Status     | PR    | Note                                                                                                                                                                                                                                                                   |
| ----- | -------------------------------------------------------------------------------- | ---------- | ----- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| P1    | History: full snapshot + JSON size on every run, also in a turn and past the hard limit | fixed      | #2008 | `SerializedTrees` (checkpoint.ts:99) serializes and measures each transform object once and shares it: size 0.04 ms, checkpoint 1.7 ms at 4201 cells / 3.2 MB (full stringify 14 ms). Size counted in ms (history.ts:654-657), 30 s pause (:650), cap of 10 (`#setEntries` :627). |
| P2    | Inspector rebuilds the running turn's column every settle; trace not virtualized | fixed      | #PRNUM | trace.tsx:151 identity-keyed `For`, `turn.steps.map` (:193); trace.tsx:95 mounts every row.                                                                                                                                                                            |
| P3    | `Show`/`Switch`/`Dynamic` built on rxjs chains                                   | fixed      | #PRNUM | control-flow.ts unchanged since v1.                                                                                                                                                                                                                                    |
| P4    | `handleKey` maps every command to an object per keydown                          | fixed      | #PRNUM | manager.ts:739-749; 27 µs at 437 commands (re-measured shape).                                                                                                                                                                                                         |
| P5    | `shallowEqual` allocates four arrays per compare                                 | fixed      | #PRNUM | shallow-equal.ts:13-14 unchanged.                                                                                                                                                                                                                                      |
| P6    | Delegated dispatch calls `composedPath()` per root for every pointermove         | fixed      | #PRNUM | render.ts:243-251, :269, :293 unchanged.                                                                                                                                                                                                                               |
| P7    | ScrollArea per scroll event                                                      | dropped    | –     | Dropped in v1: scroll events already come once per frame per target.                                                                                                                                                                                                  |
| P8    | `presetIndex` JSON compare                                                       | dropped    | –     | Dropped in v1: 1.5 µs on the listened path; the rest runs only while the popover is open.                                                                                                                                                                              |
| P9    | `TreeNavigation.#indexOf` scans all rows per key                                 | fixed      | #PRNUM | tree-keys.ts:118-123 unchanged.                                                                                                                                                                                                                                        |
| B1    | `key=` on `<For>` is dropped by the JSX runtime                                  | fixed      | #PRNUM | control-flow.ts:369 still declares `key?`; 49 `key=` sites vs 2 `by=` (multi-line count, incl. stories and tests).                                                                                                                                                    |
| B2    | Past the hard limit, every agent tool call still snapshots and measures          | fixed      | #2008 | history.ts:650: within 30 s of a snapshot over the hard limit, no snapshot; reset by `clear()`, `setCheckpointing(true)`, any checkpoint under the limit. Tests in history.test.ts.                                                                                   |

## Performance findings (all fixed in #PRNUM)

| Rank | Id  | vs v1      | Title                                                                    | Where                                                                                      | Kind        | Impact | Estimate                                                                         | Effort | Verdict   | Frozen |
| ---- | --- | ---------- | ------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------ | ----------- | ------ | -------------------------------------------------------------------------------- | ------ | --------- | ------ |
| 1    | P10 | new        | Non-passive wheel listener on the shell root blocks every panel's scroll | [shell-view.tsx:385](../../packages/prism/ui/src/views/shell/shell-view.tsx#L385)          | time        | medium | panels scroll on the compositor during long main-thread tasks (0.5-20 s)         | S      | fixed #PRNUM | no     |
| 2    | P11 | new        | Log tab redoes its list on every log flush, unvirtualized                | [log-view.tsx:103](../../packages/prism/ui/src/views/shell/log-view.tsx#L103)              | time        | low    | 0.03 / 0.5 / 4.3 ms per flush of 1 / 10 / 100 lines (jsdom) → O(new rows)       | M      | fixed #PRNUM | no     |
| 3    | P3  | still-open | `Show`/`Switch`/`Dynamic` built on rxjs chains                           | [control-flow.ts:75](../../packages/prism/ui/src/dom/control-flow.ts#L75)                  | memory      | low    | −1.3 KB per mounted `Show` (0.3-3 MB total)                                      | M      | fixed #PRNUM | no     |
| 4    | P2  | still-open | Inspector rebuilds the running turn's column every settle                | [trace.tsx:151](../../packages/prism/assistant/src/inspector/trace.tsx#L151)               | time+memory | low    | O(S) block builds per settle → O(new steps); trace DOM bounded to a window       | M      | fixed #PRNUM | no     |
| 5    | P4  | still-open | `handleKey` maps, filters and sorts every command per keydown            | [manager.ts:739](../../packages/prism/ui/src/commands/manager.ts#L739)                     | time+memory | low    | 27 µs → O(commands on the chord) per keydown, no allocation                      | S      | fixed #PRNUM | no     |
| 6    | P5  | still-open | `shallowEqual` allocates four arrays per compare                         | [shallow-equal.ts:13](../../packages/prism/ui/src/reactive/shallow-equal.ts#L13)           | time+memory | low    | ~20x faster compares, zero allocation                                            | S      | fixed #PRNUM | no     |
| 7    | P6  | still-open | Delegated dispatch calls `composedPath()` per root per pointermove       | [render.ts:293](../../packages/prism/ui/src/dom/render.ts#L293)                            | time+memory | low    | no allocation per pointer move while no handler of that type is mounted          | S      | fixed #PRNUM | no     |
| 8    | P9  | still-open | `TreeNavigation.#indexOf` scans all rows per key                         | [tree-keys.ts:118](../../packages/prism/ui/src/dom/tree-keys.ts#L118)                      | time        | low    | O(rows) → O(1) per arrow key (µs today)                                          | S      | fixed #PRNUM | no     |

### P10: pinch guard makes every panel's wheel scroll wait for the main thread

- **Mechanism**: #2021 added `on(root, 'wheel', keepPinch, { passive: false })` on the shell root (shell-view.tsx:382-385)
  so a trackpad pinch (ctrl+wheel) does not zoom the page. A non-passive wheel listener on an ancestor makes the
  compositor send each wheel gesture to the main thread and wait for it (Chrome: the first event of a gesture is
  blocking, the rest queue behind it). That covers every ScrollArea in the shell: the scene tree's `RecycledList`, the
  sequence panel, the Log, forms. The lists' own listeners are passive on purpose (recycled-list.tsx:444,
  virtual-list.tsx:225) so the compositor can scroll the already-made rows ahead of the main thread, which
  `RecycledList`'s margin logic counts on. Not a body/document listener, so Chrome's passive-by-default intervention does
  not apply.
- **Evidence**: code. The main thread is busy 0.5-20 s at a time during loads and rebuilds in the profile
  (molecular-surface 0.9-20 s, ball-and-stick 4.2 s, cartoon 0.86-5.1 s at flat-192): a wheel over any panel then does
  nothing until the task yields. Not measurable in jsdom.
- **Fix**: attach the non-passive guard only where nothing scrolls: the canvas area (which already blocks, Canvas3D's
  input) and the splash/empty canvas, optionally the rail and canvas chrome; not the dock panels. Over a panel either let
  the browser have the pinch (page zoom) or accept the trade-off explicitly. Keep `gesturestart`/`gesturechange` on root
  (not scroll-blocking). Adding the listener lazily on the first ctrl+wheel does not work: that first event already went
  to a passive path. A maintainer decision: #2021 fixed page zoom by a pinch over a side panel.
- **Risk**: a pinch over a panel may zoom the page again. Check Chrome and Safari on a trackpad.

### Low

- **P11**: each `#flush` (notifications.ts:168) sets a new `logged`; with the Log shown, `sources` (Set + sort over all
  entries, log-view.tsx:48), `shown`, and a `For` over up to 300 rows (200 messages/Plugin + 100 Transforms) run, plus
  one row build (Icon, 4 spans, `time.format`, `toISOString`) per new entry and one dispose per evicted one. Rows are
  kept (the list holds the same objects), so it is O(new) builds + O(rows) reconcile, not a full rebuild.
  `log-view.ts`: 0.03 ms per 1-line flush, 0.5 ms per 10, 4.3 ms per 100 (jsdom, net of a 1.3 ms macrotask); the
  manager's flush alone is ~10 µs. Only while the Log is shown. Fix: `sources` kept incrementally in the manager;
  `RecycledList` for the entries (fixed row height); throttle the Log's view of `logged` to the frame lane. Risk: the
  `role="log"` semantics and the sticky filter line; test log-view.test.tsx.
- **P3**: each `Show` wraps `toObservable(when).pipe(distinctUntilChanged())` in a new Observable (control-flow.ts:75);
  `Switch` `combineLatest` (:221), `Dynamic` `map` (:188). v1 `bench-show2`: 4.1 vs 2.8 KB retained per instance. Fix: a
  `derive`-based region when `when` is a state; rxjs only for foreign Observables. Risk: pause/resume and the closing
  animation; tests control-flow, hidden-quiet, leaks.heap, waste budget.
- **P2**: `model.turns` changes every 80 ms settle while a turn runs; the identity-keyed `For` (B1) rebuilds the running
  turn's `Column` and all its `StepBlock`s (trace.tsx:193), O(S²) per turn; the trace list (:95) mounts every row.
  v1: < 0.3 ms per settle at S 5-30, 0.2 ms reconcile at 2000 rows. Diagnostic page. Fix: inner `For` over step ids with
  `keyedState`, the outer one by turn id (note: `by=` alone does not help: `For` re-renders an item whose key stays but
  value changes); `RecycledList` for the rows. Risk: `reveal`/`scrollToIndex`, j/k; views.test.tsx, rows.scaling.ts.
- **P4**: `map` (one object per command) + `filter` + `sort` before the `typing` early-out (manager.ts:739); `#boundKey`
  (:336) scans every binding per key while a lazy contribution is stale. `bench-handlekey-437.mjs`: 18.5 / 26.7 /
  194 µs at 300 / 437 / 3000 commands. Fix: `Map<chordId, CommandDefinition[]>` invalidated on `#commands` or
  `keymap.layers`; chord normalization must equal `chordMatches` (mac variants, `code` fallback).
- **P5**: `Reflect.ownKeys(...).filter(...)` on both sides (shallow-equal.ts:13-14), the default compare of every record
  state. v1: 920 → 46 ns (12 keys). Fix: `for…in` + `Object.hasOwn` count; symbols only if `getOwnPropertySymbols` is
  non-empty.
- **P6**: pointer/mouse move/over/out delegated (render.ts:243-251); each root calls `composedPath()` and two `indexOf`
  per bubbling event (:293), plus a capture listener (:269), with no handler mounted outside the open command menu. Fix:
  a live per-type handler count, return at 0; or attach those types natively.
- **P9**: linear scan of `rows()` per arrow/focus/expand (tree-keys.ts:118); state-tree tab only, µs. Fix: check
  `rows[#index]` first.

## Bugs (all fixed in #PRNUM)

| Id  | vs v1      | Severity | Title                                                                         | Where                                                                                      | Verdict            |
| --- | ---------- | -------- | ----------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------ | ------------------ |
| B3  | new        | medium   | Task pills rebuilt on every progress event: spinner restarts, Cancel misses   | [canvas-chrome.tsx:1130](../../packages/prism/ui/src/views/shell/canvas-chrome.tsx#L1130)  | fixed #PRNUM |
| B5  | new        | low      | Log: Plugin warnings share the messages' bound and push the messages out      | [notifications.ts:188](../../packages/prism/ui/src/managers/notifications.ts#L188)         | fixed #PRNUM |
| B1  | still-open | low      | `key=` on `<For>` is dropped by the JSX runtime: 49 lists keyed by identity   | [control-flow.ts:369](../../packages/prism/ui/src/dom/control-flow.ts#L369)                | fixed #PRNUM |
| B4  | new        | low      | History tab 'just now' / 'N min ago' labels never refresh                     | [history-tab.tsx:97](../../packages/prism/assistant/src/views/history-tab.tsx#L97)         | fixed #PRNUM |

- **B3**: `StatusLine` renders `<For each={derive(status.tasks, filter)}>` keyed by identity; `TaskPill` (:1263) builds
  its DOM once from a plain object. Each Mol\* progress event (task-manager.ts:59, 100 ms rate) makes
  viewer-sources.ts:134 a new `TaskInfo`, and viewer.ts:509 maps every task to a new object, so every pill (also the
  unchanged ones) is disposed and rebuilt. The new spinner Icon restarts its 1 s CSS animation (canvas.module.css:581)
  up to 10 times a second; the Cancel button is replaced too, so a press that spans an update (a ~100 ms click against a
  100 ms update) has its `mousedown` on a removed button and no `click` on the new one: Cancel can take several tries
  during a busy load (inferred from the event model, not run in a browser). `task-pills.ts`: 2 tasks, 100 viewer-shaped
  `setTasks`, 200/200 spinners recreated, 1.5 ms per event incl. a macrotask. Fix: follow tasks by id as `RunList` does
  (agent-turn.tsx:110-125): ids derived with `sameIds`, a `keyedState` by id, `TaskPill` a `watch` writing label, bar
  width and the abort button's visibility, reading `abort` at click time. Optionally keep the previous `TaskInfo` in
  `scanList` when nothing changed.
- **B5**: `#bounded` counts a source in its own bucket only if it has a `keep` (notifications.ts:183-193); any other
  log-only source shares `historySize` (200) with the status messages. viewer.ts:505 logs plugin warnings with no
  `keep`. `LogOptions.keep` says "default: the history's size" per source, and `#bounded`'s comment says a burst "does
  not push the messages out". `log-bound.ts`: 5 messages, 300 Plugin warnings → Log holds 0 messages, 200 Plugin, 1
  Transforms (the messages stay in `history`, but the Log's Messages filter shows nothing). Before any source with a
  `keep`, the fast path caps all entries at 200 together. Fix: key = `m.source`, limit = `#keep.get(m.source) ?? size`,
  also in the fast path; a test with 300 keep-less log() calls.
- **B1**: `<For key={f}>` compiles to `jsx(For, props, f)`; `jsx()` (jsx-runtime.ts:16) ignores `_key`, so `keyOf` is
  identity. TypeScript accepts it (control-flow.ts:369 declares `key?`). Today it differs only for duplicate keys, and
  any future in-place update by key silently does not apply. Fix: drop `key?` from the props (TS then flags every
  site), rename to `by=`.
- **B4**: `ago(time)` is evaluated once when a row is made (history-tab.tsx:97 runs, :135 turns); rows are kept by id
  and nothing re-reads the time, so 'just now' stays 'just now' for an hour while newer rows say '12 min ago'. Latent:
  `HistoryTab` is exported by `@molstar/prism-assistant` but no app in the repo mounts it. Fix: absolute times (as the
  Log), or one scoped minute tick while the tab is shown.

## Fixed in #PRNUM

Each finding fixed with a test in its lane; numbers before → after on the same machine (Node 25, jsdom, Apple M4), the
audit's bench scripts re-pointed at the branch.

| Id  | Fix                                                                                                                       | Before → after                                                                  | Test                                                     |
| --- | ------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- | -------------------------------------------------------- |
| P10 | The non-passive wheel guard on the canvas box (the splash lets the pointer through) and its bottom bar; `gesture*` stay on the root     | code (not measurable in jsdom): no wheel listener that is not passive over a panel | shell-watch.test.tsx                                   |
| P11 | The Log a `RecycledList` of one-line recycled rows; its sources kept by the manager (`logSources`)                       | per flush of 1 / 10 / 100 lines, net of a macrotask: 0.16 / 0.58 / 4.5 → 0.11 / 0.25 / 0.39 ms; 300 → 54 rows in the DOM | log-view.test.tsx                     |
| P3  | `Show`, `Switch`, `Dynamic` follow their sources with `follow()` and write through `regionWriter()`: no rxjs on a state | 4.1 → 2.8 KB retained a `Show` (4000 over states), made in 55 → 39-49 ms        | activity.test.tsx, control-flow.test.tsx, leaks.heap     |
| P2  | Overview: a column per turn id, a block per step id (keyed states); list: a `For` per turn over the model's `rowGroups`   | blocks made over 30 steps of the running turn: 495 → 30; a settle at 300 turns (2432 rows): 1.26 → 0.64 ms | views.test.tsx               |
| P4  | Commands indexed by their key's letter once per list of commands and keymap; a keydown reads its key's bucket            | 17 → 0.06 µs a keydown at 437 commands (125 → 0.05 µs at 3000)                  | keymap.test.ts                                           |
| P5  | `for…in` with an own-key count; symbols by `getOwnPropertySymbols`                                                         | 925 → 120 ns, 656 → 32 B a compare (12 keys; the 32 B are the two empty symbol arrays) | lib.test.ts, primitives.quick (`recordSetSame`)    |
| P6  | Pointer and mouse move, over and out handlers native, not delegated                                                      | `composedPath()` per move: 4 → 0 (two roots) while no handler is under the pointer | render.test.tsx                                       |
| P9  | `#indexOf` checks the last index and its neighbours before a scan                                                         | 9,991 → 5 key lookups an arrow key at the end of a 5,000-row tree               | tree-keys.test.tsx                                       |
| B1  | `For` has no `key` prop (TypeScript rejects `key=`); the 49 sites pass `by=` (ids or indices, unique in their lists)      | –                                                                               | control-flow.test.tsx (`@ts-expect-error`)               |
| B3  | `NotificationsManager.lineTaskIds` / `lineTask(id)`; a `TaskPill` made once, a watch writes it, Cancel reads `abort` at the click | spinners made anew over 100 progress events of 2 tasks: 200 → 0            | status-pill.test.tsx                                     |
| B4  | The time of day (as the Log), not "N min ago"                                                                             | –                                                                               | history-tab.test.tsx                                     |
| B5  | One bound for the whole Log (`logSize`, 300; the per-source `keep` and the 100 transforms go); warnings and errors dropped last | 5 messages + 300 Plugin warnings: the bound holds 300, all problems             | notifications.test.ts                                    |

Decisions and departures from the suggested fixes:

- **B5** (the user's direction, 2026-10-22): the Log keeps the last 300 entries of every source together, as one bound,
  in a virtual list. Warnings and errors stay visible inside it by being dropped last: past the bound the oldest info and
  success entries go first, and problems only when they alone fill it. So a load's burst of transforms never pushes a
  problem out; a burst of 300 plugin warnings does fill the Log with warnings (the messages stay in `history`). The
  suggested per-source key (`#keep.get(source) ?? size`) is not done: there is one bound now.
- **P2**: the trace list is not virtualized. Its rows differ in height and a turn's header is sticky over its rows; a
  `RecycledList` needs one row height and positions rows absolutely (no sticky). What cost per settle was the `For`
  reconciling every row (1 ms at 2,400 rows, every 80 ms): grouping by turn makes it the running turn's rows only. The
  first mount of a long trace stays O(rows), time-sliced by `For`.
- **P3**: a region source is followed by `follow()` (a `bind` told whether a value is the first since it started, and
  when it stops), not by a `derive`-based region as suggested: `Show` keeps per-subscription state (a closing animation,
  the value it resumed with) that a `derive` does not have.
- **P10**: over a panel, a card, the sequence band or the rail a pinch zooms the page again (before #2021); ⌘0 or a pinch back over a panel
  undoes it. Not checked on a trackpad in Chrome and Safari.
- **B4**: absolute times, not a minute tick: a clock re-rendering every row would keep an idle app busy (Prism's "idle
  is idle").

## Checked and dropped

- **P1 (residual) history checkpoint per agent call**: the auditor's "a second full JSON.stringify, 20-200 ms per call
  at 5-50 MB" no longer holds: #2008's `SerializedTrees` shares unchanged transforms and memoizes their sizes.
  `bench-ckpt2.ts` (headless viewer, N × 1cbs, tree unchanged): checkpoint 0.09 / 0.49 / 1.65 ms and size 0.01 /
  0.02 / 0.04 ms at 281 / 1401 / 4201 cells (0.2 / 1.1 / 3.2 MB), vs a full stringify of 1.1 / 6.0 / 14 ms. ~0.4 µs per
  cell per call; at the 50 MB hard limit (~65k cells) about 25 ms, on an agent call that takes seconds. Not worth a
  version cache.
- The command menu ranking per keystroke: 0.14 ms at 437 commands (`bench-rank.mjs`).
- #2021's own-form code (`setParam`, `intoForm`, `formControls`) and #2008's markdown entity decoding: per key or per
  render, no repeated work found.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/1f4e2fbb-da6c-4948-aacb-385a4dab465c/scratchpad/perf-audit-2/prism-ui-assistant/`,
bundled with `node build.mjs X.ts X.mjs` (esbuild, `molstar-src` condition; `dom-env.ts` installs jsdom).

| Script                    | What                                                                                    |
| ------------------------- | --------------------------------------------------------------------------------------- |
| `bench-ckpt2.ts`          | checkpoint port per call (`checkpoint`, `size`, full stringify) on N × 1cbs, headless   |
| `task-pills.ts`           | `StatusLine`, 2 tasks, 100 viewer-shaped `setTasks`: spinners recreated, ms per event   |
| `log-view.ts`             | `LogView` with 300 rows, flushes of 1 / 10 / 100 Transforms lines                       |
| `log-bound.ts`            | B5: 5 messages + 300 keep-less Plugin warnings, what the Log keeps                      |
| `bench-handlekey-437.mjs` | v1's `handleKey` shape at 300 / 437 / 3000 commands                                     |
| `bench-rank.mjs`          | the command menu's ranking shape at 437 commands                                        |

v1's scripts (`bench-show2`, `bench-shallow`, `bench-for`, `bench-ckpt`) are in the v1 scratchpad,
`/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/ac6cef1a-5954-43a5-ae83-c22fef35f065/scratchpad/perf-audit/prism-ui-assistant/`.
