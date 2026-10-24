# Slice 11: `prism/ui`, `prism/assistant`, and a light pass over `packages/ai` (v3)

Scope: `packages/prism/ui/src`, `packages/prism/assistant/src`; new this time, a light pass over the `packages/ai` hot
paths a session meets in the UI (the streamed text, transcript and tool-result retention, the Inspector's reads).
Audited at `16277d822` (= master, 2026-10-25), v2 base `0944823d4`. Node 26 + jsdom, Apple M4, shared machine at load
4–5 with ~13 other auditors: wall times ±20–40 % (given as ranges), allocations read over ≥ 200 steps (V8's
`total_allocated_bytes` moves in ~128 KB buffers, so a single step reads nothing useful). No live model runs, no network.
Read-only: no fix below is applied.

PRs since v2 that touched the slice: #2029 (slice 11 fixed: P2–P6, P9–P11, B1, B3–B5), #2023 (the Log's filter in the
manager, an error opens the Log at Problems), #2027 (extensions moved; one path), #2037 (periodic box menu in the
playback bar), #2043 (UI parts: row tooltips sized to the row, `drag()` waits for a release after a lost capture, the
band fades as a whole). `packages/ai` has no commits since v2.

## Hot paths now

| Path | Cost (fixture) | Where the time goes |
| --- | --- | --- |
| Wheel over the turn card (`shell-view.tsx:391` over `agent-turn.tsx:282`) | every wheel gesture over the card | the canvas bar's non-passive pinch guard covers the card's scrolling body (P12) |
| Turn card per streamed delta (`agent-model.ts:98` → `Card`) | 2.8–3.8 µs, 2.3 KB a delta, flat from 1k to 40k chars and from 3 to 30 runs (jsdom) | the turn record spread, ~10 derives, one text node write; the browser's relayout of the paragraph is O(text) and not measured here |
| Status task pills (`TaskPills`) | 14–46 µs, ~5.3 KB a progress event of 2 tasks (jsdom writes); 0 pills or spinners made anew | one keyed state `replace`, one pill's watch |
| Log tab per flush (`LogView`, 300 entries) | 0.55 / 0.76 / 1.75 ms per flush of 1 / 10 / 100 lines, manager alone 50–90 µs (jsdom, loaded) | `RecycledList` re-places the ~54 rows in the window (a prepend moves them all), one slot shown |
| Inspector, a settle while a turn runs | guarded: `reading.scaling.ts`, `rows.scaling.ts`; views.test (one new block and row a settle) | the trace's incremental read; the running turn's rows only |
| Harness `Store.head` per step (`state.ts:216`) | 0.09 / 0.14 / 0.31 ms a step at a window of 200 / 800 / 2000 entries | a copy and a canonical-JSON digest of the whole entry list a step (see Checked) |

## v2 → now

All v2 items were fixed in #2029; each re-checked in the code, the claims with numbers re-run (v2's `task-pills.ts`,
`log-view.ts`, `log-bound.ts` re-pointed at HEAD; the #2029 tests; `test:perf` of prism-ui; `test:scaling` of
packages/ai).

Counts: 12 v2 items (P2–P6, P9–P11, B1, B3–B5): fixed 12, partly 0, still open 0, regressed 0, dropped 0 (P7, P8
stay dropped from v1). One residual of P10 is new (P12). New bugs from fixes since v2: none found.

| v2 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| P2 | Inspector rebuilds the running turn's column every settle | fixed | #2029 | [trace.tsx:173](../../packages/prism/assistant/src/inspector/trace.tsx#L173) columns by turn id, blocks by step id (keyed states); [trace.tsx:126](../../packages/prism/assistant/src/inspector/trace.tsx#L126) rows a `For` per turn over `rowGroups`. views.test.tsx passes (a settle makes the new step's block and row and the turn's header only). Trace still not virtualized (decision: variable heights, sticky headers) |
| P3 | `Show`/`Switch`/`Dynamic` on rxjs chains | fixed | #2029 | [control-flow.ts:42](../../packages/prism/ui/src/dom/control-flow.ts#L42), `follow()` [scope.ts:688](../../packages/prism/ui/src/dom/scope.ts#L688): no rxjs on a state. Pause/resume semantics re-read (below): equivalent |
| P4 | `handleKey` maps, filters, sorts every command | fixed | #2029 | [manager.ts:773](../../packages/prism/ui/src/commands/manager.ts#L773) bucket by the key's letter and its code's, index rebuilt per list of commands and keymap ([manager.ts:440](../../packages/prism/ui/src/commands/manager.ts#L440)); bucket choice equals `chordMatches`' code fallback. keymap.test.ts passes |
| P5 | `shallowEqual` allocates four arrays | fixed | #2029 | [shallow-equal.ts:16](../../packages/prism/ui/src/reactive/shallow-equal.ts#L16) `for…in` + count; `recordSetSame` in `test:perf` (see below) |
| P6 | `composedPath()` per root per pointermove | fixed | #2029 | [render.ts:243](../../packages/prism/ui/src/dom/render.ts#L243) move/over/out native; the JSX sites taking them are few (command menu rows, placeholders, sequence strip) |
| P9 | `TreeNavigation.#indexOf` scans all rows | fixed | #2029 | [tree-keys.ts:119](../../packages/prism/ui/src/dom/tree-keys.ts#L119) last index ±1 first |
| P10 | Non-passive wheel listener on the shell root | fixed | #2029 | [shell-view.tsx:390](../../packages/prism/ui/src/views/shell/shell-view.tsx#L390) on the canvas box and the canvas bar only. The bar holds the assistant's turn card, which scrolls: P12 |
| P11 | Log list redone per flush, unvirtualized | fixed | #2029 | [log-view.tsx:94](../../packages/prism/ui/src/views/shell/log-view.tsx#L94) `RecycledList`, sources kept by the manager. Re-run: 54 rows in the DOM for 300 entries; per flush numbers in Hot paths (loaded machine; v2's after-fix 0.11–0.39 ms was on a quieter one) |
| B1 | `key=` on `<For>` dropped by the JSX runtime | fixed | #2029 | no `<For … key=` left (multi-line search over packages/ and apps/) but the `@ts-expect-error` test; 68 `by=` |
| B3 | Task pills rebuilt every progress event | fixed | #2029 | [canvas-chrome.tsx:1276](../../packages/prism/ui/src/views/shell/canvas-chrome.tsx#L1276). v2's `task-pills.ts` at HEAD: spinners made anew over 100 events 200 → **0**. Cancel reads `abort` at click; hidden by `[hidden]`, which canvas.module.css:704 makes win over `display: grid` |
| B4 | History tab relative times never refresh | fixed | #2029 | [history-tab.tsx:38](../../packages/prism/assistant/src/views/history-tab.tsx#L38) time of day |
| B5 | Plugin warnings push the messages out of the Log | fixed (by decision) | #2029 | [notifications.ts:234](../../packages/prism/ui/src/managers/notifications.ts#L234) one bound (300), problems dropped last. v2's `log-bound.ts` at HEAD: 300 entries, all 300 Plugin warnings, 0 messages: as the decision states (messages stay in `history`) |

### Fixes since v2, re-read for new bugs or costs

- **#2029 `Show`/`Switch`/`Dynamic` on `follow()`**: `Show`'s `last` persists across pauses but a resume bypasses it
  (`resumed`), as the per-subscription `distinctUntilChanged` did; a swap or a closing that ends after a pause does not
  write (`stops`); a pause during a closing cancels it and the resume writes the fallback. `Switch` waits for every case
  again after any case stops and decides once (as `combineLatest`); a plain case never stops. `Dynamic` skips a repeat
  only past the first value. `Binding.stop` calls `stopped` only for a binding that ran. No difference found.
- **#2029 `#keyIndex`**: rebuilt when the command list or the keymap layers change; lazy commands made by `#fresh()`
  change the list, so the next lookup rebuilds. Equivalent to the old filter (scope rank, then registration order;
  runnable preferred).
- **#2029 notifications**: `#flush` adds a source with an `includes` over the sources (a handful); `#bounded` is O(300)
  only past the bound. `#setTasks` keeps the tasks by id before the ids, so a new pill reads its task.
- **#2029 `shallowEqual`**: symbols now need to be enumerable on both sides (before: own on the right); a stricter
  equal, not a bug.
- **#2043 `drag()`**: two capture listeners on the document per `drag()` instance (6 call sites), active only during the
  100 ms lost-capture wait; the timer is cleared on pause and cleanup (`finish`). A release reaching both the document
  and the target ends the drag once.
- **#2043 tooltips**: a row tip's show reads `getComputedStyle` and two `getBoundingClientRect` after writing
  `paddingBlock` (two to three forced layouts per show), once per hover after the 350 ms delay: against styling.md §5's
  "no layout reads after writes" in letter, not measurable in jsdom, too rare to rank.
- **#2037 periodic box menu**: the menu (`openContextMenu`) is not closed when the playback bar unmounts (the structure
  removed while it is open); it closes on the next outside press, and an item run meanwhile acts on a stale model. Same
  pattern as the dock tab menu and the sequence menus; not recorded as a bug.
- **#2023 `logFilter`**: one filter shared by every Log view, set by the status line; intended.

## Open performance findings (ranked)

| Rank | Id | Status | Title | Where | Kind | Estimate | Effort | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | P12 | new | The canvas bar's pinch guard covers the turn card's scrolling body | [shell-view.tsx:391](../../packages/prism/ui/src/views/shell/shell-view.tsx#L391) | time | code: the card's reply and runs stop scrolling for as long as a tool call holds the main thread (reprs 0.5–20 s in the pipeline profile) | S | no |

Nothing else ranks: the per-delta, per-progress, per-flush and per-key paths measured above are flat in what they
scale with, and the two quadratic-looking costs (`Store.head`, the trace caches) are bounded by compaction (Checked).

### P12: the canvas bar's pinch guard covers the turn card's scrolling body

- **Mechanism**: #2029 moved the non-passive `wheel` guard (`keepPinch`) from the shell root to the canvas box and the
  canvas bar (`canvasBottom`, shell-view.tsx:391), whose comment says "nothing in it scrolls". The bar holds the
  `StatusLine`, and the `StatusLine` takes the assistant's turn card as `after` (shell-view.tsx:232; prism-viewer
  shell.tsx:788; assistant viewer/ui.tsx:138). The open card (`Card`, agent-turn.tsx:282) has a `ScrollArea` body
  (the streamed reply and the run list) under a `max-height: min(420px, 60cqh)` (agent-turn.module.css:15): it scrolls
  once a reply or the runs pass that. A non-passive wheel listener on an ancestor makes the compositor wait for the main
  thread on every wheel gesture over the card (v2's P10 mechanism).
- **When it bites**: exactly while the card matters: a turn's tool calls (load, representations, a surface) are the
  long main-thread tasks of the profile (cartoon 0.86–5.1 s, ball-and-stick 4.2 s, molecular surface 0.9–20 s at
  flat-192), and reading the reply or the run list by scrolling the card then does nothing until the task yields.
- **Evidence**: code (not measurable in jsdom, as P10). The guard is the only non-passive wheel listener above the card
  (`ScrollArea` has no wheel listener; the lists' are passive: recycled-list.tsx:444, virtual-list.tsx:225).
- **Fix**: keep the guard off the turn card: attach it to the status pills and the corner (`canvasClasses.corner`)
  rather than the whole bar, or wrap `after` in an element whose wheel listener is absent and move the guard to the
  bar's other children; a pinch over the card then zooms the page, as over a panel since #2029. Fix the comment either
  way. shell-watch.test.tsx can assert that no non-passive wheel listener is an ancestor of a `ScrollArea` in the bar.
- **Risk**: a pinch over the card zooms the page (accepted for panels in #2029); check Chrome and Safari on a trackpad.

## Bugs

No new bugs found. All v2 bugs fixed (table above).

## Already good, re-checked

- **Streamed text, end to end**: the harness emits one `model.delta` observation a provider frame
  (harness/loop/model.ts:122), not stored; the session appends to a rope and tests "said" on the delta only
  (session/session.ts:460, guarded by `streaming.scaling.ts`); the viewer's journal coalesces the reply
  (journal record.ts:61, slice 10); the assistant's `update` spreads one small record. Turn card probe
  (`turn-card.tsx`, `AgentTurnCard` with a fake model, 2000 deltas a point): 2.8–3.8 µs and 2.3 KB a delta at
  1k / 10k / 40k chars and 3 / 30 runs, flat in both; `RunList` re-renders no row on a delta (its `runs` derive is the
  same array). The browser's paragraph relayout per delta (O(text), the text node's data replaced) is not measured.
- **Inspector reads**: `readTrace` only while a view watches (inspector-sources.ts), reading only the new events and the
  calls that changed (memo by call and turn); `pnpm -C packages/ai test:scaling` passes at HEAD (2 files, 6 tests);
  `#rowsOf` memoized per trace read and per unchanged turn; `row` is a `find` over rows a settle (µs at 2.4k rows).
- **Task pills**: 0 spinners made anew over 100 viewer-shaped progress events (v2: 200); 14–46 µs a progress event.
- **Log**: 300 entries, 54 rows in the DOM; the manager's flush 50–90 µs; per flush cost is the window's re-placement,
  not the entries.
- **History items** (agent-model.ts `historyItemsOf`): items kept per History row in a WeakMap; read only while watched.
- `pnpm -C packages/prism/ui test:perf` at HEAD: every case's bytes a step equals its baseline (`recordSetSame` 31.98
  vs 32.00 B, the two empty symbol arrays; `For, the same rows anew` 0; `bindWrite` 0). One case, "typing in the command
  menu", failed on the 120 s test timeout (286 s wall on the loaded machine), not on bytes: 1,238,465 B a step against a
  1,239,331 B baseline. Its 1.24 MB a keystroke (737 GCs over the run) was not attributed here; jsdom's DOM costs are in it.
  It is baselined, so not a regression, but the largest number in the tier and worth a look in a browser profile.

## Checked, not worth a fix now

- **`Store.head` per step** ([state.ts:216](../../packages/ai/src/harness/loop/state.ts#L216)): each step that moves a
  pointer stores a full copy of the head's entry ids and digests their canonical JSON, so the heads hold
  O(steps × window) ids. `heads.ts` (Store and `digest` from source, two entries a step): 100 / 400 / 1000 steps at a
  window of 200 / 800 / 2000 → 10k / 160k / 1.0M ids held, 0.84 / 1.38 / 8.2 MB retained, `head()` 0.09 / 0.14 /
  0.31 ms a step. Compaction bounds the window (the model's context), so in practice it is linear in steps at a few KB a
  step, beside a network call of seconds. A head as parent + appended ids would make it O(new), at the cost of the
  content-addressed head layout (harness.md §2, §4): not proposed.
- **Tool-result retention**: the harness store keeps every entry and blob for the session (append-only by design,
  harness.md §4: in memory, content-addressed); the Inspector's trace cache holds the same objects (references, not
  copies) plus its Maps. No viewer tool returns images or other inline attachments, so results are text: bounded by
  what the model was shown.
- `session.ts:401` and `:496`: `events.some` per tool event and per planned call, over the turn's own events: O(calls²)
  a turn, at tens of calls.

## Not audited

- Browser layout and paint (the streamed paragraph, the Log's re-placement, the row tooltip's forced layouts): jsdom
  only; no browser runs in this audit.
- `packages/ai` providers, threads, narration, compaction planning; `prism/assistant` AI setup and key storage.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/31ee827f-4715-4967-bec3-085f5013d11f/scratchpad/pa3/11/`,
bundled with `node build.mjs X.ts X.mjs` (esbuild, `molstar-src`; `build2.mjs` adds the assistant's node_modules for
`@molstar/*` imports from the scratch folder; `dom-env.ts` installs jsdom).

| Script | What |
| --- | --- |
| `task-pills.ts` | v2's B3 probe at HEAD: spinners made anew over 100 progress events |
| `task-pills2.ts` | µs and bytes a progress event, 2000 events |
| `log-view.ts`, `log-bound.ts` | v2's P11 and B5 probes at HEAD |
| `log-view2.ts` | µs a flush (no macrotask) and bytes, view and manager alone (`node log-view2.mjs view|manager`) |
| `turn-card.tsx` | `AgentTurnCard` with a fake model: µs and bytes a streamed delta at 1k / 10k / 40k chars (`node turn-card.mjs <runs>`) |
| `heads.ts` | harness `Store.head`: ids and MB retained, ms a step (`node --expose-gc heads.mjs`) |
