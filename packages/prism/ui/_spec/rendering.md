# Rendering: the DOM runtime

How prism-ui turns models into DOM: its own fine-grained JSX runtime (`src/dom/`), custom elements, the view contract,
and the React/Preact adapter (`src/react/`). The code is the truth for how.

## 1. Why its own runtime

- **A component builds its DOM once and renders as a function of its state** (view = f(state); David, 2026-10-17,
  replacing per-property bindings, which were accidental). Its static structure is made once, with refs to the parts
  that change; one `watch` render reads the current values (`state.value`) and writes only what differs (text,
  `classList.toggle`, attributes, inline styles: `dom/writes.ts`). Renders are batched by one global microtask
  scheduler, so many changes in a task are one render per component. No virtual DOM, no compiler beyond
  TypeScript/esbuild (`jsxImportSource`), real nodes any host can adopt.
- **The reactive base is plain JS; rxjs is the boundary** (David, 2026-10-17). A state is its latest value plus a plain
  array of listeners; `derive` follows states by `listen`; `watch` and `bind` listen. rxjs stays where Mol\*'s and the
  plugin's events come in (§2.2).
- **Why**: a row bound property by property held one rxjs subscription per attribute (the scene tree's row about 16,
  with cold duplicate `pipe(map, distinctUntilChanged)` chains and a `shareReplay`), so making rows as a list scrolls
  cost memory and time per row. The same scene-tree row in three shapes, on the plain base (Internals › Performance ›
  Trees › Rows made, Chromium, 1,000 components open, far jumps of 27 rows; JS heap per row made after GC; script with
  the scroll's layout at the frame):

  | scene-tree row | heap per row | script per row made | a row step | | -------------------------------------------- |
  ------------ | ------------------- | ---------- | --- | | rxjs pipes per attribute (the JSX before) | 14.7 KB | 0.076
  ms | 0.32 ms | | the same JSX, `derive` in place of the pipes | 8.8 KB | 0.071 ms | 0.29 ms | | one `watch` render |
  3.4 KB | 0.058 ms | 0.30 ms |

  (On the rxjs base before it, the pipe row was 16.1 KB.) Converted to `watch`, other views gained alike: a state tree
  cell 8.7 → 3.4 KB, a command menu row 2.2 → 1.2 KB (a keystroke 1.34 → 1.23 ms), a params form field 3.3 → 2.4 KB (a
  slider step 0.062 → 0.043 ms); the bench's rows create 1,000 in 5.4 → 4.5 ms and relabel all ten times, unbatched, in
  4.8 → 2.2 ms. What remains per row is the DOM (making it; laying it out, about 0.1 ms a scene-tree row at the frame),
  `For`'s entry and scope, the menu's region.

- **The update path** of a mounted row, from a state's change to the DOM showing it (the same rows, 48 in the window;
  Chromium cross-origin isolated, a 5 µs clock; `__mspRows.updates`/`updateCase`, `watch-bench/update-path.mjs`). Time
  is script until the DOM is right, for `watch` after its microtask flush; heap is JS allocated per update (no GC in
  between); calls are listener calls, derived recomputes and `watch` renders per update (a burst: per burst), in
  brackets those on rows the update does not change (from 48 rows against 12):

  | case (per update)                     | rxjs pipes             | `derive` JSX            | `watch`, DOM compares    | `watch` (JS compares)   | `watch`, keyed (now)  |
  | ------------------------------------- | ---------------------- | ----------------------- | ------------------------ | ----------------------- | --------------------- |
  | a row's label                         | 0.6 µs, 85 B, 1        | 0.6 µs, 85 B, 1         | 3.0 µs, 229 B, 2         | 1.0 µs, 189 B, 2        | 1.1 µs, 189 B, 2      |
  | a row's hidden (class, title, icon)   | 3.6 µs, 466 B, 4       | 3.3 µs, 466 B, 6        | 1.9 µs, 183 B, 2         | 1.1 µs, 142 B, 2        | 1.2 µs, 142 B, 2      |
  | a row's swatch (a style)              | 0.4 µs, 71 B, 1        | 0.4 µs, 71 B, 3         | 1.2 µs, 199 B, 2         | 0.8 µs, 199 B, 2        | 0.9 µs, 199 B, 2      |
  | the pick moves between two rows       | 4.9 µs, 937 B, 97 (82) | 2.1 µs, 933 B, 101 (92) | 59 µs, 5.2 KB, 97 (82)   | 4.9 µs, 3.3 KB, 97 (82) | 2.1 µs, 1.0 KB, 6 (0) |
  | a reflow far below (places map anew)  | 86 µs, 275 B, 53 (45)  | 82 µs, 322 B, 101 (91)  | 130 µs, 4.3 KB, 101 (91) | 87 µs, 2.9 KB, 101 (91) | 82 µs, 335 B, 5 (0)   |
  | 100 label sets in one task (10 rows)  | 25 µs, 6.3 KB, 100     | 26 µs, 6.3 KB, 100      | 18.5 µs, 3.0 KB, 110     | 7.2 µs, 2.8 KB, 110     | 8.1 µs, 2.8 KB, 110   |
  | 100 sets, each in its own task, a set | 14 µs, 1               | 16 µs, 1                | 38 µs, 2                 | 19 µs, 2                | 20 µs, 2              |

  With bindings the DOM is right when `set` returns; with `watch` a read right after `set` sees the old DOM until the
  microtask (or `flush()`), and `set` itself costs 0.1–1.3 µs (marking). A row's own change costs about the same either
  way (1 µs). A state every row reads (the pick) runs every row: two cold rxjs chains a row, a `derive` per row, or a
  render per row; a `watch` render that compares its sources' values with what it last wrote, in JS, costs what the
  chains do (4.9 µs); one that compares with the DOM (`getAttribute`, `classList`) costs twelve times that. Keyed (the
  last column: a row follows `keyedSelector(picked).key(node)` and `keyedState` places, view-models.md §3), the rows an
  update does not change do nothing: a pick is 6 calls (2 renders) and 2.1 µs, a reflow below the window 5 calls and no
  render. A burst in one task is one render a row (7 µs against 25). A set alone in its task pays the microtask and the
  flush (19 µs against 14–16). The reflow's 80 µs is the tree's flattening, not its rows.

- **No frame budget for rows, measured** (a scheduler that makes rows within a few ms a frame, cancels what left the
  window and shows placeholders was asked for, to build only if a budget would bite): a thumb drag over the scene tree
  (4,000 rows, top to bottom and back in 1.5 s, the mouse on the scrollbar's thumb; `watch-bench/thumbdrag.mjs`) costs
  0.17–0.18 ms of UI script and 0.19–0.21 ms of native DOM work a frame with recycled slots (0.30–0.52 and 0.44–1.07
  with keyed rows); no frame over 20 ms after the first run's warm-up, no frame showing a blank row. A far jump is
  0.6–0.9 ms of script (p90 3.1–3.4 ms with its layout): under any budget worth having. The window moves in the frame
  lane (the scroll position is a frame state), so a drag's events make one window a frame; nothing more is scheduled.
- **While a structure animates, the UI is 2–6% of the main thread's work** (0.01–0.02 of 0.4–0.5 ms a frame; the table
  and the hidden-UI checks are in [styling.md](styling.md) §5).
- **Recycled rows for the scene tree** (`RecycledList`, the default there; `VirtualList` keeps keyed rows elsewhere): a
  pool of slots, each with its DOM and scope for good, each placed by its own transform; a window that moves points
  slots at the items that came in (`VirtualSlot.show`, a `watch` re-pointed by `retarget`) and makes, removes and
  reorders no node. Measured against keyed rows, both `watch` rows (Rows made, 1,000 components open, Chromium on an M4,
  two runs; `watch-bench/recycle.mjs`):

  | scene tree                               | keyed rows (`VirtualList`)                                  | recycled slots (`RecycledList`) |
  | ---------------------------------------- | ----------------------------------------------------------- | ------------------------------- |
  | a thumb-drag jump, median / p90          | 4.7–4.9 / 5.0–5.4 ms                                        | 1.7–2.4 / 3.1–3.4 ms            |
  | of it script / layout                    | 1.8–1.9 / 2.8–2.9 ms                                        | 0.6–0.9 / 1.2–1.7 ms            |
  | a small step (a row), mean / p95         | 0.31–0.39 / 1.8–2.5 ms                                      | 0.12–0.15 / 0.7–0.8 ms          |
  | 120 jumps: allocated / kept after GC     | 2.26 MB / 0.52 MB                                           | 1.10–1.17 MB / 0.43 MB          |
  | 120 jumps: GC / native (DOM) / UI script | 22.8 / 69–70 / 63–67 ms                                     | 12.7–13.0 / 18–23 / 45–47 ms    |
  | a 2 s wheel scroll: UI script + native   | 21–24 + 6–7 ms                                              | 18–20 + 1 ms                    |
  | a real wheel (60 + 30 steps)             | no long animation frame; one 33 ms frame in one of two runs | the same                        |

  The costs, and how they are met: while scrolling, the DOM order is not the screen order; **at rest** (the margin's
  rest signal, 120 ms after the last scroll) one pass moves the slots' nodes into screen order, so a screen reader
  browsing a still list reads its rows in order (a focused row gets its focus back; measured after far jumps: 0.73–0.95
  ms median, 1.7 ms at most, about 69 nodes moved, once a rest; `__mspRows.settleCost`). Otherwise (a row's place is its
  `aria-posinset`, `aria-setsize` and `aria-level`; the keyboard moves by the rows' model, `TreeNavigation`); a
  re-pointed slot closes its menu, moves the hover highlight to the node it shows now and its keyboard registration
  (`TreeNavigation.slot`), and a focus in it goes to the tree as when a row leaves the DOM; a put-away slot is hidden
  and follows nothing. Tests: `components/virtual-list/_test/recycled-list.test.tsx`, the scene tree's.

- **Hand-rolled, not `@preact/signals-core`**: backing `State` and `derive` with signals measured the same time and more
  memory (the `watch` row 3.7 vs 3.4 KB, the `derive` row 9.5 vs 8.8 KB: an effect per listener), needed this batch's
  and queue's semantics given up, and its auto-tracking misses sources that are not signals.
- **Binding writes are immediate** (`bind`, inside `batch`), only `watch` renders are queued: queued through the
  scheduler, a burst of 100 label sets on 48 shown rows took 0.017 instead of 0.039 ms and a pick the same 0.01 ms, not
  worth a field and a closure per binding and a DOM behind its state.
- **Chosen over Preact-in-elements, Lit, Solid and React** on a bench (the state tree at 5,000 cells, a burst of
  updates, the PD form, the command menu); a Preact spike for rows was dropped (2026-10-17). What made the big lists
  fast is the runtime's and the models' to keep: `For` slices long lists (12 ms slices on Mol\*'s `setImmediate`;
  `timeSlicing`, off in jsdom), lists are virtualized (`VirtualList`), records per node
  ([view-models.md](view-models.md) §3), `batch` and `keyedSelector`.
- **The benches stay**: Internals › Runtime performance (`src/_stories/runtime-bench.stories.tsx`: the
  js-framework-benchmark operations, bound rows and `__mspWatchRows`; `__mspMenu`, `__mspForm`) and Internals ›
  Performance › Trees › Rows made (`__mspRows`: far jumps, small steps, the heap, a wheel), driven from a browser
  (Playwright with CDP for GC and the heap). They mount outside Storybook's root: its highlight addon re-scans the
  page's styles on every node added there (~13 ms). Run them after a runtime change.

## 2. The runtime

### 2.1 JSX

- **Every `.tsx` names its runtime with a pragma**: `/** @jsxImportSource @molstar/prism-ui/dom */` (React files:
  `react`). The tsconfig default is not enough: bundlers and published `src/` apply their own JSX settings.
- **`jsx(type, props)` creates real nodes.** Strings make elements (SVG tags in the SVG namespace); functions are
  components, **called once to build** their DOM: what changes afterwards is written by their `watch` (§2.3) or, for
  structure, by a region (§2.5). JSX `key` is ignored (lists key through `For`).
- **Props and children may still be reactive** (`MaybeReactive<T>`: a state, a derived value, a writable
  `observableState`; an `Observable` as the exception, §2.2): a reactive prop is a `bind` of its own, written on each
  value; `undefined` binds nothing. Fine for a one-off (a panel's title, a story); a component's own properties and
  anything repeated (rows, fields) use `watch`. A stream child is a range between two anchors: text updates one node,
  anything else is rendered anew in a fresh child scope.
- **Properties vs attributes**: SVG, `data-*` and `aria-*` are attributes; on a custom element other names are
  properties; elsewhere a name the element has is a property (except `list`, `form`, `download`, `role`, `width`,
  `height`). `attr:` and `prop:` force either; `on:any-name` for exact event names; handlers are not reactive.
- **`class` + `classList`** combine; `style` takes a string or an object of (reactive) values, `--msp-*` included.
- **`ref`** gets the element after its props and children: how a component keeps the parts its render writes.

### 2.2 States and the boundary

`reactive/state.ts`; a state is "the latest value plus listeners", and views do no rxjs on it.

- **`State<T>`** (what `observableState` returns; its call sites keep their API): `value` (the latest set, also while an
  earlier set is still being told: a listener's `set(value + 1)` or a record's `update` builds on it), `set` (compared
  with the latest value; in a `batch`, held and told once at its end), `listen(fn) → unlisten` (a plain array), a
  record's `update`, `asReadonly` (a view without `set`, made on first use).
- **Re-entrancy**: a value told while a state, a derived value or a key (`keyedState`, `keyedSelector`) is telling its
  listeners (a listener sets it again) is queued and told after (`Listened.tell`), so every listener gets the values in
  order and ends on the latest; told at once, nested, a listener after the setter would end on the older value. A
  subscriber (`asObservable`) set to the state while it gets the current value gets the new one after.
- **`derive(source | sources, fn, same = Object.is)`**: lazy and hot. Unlistened, `value` computes from the sources'
  current values and nothing is held; listened, it follows its sources by `listen`, keeps the value and tells only when
  it differs; shared by every binding that uses it. It replaces `x.pipe(map(f), distinctUntilChanged())` and
  `combineLatest([a, b]).pipe(map(...))`; `mapState` is `derive` with a shallow comparison.
- **rxjs becomes a state at the model or source boundary**, with `toState(stream, initial)`, once per value, never per
  JSX binding (a binding that wraps its own stream owns an rxjs subscription again). Only code that consumes Mol\*'s and
  the plugin's events (`plugin.events`, cell updates, the viewer's `sources/`) keeps rxjs inside, and exposes states.
  The other way, a state is an `Observable` when subscribed (`asObservable`: one listener per subscriber), for the code
  that still pipes it.
- **Views receive states and derived values.** `bind` listens, then writes `value` (a set made by that write is written
  after it, not missed); pausing a scope unlistens, resuming does the same again. An `Observable` given to JSX still
  works through the old path (a subscription that replays): the exception, for a stream with no state behind it.

### 2.3 Updates: `watch`

`watch(sources, render)` (`dom/scope.ts`) in a component's scope:

- **The first render is synchronous**, when `watch` is called (after the DOM it writes is built).
- **A change marks, a microtask renders**: a source's change puts the watch in one global dirty set; a microtask flushes
  it, so many sets in a task (or a `batch`) are one render per watch. Dirty watches render in scope depth order (parents
  first; siblings in the order they changed); disposed and paused ones are skipped. A render that marks others (it sets
  state) has them rendered in the same flush; one that keeps marking itself is dropped after 100 passes, with an error.
- **Paused, it listens to nothing** (`Scope.setActive`, `activity`); resumed, it renders once, at once. Made in a paused
  scope, its first render waits for the resume.
- **`flush()`** renders now what is waiting: for code that sets state and then reads the DOM in the same task.
  `VirtualList.scrollToIndex` and `TreeNavigation`'s keys flush; so does every delegated event's dispatch at its end (a
  click's change is rendered when the dispatch returns, as a click by code in a test is). The pointer's move, over and
  out are not delegated (`dom/render.ts`): they come at pointer rate over the whole page and few elements take them, so
  their handlers are native listeners and render in the microtask.
- **A render owns nothing**: it runs outside any scope (a reactive `bind` in it throws), reads `.value`, writes through
  `dom/writes.ts` (`setText`, `setAttr`, `setClass`, `setStyle`, `setShown`: each compares with the DOM, reads no
  layout).
- **A render compares in JS first**: it keeps what it last wrote per source (or per derived value) and touches the DOM
  only for what changed. A shared source (the picked node, a tree's places) renders every row, and most write nothing:
  compared with the DOM, a pick over 48 rows took 59 µs; compared in JS, 4.9 µs (§1, the update path). Content of a few
  plain nodes may be replaced whole (the command menu's marked title); content with bindings of its own is a region.
- **Sources are explicit**: states and derived values (their `listen`), writable states, anything with `listen`
  (`keyedSelector(...).key(k)`: told only when that key becomes or stops being the current one, so a highlight's move
  renders two rows, not all), or, as the exception, a stream (subscribed; what it emits while subscribing is not a
  change, nor is a key's first answer over a stream that replays).
- **Errors**: a render that throws is logged; the others in the flush render and the watch goes on.
- **Two lanes.** The microtask lane is for logic and model state: what a task changes renders once, at its end. The
  frame lane is for input: a trackpad sends wheel and pointer events at 120 Hz or more, several a frame, so input-fed
  state is batched at its source. A **frame state** (`frameState(initial)` in a component, `FrameState` in a source,
  `reactive/frame-state.ts`) keeps only the latest `offer` (no listener told, nothing allocated); the first offer in a
  frame asks for one `requestAnimationFrame`, shared by every frame state, which publishes each once, then runs the
  microtask lane, then the frame lane: the watches with a frame state among their sources (a watch of a frame state and
  a model state renders once in the frame, after both). Paused or disposed watches are skipped in both lanes; a
  component's frame state drops its pending offer with its scope. `flush()` runs the microtask lane only (it runs at
  every delegated event, and would publish per event); `flushFrame()` runs the frame's tick now (tests, and
  `scrollToIndex`, whose row the keyboard focuses next). Ten offers in a frame cost one publish and one render, and no
  allocation (`primitives.quick.ts`).
- **One frame request a page** (`reactive/frame-state.ts`): what a component would ask `requestAnimationFrame` for (a
  label's placement, a canvas's draw, a drag's step, the splash's simulation, Seq*'s frames) is a frame callback
  (`requestFrameCallback`, `cancelFrameCallback`; `frame(fn)` in a component, released with its scope, run after its
  resume when it is paused), run in the frame's tick before the frame states publish. A canvas resized from a resize
  observer is the exception that draws at once: its size wiped it, and the observer runs past the frame's tick, so a
  frame callback would paint a blank frame first (the splash blinked through a window's resize; `PrismField.redraw`). A
  cancel reaches a callback of the frame running now too (one callback cancelling another, a `frames().cancel()`), so a
  cancelled one never runs. A hover sweep in the viewer asked for up to 3 frames a browser frame (the hover label's,
  Seq*'s scheduler's and the frame lane's), now 1 (the in-app browser, a predicted model: 20 requests in 12 frames
  before, 13 in 14 after). A loaded viewer at rest asks for none and queues no microtask
  (`viewer/src/ui/_test/frames.quick.ts`, with controls: a loop left running, a frame per event). Mol*'s own canvas loop
  is not Prism's, but at rest it asks for none either: the viewer turns on `PluginConfig.General.PauseLoopWhenIdle`
  (`options.ui.pauseWhenIdle`), so the loop sleeps once Canvas3D has been idle for 250 ms and wakes on its `wake`
  (viewer README decision 29).
- **What is a frame state**: a value fed by input events: the virtual lists' scroll position, the canvas's hover label
  (`sources/pointer.ts`, from Mol*'s label at every pick), a drag (a dock's splitter, a slider). Not one: what an event
  decides once (a click, a key), what the model says; a resize observer's callback already comes once a frame, after
  layout. Prism's views handle no wheel events of their own (the browser scrolls).

### 2.4 Ownership

- **Every binding belongs to the current render scope** (`Scope`, `dom/scope.ts`), opened by `render`, an `MspElement`,
  each stream region (`Show`/`Switch`/`Dynamic` branch) and each `For` entry. A component opens none: its bindings live
  as long as the region it was rendered in. Disposing runs children first, then cleanups in reverse; a failing cleanup
  is logged and the rest run.
- **API**: `onCleanup`, `watch(sources, render)`, `flush`, `bind(source, fn)`, `createRoot`, `runInScope`, `getScope`,
  `requireScope`, `toStream`, `onActive`. A reactive binding or a watch outside any scope throws instead of leaking.
- **Active and paused scopes** ([styling.md](styling.md) §5): an inactive scope (or one under an inactive parent) drops
  its bindings' listeners and `watch`es and makes them again on resume: each binding writes its state's current value
  and each watch renders once; nothing is rebuilt. A source that does not replay misses what happened meanwhile. `Show`,
  `Switch` and `Dynamic` keep their content on resume; `For` reconciles and does not continue a time-sliced render while
  paused.
  - `activity(active, () => content)` renders in a child scope that follows `active` and sets `inert` on its top
    elements while inactive. Bind what does the hiding outside it, or it pauses with it.
  - `onActive(fn)` is for what bindings cannot pause (observers, rAF loops, timers, store subscriptions).
  - `render(fn, container, parent?)` renders under `parent`, so a mounted view pauses with its host's scope.
- **Imperative attachments are scoped** (David, 2026-10-17; `dom/scoped.ts`): a component attaches a listener with
  `on(target, type, fn, options)` (detached while its scope is paused, removed when it ends), an observer with
  `observeResize(el, fn)` or `observeIntersection(el, fn)` (disconnected at the end), a frame with `frame(fn)` and a
  timer with `timeout(fn, ms)` (cancelled at the end), or `restartable(fn)` for one a callback in an effect starts again
  and again (one wait at a time, held while the scope is paused and run in full after it resumes, dropped at the end); a
  component's `frameState` drops its pending offer with its scope, and the frame lane skips disposed and paused watches.
  Released by construction, not by memory: a raw `addEventListener`, `new ResizeObserver`/`IntersectionObserver`,
  `setInterval`, `requestAnimationFrame` or a frame callback by hand in a component file fails
  (`ui/src/_test/scoped-attach.test.ts`, strict since every view moved to the helpers: 93 in 30 files at the start):
  every `.tsx`, and the `.ts` view code (`components/`, `views/`, and helpers a view calls in its scope: the tooltips,
  `TreeNavigation`, `rovingTabs`), a `.ts` that attaches for a reason listed in the test with it (a press-attached
  `dragToReorder` for code with no scope, a splash simulation's loop). A listener added in an event handler (no scope
  there) is attached once instead, scoped, and acts while a drag is on; a frame asked for from a handler or a binding is
  `frames(fn)`'s `request()` (made with the scope: once a frame, put off while paused, cancelled at the end).
- **Gestures** (`drag(target, { start, move, end })`, `dom/scoped.ts`): the press, the moves and every end are listeners
  attached once with the scope, acting while a drag is on; the pointer is captured at the press. Chrome sometimes drops
  the capture just before the release and sends the `pointerup` a moment later to the element under the pointer
  (measured, a PAE box drag): a capture lost with no button held waits a moment (100 ms) for that release, read on the
  document (capture phase, attached once), and it is the release (`up`); a press anywhere meanwhile, or no release, is a
  cancel. Outside that wait the document's releases are not the drag's (another gesture's never commits it). A drag ends
  on `pointerup`, `pointercancel`, `lostpointercapture`, a move with no button held (a release that never came), a new
  press, and the scope's pause and end, so a release lost to a window switch or a native dialog leaves no listener and
  no drag on. **Only a release commits** (David, 2026-10-17): `end`'s `how` is `up` for the drag's own pointer let go
  and `cancel` for every other end, and what a drag decides on its release (a row's move, a seek) happens on `up` only;
  a cancel leaves things as they were (what a drag applied live as it moved, a crop or a level, stays). **It follows the
  pointer that started it**: another pointer's moves, release and cancel are not its own, and its press starts nothing
  while the first still holds the capture (a second finger); a press after the first pointer ended unseen (its capture
  gone) starts anew. The scroll area's thumb is the first; each gesture's test presses, moves and loses its release
  twenty times and counts the listeners on every target (`testing/gestures.ts`). **A move reads no layout**: what a drag
  hit-tests against is read at the press, before its first write, and again only after a scroll or a resize (a row
  reorder's midpoints, then a binary search a move; `reorderRows`, guarded in `test:scaling`).
- **Views leave nothing behind** (`ui/src/dom/_test/leaks.heap.tsx`, `test:heap`): a keyed list with rows, a recycled
  list, the command menu, the mock shell with a tab shown and a context menu, mounted and disposed again and again, keep
  no heap outside code (n and 4n cycles, a frame between them) and leave no listener on the window, the document or its
  body (a patched `addEventListener` counts them); a view that leaves a document listener is caught. It found the mock
  shell's `FullscreenModel` left undisposed (a `fullscreenchange` listener a mount).

### 2.5 Control flow

`dom/control-flow.ts`; each renders into its own range and disposes what it leaves.

- **`Show`** renders again, in a fresh scope, whenever the value changes identity (a flag: when it flips). Narrow the
  value first (`derive`) when the content should stay.
- **No rxjs inside the regions**: `Show`, `Switch` and `Dynamic` follow their sources with `follow` (`dom/scope.ts`: a
  binding told whether a value is the first since it started, and when it stops) and write through `regionWriter`
  (`dom/render.ts`); a state is followed by its `listen`, only a foreign stream by a subscription. Why: a `Show` is in
  most views, and its `pipe(distinctUntilChanged())` chain kept 4.1 KB a `Show` against 2.8 KB without (perf audit 2,
  11/P3). `Switch` decides as `combineLatest` did: once every case has a value, and once on a resume.
- **`For each by`**: an entry keeps its nodes and scope while its key (`by`, else the item) stays **and** its item is
  the same value; a new value under the same key renders it anew. Not `key=`: the JSX runtime keeps `key` for itself, so
  `For` has no such prop and TypeScript rejects it (it was silently dropped). `index` is a state. A duplicate key warns
  once. A long list renders over time slices (`timeSlicing`): a slice makes only the entries the last one deferred, each
  before the next entry placed, never a reconcile of the whole list again (that was quadratic over the slices;
  `dom/_test/for-slices.scaling.tsx`).
- **`Switch`/`Match`**, **`Dynamic`**: the first truthy case; the current component.
- **Views follow records by id.** Where a list's records are replaced on every change (an agent turn's runs, History),
  run `For` over the ids (`distinctUntilChanged` on their contents) and give each row a stream of its record, so a row
  does not restart (a spinner, an open fold). Lists edited by position iterate indices likewise.

### 2.6 Rules

- **A component's own properties and repeated views are one `watch`**: build the DOM once with refs, list the states it
  shows as sources, read `.value` in the render, write through `dom/writes.ts`. What is fixed for the component's life
  (a row's node, its depth, whether it has a menu) is built in, not written.
- **Structure stays in regions** (`Show`, `For`, `Switch`, a stream child): content that appears, goes or repeats.
- **Outside a watch, read states, not values**: `{model.label}` binds, `{model.label.value}` renders once.
- **Derive in the model**, or with `derive` in the view; never a pipe on a state.
- **No rxjs in views** (`ui/src/_test/views-without-rxjs.test.ts`): `views/` in ui, viewer and assistant, ui's
  `components/` and the extensions' Prism UI import no `rxjs` and call no `.pipe(`; any does fails (strict since every
  view moved: 357 uses in 42 files at the start). Models, sources and the boundary helpers (`toState`, `asObservable`)
  are outside it.
- **A derive computes when read unlistened**: with no listener it holds nothing, so each read computes again (and a
  subscription reads, then listens: twice). A stateful function in one (`keepRows`) must be keyed by its inputs, or the
  derive listened to before it is read.
- **Keep identities stable where content should stay** (§2.5).
- **Bursts render once by themselves** in a `watch` (a microtask's worth of sets is one render); `batch` makes a state
  tell its listeners once (bindings write at once), and a short debounce helps where a burst spans tasks.
- **No component-local state that should survive a remount**: put it in a model.
- **Work a change asks for is wanted, not done** (§2.7): a recomputation that many changes in a task ask for (a tree
  flattened again) is a work queue's key; an answer that arrives late (a search, a count, a fetch) is `latest`.

### 2.7 Work: `workQueue` and `latest`

`reactive/work-queue.ts`; `dom/scoped.ts` ties each to the current scope.

- **`workQueue({ lane, budgetMs })`**: `want(key, run, { priority, stillWanted })` asks for `run` in the lane's next
  run, and a later want of the same key replaces it: the last write wins, so a burst of changes costs one run. Each run
  asks `stillWanted` before running (false drops it), orders by `priority` (higher first; a function is asked again each
  run), and with a `budgetMs` stops when the budget is spent, at least one done: the rest carries to the next run.
  Paused with its scope (what waited runs on resume), dropped when it ends; `flush()` on the queue runs all of it now.
  Nothing is allocated a want of a known key or a run (`primitives.quick.ts`, `primitives.heap.ts`: `workWant`,
  `workFrameRun`, 0 B a step): each key keeps its entry, as a map's deletes and sets allocate when it rehashes.
- **Errors and loops**: a piece of work that throws (or its `priority`, `stillWanted`) is logged, and the rest of the
  run and the lane's other queues go on, as a failing render does: one failure must not leave the others enlisted and
  never run. The task lane runs in microtasks, which let no timer, event or frame in, so it yields: a queue's run wanted
  by a task-lane run 100 times in a row (work wanting itself, two queues wanting each other), a task-lane budget's
  carry-over (at once: a budget is there to yield), and the work left when `flush()` drops watches that kept marking,
  wait for a macrotask; the lane's other queues go on.
- **Lanes**: `task` runs in the runtime's microtask flush before the watches it feeds (`flush()` runs it: code that
  reads the DOM after a change sees the work done); `frame` runs in the frame's tick, after the frame states publish and
  before the renders; `idle` when the page is idle (`requestIdleCallback`, its deadline the budget).
- **Which tool**: a value fed by input is a frame state (§2.3: the latest value, once a frame); work that changes ask
  for is a queue's key; a question whose answer may come late is `latest`. The virtual lists need no budget and no queue
  of rows (§1: a jump's rows cost 1.7 to 2.4 ms in Chromium, a thumb drag 0.14 to 0.17 ms a frame): they place their
  window once a frame and make what it needs then.
- **The users**: the scene tree's flattening (a `task` queue, below); the command menu's filter (its query a frame
  state: keys typed faster than frames filter once, for the latest text; a key acting on the rows, or a pick, runs the
  frame first; a row's highlight renders at once, its title's marks in the frame); the canvas's hover label
  (`sources/pointer.ts`, a frame state over Mol*'s labels). Already latest-wins, and left as they are: the command
  menu's preview and count of the highlighted row (`CommandManager.browse`: an `AbortController` a highlight, after 80
  ms), the Select group's counts among them, each command's targets resolved once. To move to `latest` as they are
  touched: fetches (the loaders') and exports.
- **The scene tree flattens once a burst** (`task` lane): each shown node's `isExpanded` and `children` asked for a
  flattening of every row, so expanding 1,000 components at once flattened 1,000 times (401 ms in jsdom; collapsing
  them, 393 ms); wanted once a task, 19.5 and 24 ms. Its first flattening is at the mount (the rows are there when it
  returns). Tests that change the scene and read the rows flush first, as for any watch.
- **The margin ahead of the scroll grows with its speed** (`components/virtual-list/window.ts`): the compositor scrolls
  a wheel or a flick before the list hears, so what it uncovers that was not made shows the skeleton for a frame. Ahead,
  the rows where the scroll will be two frames on at its speed now, at least what it was (half a screen), at most three
  screens; only a wheel's or a touch's speed counts (a thumb's drag and `scrollToIndex` set the position from script, in
  the same frame as the rows). Real wheel steps in Chromium on the scene tree (1,000 components, a screen of 23 rows, 40
  steps down and up), at the scroll event before the list moves (what the compositor scrolled onto):

  | Wheel step | Scrolls onto a skeleton, before | after    | Rows a skeleton (mean), after | UI ms a frame | Slots at most |
  | ---------- | ------------------------------- | -------- | ----------------------------- | ------------- | ------------- |
  | 600 px     | 1 of 40                         | 1 of 40  | 0.4                           | 0.08          | 62 → 79       |
  | 1,200 px   | 40 of 40 (22 rows each)         | 2 of 40  | 1.1                           | 0.15          | 62 → 104      |
  | 2,400 px   | 40 of 40                        | 40 of 40 | 20 (past three screens)       | 0.21          | 62 → 104      |

  After the frame's paint, no row place was empty in any run. A thumb's drag stays at 0.14 ms of UI a frame, no blank
  frame. The first step of a fast wheel still lands past the window (the speed is not known yet).

- **`latest(fn)`**: `fn(signal, ...args)`; each call aborts the one before it, and a call a newer one (or `abort()`)
  superseded resolves to undefined, so an old answer never lands after a newer one. The scoped one aborts its call in
  flight when its scope ends.

## 3. Custom elements

- **`MspElement<M>`** (`dom/element.ts`) renders a view of its `model` into its own light DOM once, when connected with
  a model; another model renders again; a model set before upgrade is taken over.
- **Moves keep state**: disconnecting disposes on a microtask, cancelled by reconnecting.
- **Styles**: on connect it adopts the registered styles into its document or shadow root and releases them when gone; a
  top-level element becomes `.msp-root` ([styling.md](styling.md) §3).
- **An element does not own its model** unless it creates it (`createModel()`), then mounts it while rendered.
- **Importable outside a browser** (extends `HTMLElement` only where it exists). Tags are `msp-`, registered by an
  explicit `defineElements([...])` (no import side effects; a registered tag is left alone). JSX typing by augmenting
  `MspElements` in `dom/jsx-types.ts`.

## 4. The view contract

Views a shell hosts (tabs, sections, cards, pages) implement `ViewRenderer<Ctx>` (`dom/view.ts`):
`mount(host, ctx, scope) => dispose`, where `scope` is the hosting scope the shell pauses while the view is hidden.

- `domView(Component)`: the built-in runtime; follows the scope by itself.
- `reactView(Component, { createRoot })` (`react/react-view.ts`): a `display: contents` container with a React or Preact
  root; passes `active`. The caller injects `createRoot`, so prism-ui has no React dependency.
- `elementView(tag, property = 'ctx')`: a custom element given the context.
- The shell mounts renderers with `ViewHost` (`views/shell/view-host.tsx`), disposed with the region that shows it.

## 5. React and Preact

- **`react/`** (one module each, no index): `useObservable`, `useObservableProp`, `useObservableSelection`,
  `useReactiveModel`, `reactView`, `MspElementHost` (sets properties on an element for React 18 and `preact/compat`).
- **One build serves both**: `useSyncExternalStore`; a Preact app aliases `react` to `preact/compat`. React is an
  optional peer; nothing outside `react/` imports it. Memoize sources: a new observable per render starts over (the
  hooks warn).
- **Tested on both**: `react/_test/` on React 18, `react/_test/preact.test.tsx` with `react` mocked to `preact/compat`;
  Storybook's Internals › React and Preact.
