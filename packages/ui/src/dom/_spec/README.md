# Mol\* DOM runtime: the specification

`@molstar/ui/dom` (`packages/ui/src/dom/`): a JSX runtime that builds real DOM once and updates it from states
(`@molstar/core/reactive`), with no virtual DOM. This spec is the rules the code relies on and why; code comments cite
its sections (`dom §2.2`). The code is the truth for how.

| File                                                             | What                                                                                                               |
| ---------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ |
| `scope.ts`                                                       | `ViewScope`, `getScope`, `requireScope`, `runInScope`, `onCleanup`, `createRoot`, `onActive`, `MaybeReactive` (§1) |
| `watch.ts`                                                       | `watch`, `bind`, `follow`, `effect`, `flush`, `frameState` (§2)                                                    |
| `frame.ts`, `work.ts`, `host.ts`                                 | the frame lane, `FrameState`, frame callbacks; work queues, `latest`; `DomHost` (§3, §8)                           |
| `jsx-runtime.ts`, `jsx-types.ts`, `render.ts`, `writes.ts`       | JSX, props, regions, event delegation, `render`; the writes (§4)                                                   |
| `control-flow.ts`, `motion.ts`                                   | `Show`, `Switch`/`Match`, `Dynamic`, `For`, `activity`, time slicing; `Motion` (§5)                                |
| `scoped.ts`                                                      | `on`, observers, `frame`, `frames`, `timeout`, `restartable`, `workQueue`, `latest`, `drag` (§6)                   |
| `element.ts`, `view.ts`, `focus.ts`, `tree-keys.ts`, `events.ts` | `ViewElement`, the view contract, focus and keyboard helpers (§7)                                                  |
| `probe.ts`                                                       | `DomProbe`, `DomCounts` (§9)                                                                                       |

## 1. Scopes

- **Every binding belongs to the current scope** (`ViewScope`), opened by `render`, a `ViewElement`, each region's
  content and each `For` entry. A component opens none: its bindings live as long as the region it was rendered in. A
  reactive binding or a watch outside any scope throws instead of leaking.
- **A `ViewScope` is core's `Scope`** (`reactive §7`: cleanups last added first, every one attempted) plus a parent, a
  depth, children and activity. Children are an intrusive list (no allocation per child, O(1) removal). Disposing
  disposes the children first (newest first), then the cleanups; a failing cleanup is reported
  (`ReactiveConfig.onUnhandledError`) and the rest run. **Walks of the children read the list again after each child**
  (disposing unlinks each one first; a sync reads its next after it, and starts over when it was disposed), never
  through a saved link: the scope's own bindings are live meanwhile, so a child's cleanup or resume can make a region
  dispose a sibling.
- **Activity**: a paused scope (`setActive(false)`, or a paused parent) holds no listener for its bindings and watches;
  resumed, each binding writes its source's current value and each watch renders once. Nothing is rebuilt. The parts
  that pause are kept apart from the cleanups, so a pause visits only them. Resuming starts the scope's own bindings,
  then its `onActive` hooks, then the children; pausing runs the hooks first. A child made while paused starts paused.
  Why: hidden UI that stays mounted does no work.
- `onActive(fn)` is for what bindings cannot pause (observers, frame loops, timers); `scope.follow(state)` makes a scope
  active while a state is true (its own subscription never pauses).

## 2. Updates

### 2.1 Bindings

- `bind(source, fn)`: a plain value is applied once; a reactive one is written now and after each change while the scope
  is active. **On a graph node** (a `State`, a derived value, a keyed value, a `FrameState`) the binding is a **sink**
  (`reactive §4.1`): a dependent slot in the node, run in the state flush once per change with a consistent value; its
  first write runs in a `batch`, so a set made by it is written after it, not nested. **Parents first**: a sink's run
  only queues the binding, and the queue is written by one sink of the runtime's own, marked behind the first binding
  queued, by depth (a stable counting sort, nothing when in order), skipping those stopped meanwhile. Why: core runs
  sinks in their sources' dependents order, which a removal reorders; a binding inside a `Show` must not run (and read a
  value its content no longer has) before the `Show` removes it. **On a foreign source with a value** (a
  `BehaviorSubject`) a subscriber in its listener list, then its current value (written again if it moved meanwhile).
  **On any other stream** a subscription, which replays its current value if it has one.
- `follow(source, fn, stopped)` tells `fn` whether a value is the first since the binding started (made or resumed) and
  calls `stopped` when it stops; `effect` is `bind` for a side effect (the probe does not count it as a write).

### 2.2 `watch`

- `watch(sources, render)`: **the first render is synchronous**; **a change marks, a microtask renders**, so many sets
  in a task (or a `batch`) are one render per watch.
- **A watch is a sink**: one dependent slot per graph source, no listener object; only a foreign source (a subject, a
  stream, a `Listenable` that is not a node) costs a subscriber. Why: a watch of N states is one object (and, for more
  than one, an array of their versions), and marking and rendering it allocate nothing (the state flush's queue, then
  the dirty array, both reused; `test:perf`). A diamond renders once with consistent values (the render reads `.value`,
  which pulls).
- **A marked watch renders only when a version moved**: it keeps its graph sources' versions as of its last render; a
  state set and set back in a batch keeps its version (`reactive §4`), so neither a watch nor a binding runs for it. A
  foreign source has no version: its change always renders.
- **Order**: dirty watches render parents first (depth, a stable counting sort only when out of order), siblings in the
  order marked; paused and disposed ones are skipped. A render that marks others has them rendered in the same flush;
  one that keeps marking itself is dropped after 100 passes, reported.
- **A render owns nothing**: it runs outside any scope (a binding made in it throws); it reads `.value` and writes
  through `writes.ts` (`setText`, `setAttr`, `setClass`, `setStyle`, `setShown`: each compares with the DOM, reads no
  layout). A render that throws is reported and the watch goes on.
- `flush()` renders now (the task lane's work first): for code that sets state and then reads the DOM. Every delegated
  event's dispatch flushes at its end. `retarget(sources)` follows other sources (a recycled row).

## 3. Scheduling

- **Microtask lane**: marked watches and the task lane's work, one microtask a burst.
- **Frame lane**: a `FrameState` keeps only the latest `offer` (nothing allocated, no listener told); the first offer in
  a frame asks for the page's one animation frame (asked again of a new `DomHost.requestAnimationFrame`: a tick pending
  on a swapped host, a test's stub dropped with its queue, would never run and the lane would stall). The tick runs the
  frame callbacks (`requestFrameCallback`, in the order asked; a cancel reaches one of the running frame too), publishes
  the offered frame states in one `batch` (each a `set`: compared, marked in the graph, so `derive` follows it as any
  state; a value over several sees none torn), runs the frame lane's work, then the microtask lane and the watches with
  a frame state among their sources. Why: input comes several times a frame; it is batched at its source, once a frame
  for the whole page. `flushFrame()` runs the tick's publish and renders now; `flush()` never publishes.
- **`DomHost`**: the microtask, frame, idle, timer and clock APIs are read through one object, whose defaults read the
  globals at each call (fake timers drive them). A test or a host without a browser replaces them.

## 4. JSX and rendering

- The automatic runtime at `@molstar/ui/dom/jsx-runtime` (`"jsxImportSource": "@molstar/ui/dom"`). A component runs once
  and returns DOM. Custom elements add their JSX props to `CustomElements` by declaration merging.
- **Namespace**: SVG tags in the SVG namespace. The tags HTML shares (`a`, `title`, `style`, `script`) are SVG when made
  under an SVG element (a function child of one, a region's or a `For` entry's render whose place is in one: not a
  `foreignObject`; a render into a fragment keeps the namespace of what renders it), when they hold SVG content (an
  `<a>` around a `<text>`), or as a static child of an SVG element: JSX makes children before their parent, so such a
  child is made again there (attributes, delegated handlers, children moved; a bound prop, listener or `ref` on it stays
  with the first element, and warns). What a component or a control flow placed directly in an SVG element renders first
  is made before its place exists: give it as a function child (`{() => <For …/>}`).
- **Props**: names with a dash, SVG attributes and a few read-only properties as attributes; other names an element has
  as properties; on a custom element, properties. `attr:`, `prop:` force either; `on:type` is a native listener. `class`
  and `classList` combine. `style` is a string or an object (numbers are pixels except unitless properties; custom
  properties by `setProperty`). `ref` gets the element after its props and children. A static value is written once, no
  binding.
- **Regions**: a reactive child writes between anchors. Text updates one node in place (anchors come only when the value
  is nodes); anything else renders anew in a child scope, the previous one disposed. A control-flow region writes the
  same function while its content stays, and the region does not render it again (a resume keeps the content).
- **Events**: bubbling events are delegated, one listener per root and type, walking the composed path to the root and
  calling handlers with `currentTarget` as their element; `stopPropagation` ends the walk; nested roots hand over to the
  next; the walk ends with `flush()`. A non-bubbling dispatch reaches its target's handler. Pointer move, over and out,
  focus, blur, enter, leave, touch and wheel stay native. Why: one listener per root instead of one per handler, and a
  click's change is on screen when its dispatch returns.

## 5. Control flow

- **`Show`** renders anew, in a fresh scope, when the value changes identity (a flag: when it flips); resumed with a
  value shallow-equal to the shown one, it keeps the content. `reveal` marks opening content with `data-reveal` and
  plays a closing back through `Motion.close` (the content is `inert` while it closes; when it plays for only some top
  elements, the content goes at once and the started animations are cancelled); `swap` runs a replacement through
  `Motion.swap`. Motion is the application's (`motion.ts`): unset, content goes at once.
- **`Switch`/`Match`** decides once every case has a value, then on each change; **`Dynamic`** renders the current
  component.
- **`For each by`**: an entry keeps its nodes and scope while its key stays **and** its item is the same value; a new
  value under a key renders anew. Entries are placed with as few moves as possible (the longest run of kept places
  stays). `index` is a state, made only when `children` takes it. A duplicate key warns once. Not `key=`: the JSX
  runtime keeps `key`.
- **Time slicing** (on in a browser, off in jsdom; `timeSlicing`): past 64 new entries in a 12 ms slice, the rest wait
  for a fresh task (a `MessageChannel` task, not a clamped timer). A slice makes only the deferred entries, never a
  reconcile of the whole list (that was quadratic). Existing entries update at once.
- **`activity(active, render)`** renders in a child scope that follows `active` and sets `inert` on its top elements
  while inactive.

## 6. Scoped attachments and gestures

- What a component attaches imperatively is released by construction: `on` (detached while paused, removed at the end),
  `observeResize`, `observeIntersection` (disconnected at the end), `frame` (the next tick; after a resume when paused),
  `frames` (a request from anywhere, once a frame, put off while paused), `timeout`, `restartable` (one wait at a time,
  held while paused; a start after the end does nothing), `workQueue` (paused with the scope), `latest` (aborted at the
  end). Each asks for its scope before it makes anything (outside one it throws, and leaks nothing); so does `held`.
- **`drag(target, { start, move, end })`**: listeners attached once with the scope; the pointer captured at the press.
  It follows the pointer that started it. It ends on a release (`up`, the only commit), and as a cancel on
  `pointercancel`, a lost capture, a move with no button held, a new press, and the scope's pause or end. A capture lost
  with no button held waits 100 ms for its release on the document (Chrome drops the capture just before it). Why: a
  release lost to a window switch leaves no listener and no drag on.

## 7. Elements, views, focus

- **`ViewElement<M>`**: renders a view of its `model` into its light DOM when connected with a model; a disconnect
  disposes the rendering on the next microtask unless it is connected again (a move keeps it). `onConnected(root)` and
  `onReleased(root)` are the application's hooks (its styles, a root class). A model it creates itself (`createModel`)
  is mounted while rendered. `defineElements` registers once.
- **The view contract** (`ViewRenderer`, `domView`, `elementView`): render into a host under the host's scope, return a
  disposer.
- **Focus and keys**: `trapTab`, `rememberFocus`, `menuKeys`, `rovingFocus` (one Tab stop, arrows on an axis; `dir=rtl`
  swaps ← and →; `rovingTabs` is it over tabs), `focusables`, `deepActiveElement`, `TreeNavigation` (an ARIA tree's
  keyboard by its rows' model, roving focus), across shadow roots; `events.ts` for composed targets. A menu's
  `aria-disabled` items take the focus (APG: their reason is read), only `disabled` or hidden ones are skipped.

## 8. Work queues

- `workQueue` wants work by key: a later want of a key replaces it (one run a burst), in the `task` (before the
  watches), `frame` or `idle` lane, ordered by `priority` (asked each run), dropped when no longer `stillWanted`, cut by
  a `budgetMs` (the rest carried; in the task lane after a macrotask). A task-lane chain of wants past 100 runs waits
  for a macrotask. A queue disposed by its own work runs nothing more; paused by it, the rest stays wanted and runs on
  the resume (`flush` runs it regardless). `latest(fn)` lets only the newest call land.

## 9. Probe

- `DomProbe` (unset: one branch a write and a render): `on(kind, scope)` is told each binding write and watch render;
  `run(kind, site, first, fn)` runs each with where it was made and whether it is a binding's first write since it
  started (made, or resumed after a pause; false for a render). `DomCounts` counts renders and writes.

## Tests

`_test/` (jsdom): the rules above in order (`render`, `watch`, `control-flow`, `scoped`, `element`). The frame tick is
run by the tests through `DomHost.requestAnimationFrame`. `pnpm -C packages/ui test:perf` (opt in, by hand): the quick
allocation tier (`*.quick.ts`, `src/_test/perf-baseline.json`): a binding's write, a watch of five states, a set back in
a batch, ten frame offers: 0 bytes a step, each with its positive control.
