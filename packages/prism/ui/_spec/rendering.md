# Rendering: the DOM runtime

How prism-ui turns models into DOM: its own fine-grained JSX runtime (`src/dom/`), custom elements, the view contract,
and the React/Preact adapter (`src/react/`). The code is the truth for how.

## 1. Why its own runtime

- **A component builds its DOM once and renders as a function of its state** (view = f(state); David, 2026-10-17,
  replacing per-property bindings, which were accidental). Its static structure is made once, with refs to the parts
  that change; one `watch` render reads the current values (`state.value`) and writes only what differs (text,
  `classList.toggle`, attributes, inline styles: `dom/writes.ts`). Renders are batched by one global microtask
  scheduler, so many changes in a task are one render per component. No virtual DOM, no dependency beyond rxjs (kept at
  the models' edges), no compiler beyond TypeScript/esbuild (`jsxImportSource`), real nodes any host can adopt.
- **Why**: a row bound property by property held one rxjs subscription per attribute (the scene tree's row about 16,
  with cold duplicate `pipe(map, distinctUntilChanged)` chains and a `shareReplay`), so making rows as a list scrolls
  cost memory and time per row. Measured (Internals › Performance › Trees › Rows made and the runtime bench's probes,
  Chromium, JS heap per row made after GC; script with the scroll's layout moved to the frame):

  | view                                   | heap per row (bound → watched) | script                                                         |
  | -------------------------------------- | ------------------------------ | -------------------------------------------------------------- |
  | scene tree row (1,000 components open) | 16.1 → 4.0 KB                  | 0.077 → 0.067 ms per row made                                  |
  | state tree cell (5,025 cells)          | 8.7 → 3.4 KB                   | 0.071 → 0.066 ms per row made                                  |
  | command menu row (1,000 commands)      | 2.2 → 1.2 KB                   | a keystroke 1.34 → 1.23 ms                                     |
  | params form field (500 fields)         | 3.3 → 2.4 KB                   | a slider step 0.062 → 0.043 ms                                 |
  | bench row (js-framework-benchmark, 1k) | —                              | create 5.4 → 4.5 ms; 10 unbatched relabels of all 4.8 → 2.2 ms |

  What remains per row is the DOM (making it, laying it out: about 0.1 ms a scene-tree row at the frame), `For`'s entry
  and scope, the menu's region.

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
  components, **called once to build** their DOM: what changes afterwards is written by their `watch` (§2.2) or, for
  structure, by a region (§2.4). JSX `key` is ignored (lists key through `For`).
- **Props and children may still be reactive** (`MaybeReactive<T>`: an `Observable`, a `ReadonlyState`, a writable
  `observableState`): a reactive prop is a `bind` of its own, written on each value; `undefined` binds nothing. Fine for
  a one-off (a panel's title, a story); a component's own properties and anything repeated (rows, fields) use `watch`. A
  stream child is a range between two anchors: text updates one node, anything else is rendered anew in a fresh child
  scope.
- **Properties vs attributes**: SVG, `data-*` and `aria-*` are attributes; on a custom element other names are
  properties; elsewhere a name the element has is a property (except `list`, `form`, `download`, `role`, `width`,
  `height`). `attr:` and `prop:` force either; `on:any-name` for exact event names; handlers are not reactive.
- **`class` + `classList`** combine; `style` takes a string or an object of (reactive) values, `--msp-*` included.
- **`ref`** gets the element after its props and children: how a component keeps the parts its render writes.

### 2.2 Updates: `watch`

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
  click's change is rendered when the dispatch returns, as a click by code in a test is).
- **A render owns nothing**: it runs outside any scope (a reactive `bind` in it throws), reads `.value`, writes through
  `dom/writes.ts` (`setText`, `setAttr`, `setClass`, `setStyle`, `setShown`: each compares with the DOM, reads no
  layout). Content of a few plain nodes may be replaced whole (the command menu's marked title); content with bindings
  of its own is a region.
- **Sources are explicit**: states (their `listen`: a plain callback, no rxjs Subscriber), writable states, anything
  with `listen` (`keyedSelector(...).key(k)`: told only when that key becomes or stops being the current one, so a
  highlight's move renders two rows, not all), or any stream (subscribed; what it emits while subscribing is not a
  change). A `mapState`'s `listen` is its source's ("may have changed": the render's compares keep it cheap). rxjs stays
  at the models' edges: `asReadonly` remains an `Observable`.
- **Errors**: a render that throws is logged; the others in the flush render and the watch goes on.

### 2.3 Ownership

- **Every binding belongs to the current render scope** (`Scope`, `dom/scope.ts`), opened by `render`, an `MspElement`,
  each stream region (`Show`/`Switch`/`Dynamic` branch) and each `For` entry. A component opens none: its bindings live
  as long as the region it was rendered in. Disposing runs children first, then cleanups in reverse; a failing cleanup
  is logged and the rest run.
- **API**: `onCleanup`, `watch(sources, render)`, `flush`, `bind(source, fn)`, `createRoot`, `runInScope`, `getScope`,
  `requireScope`, `toStream`, `onActive`. A reactive binding or a watch outside any scope throws instead of leaking.
- **Active and paused scopes** ([styling.md](styling.md) §5): an inactive scope (or one under an inactive parent) drops
  its `bind` subscriptions and `watch` listeners and makes them again on resume: each state replays its value (one write
  per binding) and each watch renders once; nothing is rebuilt. A source that does not replay misses what happened
  meanwhile. `Show`, `Switch` and `Dynamic` keep their content on resume; `For` reconciles and does not continue a
  time-sliced render while paused.
  - `activity(active, () => content)` renders in a child scope that follows `active` and sets `inert` on its top
    elements while inactive. Bind what does the hiding outside it, or it pauses with it.
  - `onActive(fn)` is for what bindings cannot pause (observers, rAF loops, timers, store subscriptions).
  - `render(fn, container, parent?)` renders under `parent`, so a mounted view pauses with its host's scope.

### 2.4 Control flow

`dom/control-flow.ts`; each renders into its own range and disposes what it leaves.

- **`Show`** renders again, in a fresh scope, whenever the value changes identity (a flag: when it flips). Narrow the
  stream first (`map`, `distinctUntilChanged`) when the content should stay.
- **`For each key`**: an entry keeps its nodes and scope while its key stays **and** its item is the same value; a new
  value under the same key renders it anew. `index` is a state. A duplicate key warns once.
- **`Switch`/`Match`**, **`Dynamic`**: the first truthy case; the current component.
- **Views follow records by id.** Where a list's records are replaced on every change (an agent turn's runs, History),
  run `For` over the ids (`distinctUntilChanged` on their contents) and give each row a stream of its record, so a row
  does not restart (a spinner, an open fold). Lists edited by position iterate indices likewise.

### 2.5 Rules

- **A component's own properties and repeated views are one `watch`**: build the DOM once with refs, list the states it
  shows as sources, read `.value` in the render, write through `dom/writes.ts`. What is fixed for the component's life
  (a row's node, its depth, whether it has a menu) is built in, not written.
- **Structure stays in regions** (`Show`, `For`, `Switch`, a stream child): content that appears, goes or repeats.
- **Outside a watch, read states, not values**: `{model.label}` binds, `{model.label.value}` renders once.
- **Derive in the model**; an inline `.pipe(map(…))` is fine for a one-off prop, not per row.
- **Keep identities stable where content should stay** (§2.4).
- **Bursts render once by themselves** (a microtask's worth of sets is one render); `batch` still makes a state emit
  once to its subscribers, and a short debounce helps where a burst spans tasks.
- **No component-local state that should survive a remount**: put it in a model.

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
