# Rendering: the DOM runtime

How prism-ui turns models into DOM: its own fine-grained JSX runtime (`src/dom/`), custom elements, the view contract,
and the React/Preact adapter (`src/react/`). The code is the truth for how.

## 1. Why its own runtime

- **Fine-grained bindings, no virtual DOM.** Models expose many small `ReadonlyState`s; each binds directly to the text
  node, attribute or child range that shows it, so a change touches only that DOM. No dependency beyond rxjs, no
  compiler beyond TypeScript/esbuild (`jsxImportSource`), real nodes any host can adopt.
- **Chosen over Preact-in-elements, Lit, Solid and React** on a bench: the state tree at 5,000 cells, a burst of
  updates, the PD form, the command menu. What made the big lists fast is the runtime's and the models' to keep: `For`
  slices long lists (12 ms slices on Mol\*'s `setImmediate`; `timeSlicing`, off in jsdom), off-screen rows skip layout
  (`content-visibility`; trees are virtualized now), records per node ([view-models.md](view-models.md) §3), `batch` and
  `keyedSelector`, `VirtualList`.
- **The bench stays**: Internals › Runtime performance (`src/_stories/runtime-bench.stories.tsx`) runs the
  js-framework-benchmark operations on the runtime, rows written as the app writes them; run it after a runtime change.

## 2. The runtime

### 2.1 JSX

- **Every `.tsx` names its runtime with a pragma**: `/** @jsxImportSource @molstar/prism-ui/dom */` (React files:
  `react`). The tsconfig default is not enough: bundlers and published `src/` apply their own JSX settings.
- **`jsx(type, props)` creates real nodes.** Strings make elements (SVG tags in the SVG namespace); functions are
  components, **called once**, no re-render; JSX `key` is ignored (lists key through `For`).
- **Props and children may be reactive** (`MaybeReactive<T>`: an `Observable`, a `ReadonlyState`, a writable
  `observableState`); a reactive prop writes on each value; `undefined` binds nothing. A stream child is a range between
  two anchors: text updates one node, anything else is rendered anew in a fresh child scope.
- **Properties vs attributes**: SVG, `data-*` and `aria-*` are attributes; on a custom element other names are
  properties; elsewhere a name the element has is a property (except `list`, `form`, `download`, `role`, `width`,
  `height`). `attr:` and `prop:` force either; `on:any-name` for exact event names; handlers are not reactive.
- **`class` + `classList`** combine; `style` takes a string or an object of (reactive) values, `--msp-*` included.
- **`ref`** gets the element after its props and children; with `bind` it is the escape hatch for imperative DOM.

### 2.2 Ownership

- **Every binding belongs to the current render scope** (`Scope`, `dom/scope.ts`), opened by `render`, an `MspElement`,
  each stream region (`Show`/`Switch`/`Dynamic` branch) and each `For` entry. A component opens none: its bindings live
  as long as the region it was rendered in. Disposing runs children first, then cleanups in reverse; a failing cleanup
  is logged and the rest run.
- **API**: `onCleanup`, `bind(source, fn)`, `createRoot`, `runInScope`, `getScope`, `requireScope`, `toStream`,
  `onActive`. A reactive binding outside any scope throws instead of leaking.
- **Active and paused scopes** ([styling.md](styling.md) §5): an inactive scope (or one under an inactive parent) drops
  its `bind` subscriptions and makes them again on resume, when each state replays its value: one write per binding,
  nothing rebuilt. A source that does not replay misses what happened meanwhile. `Show`, `Switch` and `Dynamic` keep
  their content on resume; `For` reconciles and does not continue a time-sliced render while paused.
  - `activity(active, () => content)` renders in a child scope that follows `active` and sets `inert` on its top
    elements while inactive. Bind what does the hiding outside it, or it pauses with it.
  - `onActive(fn)` is for what bindings cannot pause (observers, rAF loops, timers, store subscriptions).
  - `render(fn, container, parent?)` renders under `parent`, so a mounted view pauses with its host's scope.

### 2.3 Control flow

`dom/control-flow.ts`; each renders into its own range and disposes what it leaves.

- **`Show`** renders again, in a fresh scope, whenever the value changes identity (a flag: when it flips). Narrow the
  stream first (`map`, `distinctUntilChanged`) when the content should stay.
- **`For each key`**: an entry keeps its nodes and scope while its key stays **and** its item is the same value; a new
  value under the same key renders it anew. `index` is a state. A duplicate key warns once.
- **`Switch`/`Match`**, **`Dynamic`**: the first truthy case; the current component.
- **Views follow records by id.** Where a list's records are replaced on every change (an agent turn's runs, History),
  run `For` over the ids (`distinctUntilChanged` on their contents) and give each row a stream of its record, so a row
  does not restart (a spinner, an open fold). Lists edited by position iterate indices likewise.

### 2.4 Rules

- **Read states, not values, in render**: `{model.label}` binds, `{model.label.value}` renders once.
- **Derive in the model**; an inline `.pipe(map(…))` is fine for trivial formatting (components run once).
- **Keep identities stable where content should stay** (§2.3).
- **Combine, then debounce**: a row that depends on several streams combines them; where bursts come, a short debounce
  makes one change render once.
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
