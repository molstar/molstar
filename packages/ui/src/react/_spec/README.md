# React bindings

`@molstar/ui/react`: core's reactive sources in React components, React views where the UI mounts views, and the UI's
elements in a React tree. React (from 18) is an optional peer dependency; `preact/compat` works as written (`react`
aliased to it). Nothing else in the package imports React.

## 1. Reading a source: `useReactive`

- A source is a state (`State`, `derive`, a `BehaviorSubject`, ours or rxjs's: it has a current `value`), a stream (an
  `Observable`), or nothing. The hook is `useSyncExternalStore` over a small store per source.
- **A source with a current value**: `getSnapshot` reads `value`; `subscribe` is core's `listen`, which tells changes
  only. Why: subscribing to a state replays its value; delivered to React's store, that replay would make React check
  (and, under a comparer that says "changed", render) a second time on every mount. Reading `value` gives the first
  render its value, and React reads the snapshot again after it subscribes, so a change between the render and the
  subscription is not lost.
- **A stream without one**: subscribed; what it emits while subscribing (a cold `of(x)`) is kept without telling React
  (React reads the snapshot after subscribing), later values tell it. `defaultValue` seeds a source's lifetime until its
  first value; it never replaces a value, `undefined` included.
- A current-value stream that ended in an error throws on `value`: the last snapshot stays.
- **The source's identity is its lifetime**: a new source (an inline `derive(...)` or `.pipe(...)`) starts over and
  shows only its default. A source new on renders in a row warns once. Memoize it.

## 2. Selecting a part: `useReactiveProp`

- `selector` picks the part; `compare` (shallow by default, core's `sameShallow`) decides whether it changed. An equal
  selection keeps its previous reference and renders nothing. Why: a component that shows one field of a record state is
  not rendered for the others.
- The selector and the comparer may be new each render: the latest is used, on the next value and on the render itself
  (a selector that depends on a prop selects again), and the subscription stays.

## 3. Models: `useReactiveModel`

`mount` while the component is mounted, `unmount` when it is not (Strict Mode's remount mounts again); never `retire`:
retiring is the creator's.

## 4. Views and elements

- `reactView(component, { createRoot })`: a `ViewRenderer` (dom `view.ts`) that renders a React component into the host,
  with `active` (false while the view is hidden but kept mounted, after its scope) so it can stop following sources
  (`useReactive(active ? x : undefined)`). The caller passes `createRoot` (`react-dom/client` or
  `preact/compat/client`).
- `UiElementHost`: one of the UI's elements in a React tree, its properties (models, callbacks) set through a ref, since
  React 18 would set them as attributes.
