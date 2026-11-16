# React bindings

`@molstar/ui/react`: core's reactive sources in React components, React views where the UI mounts views, and the UI's
elements in a React tree. React (from 18) is an optional peer dependency; `preact/compat` works as written (`react`
aliased to it). Nothing else in the package imports React.

## 1. Reading a source: `useReactive`

- A source is a state (`State`, `derive`, a `BehaviorSubject`, ours or rxjs's: it has a current `value`), a stream (an
  `Observable`), or nothing. The hook is `useSyncExternalStore` over a store per source that holds only the raw value;
  the selection is made by the render's `getSnapshot` (§2).
- **A source with a current value**: the store reads `value`; it is followed through core's `listen`, which tells
  changes only. Why: subscribing to a state replays its value; delivered to the host's store, that replay makes it check
  (and, under a comparer that says "changed", render) a second time on every mount. Reading `value` gives the first
  render its value, and a change between the render and the subscription is caught by the host's check.
- **A stream without one**: subscribed; every new value tells the host, including one emitted while subscribing (a cold
  `of(x)`). Why: React reads the snapshot again after it subscribes, but Preact's `useSyncExternalStore` (and React's
  shim) checks before subscribing and never after, so a value emitted then would be lost. `defaultValue` seeds a
  source's lifetime until its first value; it never replaces a value, `undefined` included.
- **Errors**: a stream that ended in an error (`hasError`) keeps its last value; any other throw from `value` (a derived
  value with a bug) reaches the host (an error boundary). A throw while a set is told is core's to report.
- **The source's identity is its lifetime**: the store is kept per source in refs (the committed render's, and one a
  render not committed yet made), not in `useMemo`, which React may drop; a render React discards never replaces the
  committed source's store. A new source (an inline `derive(...)` or `.pipe(...)`) starts over and shows only its
  default; one new on committed renders in a row warns once. Memoize it.

## 2. Selecting a part: `useReactiveProp`

- `selector` picks the part; `compare` (shallow by default, core's `sameShallow`) decides whether a new selection is the
  last one: then the last one's reference is kept and nothing renders. Why: a component that shows one field of a record
  state is not rendered for the others.
- The selector and the comparer are the render's own (no writes during render): the last selection is kept in a cell
  keyed by what made it (the store, the selector, the raw value), so a read from a discarded render only misses. A new
  selector (inline: new each render) selects again on its render; memoize an expensive one (`useCallback`). A comparer
  that changed alone selects nothing again.

## 3. Models: `useReactiveModel`

`mount` while the component is mounted, `unmount` when it is not (Strict Mode's remount mounts again); never `retire`:
retiring is the creator's.

## 4. Views and elements

- `reactView(component, { createRoot })`: a `ViewRenderer` (dom `view.ts`) that renders a React component into the host,
  with `active` (false while the view is hidden but kept mounted, after its scope) so it can stop following sources
  (`useReactive(active ? x : undefined)`). The caller passes `createRoot` (`react-dom/client` or
  `preact/compat/client`). Disposed, its container leaves at once and the root unmounts in a microtask: never
  synchronously inside a React render or commit (a React app hosting a UI element whose view is a React view).
- `UiElementHost`: one of the UI's elements in a React tree, its properties (models, callbacks) set through a ref, since
  React 18 would set them as attributes. Only values that changed are written, and a property no longer given is set
  back to `undefined` (a stale callback is not kept).
