# Rendering: web components, JSX, React and Preact

How prism-ui turns view models into DOM. The model layer ([view-models.md](view-models.md)) is renderer-free; this
document proposes the renderer for the built-in UI and defines how other renderers plug in.

State (mocks): every mock is built on option A (the `dom/` runtime, §4), the one decided on (§3,
[../spike.md](../spike.md)).

## 1. Requirements

1. **Granular updates.** View models expose many small `ReadonlyState`s. A change should touch the DOM it affects, not
   re-render a panel.
2. **TSX authoring** with type-checked props.
3. **Embeddable pieces.** Hosts (React apps, plain pages, MVS stories, PDBe/RCSB wrappers) place parts of the UI
   (viewport, sequence, PAE plot, snapshot stepper, a PD form) without adopting the whole shell.
4. **Framework choice for contributors.** Extension and app authors can write views in React or Preact on top of the
   same view models.
5. **No compiler beyond TypeScript/esbuild.** The library is built with `tsc -b` (ESM `lib/`) and apps with esbuild; no
   Babel transforms.
6. **Survive moves.** Docking moves panels between regions; moving a DOM subtree must not reset its state.

## 2. Options

| Option                                                 | Granularity                                      | Embeddable               | Deps           | Toolchain         | Risk                                                      |
| ------------------------------------------------------ | ------------------------------------------------ | ------------------------ | -------------- | ----------------- | --------------------------------------------------------- |
| A. Own fine-grained JSX runtime + light-DOM elements   | Per binding (no diff)                            | Yes (custom elements)    | rxjs only      | `jsxImportSource` | Own small framework to maintain; "runs once" mental model |
| B. Preact (or React) components inside custom elements | Per component (VDOM re-render of the subscriber) | Yes (a root per element) | preact (~5 kB) | `react-jsx`       | Two runtimes to test if both React and Preact are targets |
| C. Lit                                                 | Per template part                                | Yes                      | lit (~6 kB)    | none              | Tagged templates, not TSX                                 |
| D. Solid                                               | Per binding                                      | Via solid-element        | solid          | Babel plugin      | Breaks requirement 5                                      |
| E. React only (as today, with view models)             | Per component                                    | Only into React hosts    | react (~45 kB) | `react-jsx`       | Keeps React a hard dependency of every consumer           |

C and D fail requirements 2 and 5. E fails 3 and keeps the heaviest dependency. The real choice is A or B.

## 3. Recommendation: A for the built-in UI, framework-neutral view contract for contributions

Because the view models are already granular, a virtual DOM adds work without adding value: every `ReadonlyState` can
bind directly to the text node, attribute, or child range that shows it. A runtime that does exactly that is small, has
no dependency beyond rxjs, and produces real DOM nodes that any host can adopt. As built, `dom/` is about 1,000 lines
(about 750 without comments and blank lines, the SVG tag list included) plus about 520 lines of tests.

The built-in UI uses this runtime (`prism-ui/dom`). Contributions do **not** have to: the shell's view contract is a
mount function (§6), so a contribution can render with React, Preact, Lit, or nothing.

Gate: before committing, a Storybook spike builds three representative pieces in both A and B (Preact): the PD form
(group, mapped, conditioned, color list), the state tree with 5,000 cells and a burst of updates, and the command menu
with keyboard navigation. Compare LOC, update cost (DOM mutations and time per burst), and authoring friction. If A is
not clearly simpler, fall back to B with the same element and contribution contracts; the model layer is unaffected
either way.

**Decided: option A** ([../spike.md](../spike.md)). The state tree at 5,025 cells renders without blocking the page once
`For` slices long lists (12 ms slices, continued with Mol*'s `setImmediate`; off in jsdom, `timeSlicing`), rows skip
off-screen layout (`content-visibility`), sources hand out per-node records (a row follows its record by id; an
immutable tree rebuilds subtrees, `For`'s identity rule, §4.3), and long lists are virtualized
(`components/virtual-list`). None of that is the renderer's to remove; the Preact side was not built.

## 4. The `prism-ui/dom` runtime

### 4.1 JSX

Implemented in `src/dom/` (`scope`, `render`, `jsx-runtime`, `jsx-types`, `control-flow`, `element`, `view`).

**Every `.tsx` file names its runtime with a pragma**: `/** @jsxImportSource @molstar/prism-ui/dom */` (React files:
`/** @jsxImportSource react */`). The tsconfig sets the same default for editors and `tsc`, but bundlers do not all read
it (Vite 8 compiled with React's runtime; Jest's esbuild transform is global), and the package publishes `src/` for
consumers' bundlers, which apply their own JSX settings. `./dom/jsx-runtime` and `./dom/jsx-dev-runtime` resolve through
the package's `./*` export.

`jsx(type, props)` (`createElement` in `dom/render`) creates real nodes; the JSX `key` is ignored (lists key through
`For`, §4.3):

- **Strings** create elements; SVG-only tag names (`svg`, `path`, `g`, …) are created in the SVG namespace (children are
  created before their parent, so the namespace comes from the tag, not the parent).
- **Functions** are components, called **once** with props, in the current render scope (§4.2); they return nodes. There
  is no re-render. A component may return data for its parent component instead of DOM (marked with `Passthrough`, as
  `Match` does for `Switch`). Some components are plain functions in `.ts` files that build their element with DOM calls
  (`Icon`, `Markdown`); they compose the same way.
- **Props and children may be reactive** (`MaybeReactive<T>` in `dom/scope`): an `Observable`, a `ReadonlyState`, or a
  writable `observableState` (its readonly view is used, so local state works directly: `<div hidden={open}>`). A
  reactive prop subscribes and writes the attribute/property on each value. A prop whose value is `undefined` is skipped
  (no binding).
- **Children**: nodes, text, numbers, arrays, `null`/`undefined`/booleans (nothing), functions (called at once), and
  streams. A stream child is a range between two comment anchors (`insertRegion`): text values update one text node in
  place; any other value replaces the range's content, rendered in a fresh child scope (the previous one is disposed),
  and a function value is called inside that scope.
- **Properties vs attributes**: SVG elements and hyphenated names (`data-*`, `aria-*`, on any element) get attributes;
  on a custom element (a tag with a hyphen) every other name is a property (models, lists, callbacks); on other elements
  a name that exists on the element is a property, except `list`, `form`, `download`, `role`, `width`, `height`, which
  are always attributes. `attr:` and `prop:` force either (`attr:popover="manual"`). Attribute values: `null`,
  `undefined` and `false` remove it, `true` sets it empty, and booleans on `aria-*` become `"true"`/`"false"`. A
  property set to `null`/`undefined` is set to `''` and its attribute removed.
- **`select`**: its `value` is applied after its options are rendered (a value set before the options exist does not
  stick).
- **Events**: `onClick`, `onKeyDown`, … (camel case: `on`, a capital letter, and a function; typed from the DOM event
  map with `currentTarget` typed as the element) and `on:any-name` for exact names. Handlers are added once and are not
  reactive; they go away with the element.
- **Classes and styles**: `class` and `classList` combine (the attribute is the base plus the flags that are on);
  `style` takes a (reactive) string or an object of (reactive) values, camelCase or kebab-case or `--msp-*` custom
  properties; numbers other than 0 are pixels except for unitless properties (and custom properties, which get the bare
  number).
- `ref` receives the element after its props and children are set. It, or the element a JSX expression returns, combined
  with `bind` (§4.2) is the escape hatch for DOM that JSX does not express (`views/shell/canvas-chrome` swaps a button's
  icon face through `ref`; `views/canvas/hover-label` keeps its label element and places it each frame).

```tsx
function ComponentRow({ model }: { model: ComponentVM }) {
  return (
    <div class={classes.row} classList={{ [classes.hidden]: model.isHidden }}>
      <span class={classes.label}>{model.label}</span>
      <Button
        variant="ghost"
        size="sm"
        icon={<Icon name="eye" />}
        title={model.isHidden.pipe(map((h) => (h ? 'Show' : 'Hide')))}
        pressed={model.isHidden}
        onClick={() => model.toggleVisibility()}
      />
    </div>
  );
}
```

`model.label` is a `ReadonlyState<string>`; the text node updates when it changes and nothing else does. `classes` is a
CSS Module's generated class map ([styling.md](styling.md)); `Button` and `Icon` are in `components/`.

### 4.2 Ownership

Every binding is owned by the current **render scope** (`Scope` in `dom/scope`). Disposing a scope disposes its child
scopes first (newest first), then runs its own cleanups in reverse order; a failing cleanup is logged and the others
still run. Scopes are opened by a `render` root, an `MspElement`'s rendering, the content of each stream region (a
reactive child, a `Show`/`Switch`/`Dynamic` branch), and each `For` entry and fallback; replacing or removing that
content disposes its scope. **A component does not open a scope of its own**: its bindings and cleanups belong to the
nearest enclosing one, so they live as long as the region, list entry or root the component was rendered in.

The scope API, all in `dom/scope`:

- `onCleanup(fn)`: runs `fn` when the current scope is disposed (at once if it already is). Used for imperative work in
  component bodies: document listeners, observers, timers, `installTooltips`.
- `bind(source, fn)`: calls `fn` with a plain value once, or with each value of a reactive one while the scope lives.
  The JSX runtime is built on it; views use it with `ref` for imperative updates.
- `createRoot(fn)`, `runInScope(scope, fn)`, `getScope()`, `requireScope()`, `toStream(value)`.

A reactive binding created outside any scope throws instead of leaking its subscription.

### 4.3 Control flow

All in `dom/control-flow`; each renders into a range of its own and disposes what it leaves.

- `<Show when={value} fallback={() => ...}>{(v) => ...}</Show>`: renders the children with the value while it is truthy,
  else the fallback. The branch is rendered again, in a fresh scope, whenever the value changes identity; for a flag,
  only when it flips. A view that passes an object either wants one rendering per value (the state tree's inspector is
  shown under a `block\nref` string key, so it is rebuilt for another cell) or narrows the stream first (a `map` to a
  flag, or `distinctUntilChanged` with its own comparison, as the dock content does for its focused section).
- `<For each={list} key={(item) => item.id} fallback={() => ...}>{(item, index) => ...}</For>`: keyed. An entry keeps
  its nodes and scope while its key stays in the list **and** its item is the same value; the key defaults to the item
  itself. When the key stays but the item is a different value, the entry is disposed and rendered again: the children
  got the item as a plain value, so a kept entry would show the old record. `index` is a `ReadonlyState<number>` that
  follows the entry's position. Entries are placed from the end and one already in place is not touched; an entry of
  several nodes is kept between two anchors. A duplicate key warns once and later items with it are not rendered.
- `<Switch fallback={...}>` with `<Match when={...}>{(v) => ...}</Match>` cases: the first truthy case renders, again
  when the case or its value changes. Built and tested; no view uses it yet.
- `<Dynamic component={state} {...props} />`: renders the current component with the other props as they are, again when
  the component changes. Built and tested; no view uses it yet (dock and card content mount through `ViewHost`, §6).

**Views follow records by id.** Most view-model lists are keyed trees of stable records (scene nodes by id with reactive
`children`, sections, tabs, providers by id; [view-models.md](view-models.md) §4) and use `key` directly. The state tree
keys by ref too, but its mock replaces records (§3). Where a list's records are replaced on every change (an agent
turn's runs, History items), keying by id would still render a row anew on each change (a spinner restarts, a fold
closes). Such views run `For` over the ids (a stream of id lists with `distinctUntilChanged` on their contents) and give
each row a stream of its record, found by id in the list (`RunList` and `HistoryTab` in the agent package's `views/`).
Lists edited by position (a PD object list's items, a line graph's points) iterate indices the same way, each row
reading its item by index.

### 4.4 Rules that keep it correct

- **Read states, not values, in render.** `{model.label}` binds; `{model.label.value}` renders once. Not yet: a lint
  rule (Biome plugin or a small script in CI) that flags `.value` inside JSX expressions in `prism-ui/**/*.tsx`.
- **Derive in the model, not in JSX.** Inline `.pipe(map(...))` in JSX is acceptable for trivial formatting because
  components run once (no resubscription per render, unlike React); anything shared or expensive is a model field.
- **Keep identities stable where content should stay.** `Show`, `Switch`, `Dynamic` and `For` render again when a value
  changes identity; narrow the stream (`map`, `distinctUntilChanged`) or follow records by id (§4.3) when the DOM should
  stay.
- **Combine, then debounce.** A row that depends on several streams combines them (`combineLatest`); where bursts are
  expected (state-tree rebuilds, selection changes) a short debounce makes one change render once (plan.md, Phase 4).
  Not yet applied in the mocks.
- **No component-local mutable state that should survive moves.** Put it in a model.

## 5. Custom elements

Custom elements are used at **island** granularity: things a host could embed or a dock could move (shell, dock region,
view host, viewport, toolbar, command menu, PD form, state tree, sequence, PAE plot, snapshot stepper). Inside an
island, plain function components do the work; elements are not the re-render boundary (bindings are).

State (mocks): the element base is built and tested (`dom/element`), but no island element exists yet. Outside tests,
the only `MspElement` subclass is `MspStoryRoot` (`testing/story-root`), which hosts a story's JSX so that Storybook
removing it disposes the bindings. The shell, docks, views and components are function components rendered inside it.
The element below shows the intended shape of an island, not code that exists:

```tsx
export class MspStateTree extends MspElement<StateTreeModel> {
  static tag = 'msp-state-tree';
  protected render(model: StateTreeModel) {
    return <StateTreeTab model={model} />;
  }
}
```

- `MspElement<M>` has a `model` property. Rendering starts when the element is connected and has a model, into its own
  light DOM, once; setting another model disposes the rendering and renders again. A model set before the element was
  defined (an own property on the not-yet-upgraded element) is taken over.
- `disconnectedCallback` schedules disposal of the rendering on a microtask and cancels it if the element is reconnected
  first, so moving an element keeps its DOM and bindings.
- On connect an element adopts the registered stylesheets into its document; one that is not inside another prism-ui
  root gets `.msp-root`, which scopes the reset and the theme tokens ([styling.md](styling.md)).
- An element **does not own** its model by default (the shell or the parent view model does: one owner per model). An
  element that creates its own model (`createModel()`; planned: `<msp-viewer src=…>` for plain pages) mounts it while
  rendered and unmounts it when the rendering is disposed.
- The base class extends `HTMLElement` only where it exists, so modules that define elements import in Node and workers.
- **Light DOM, no shadow root.** Styles come from CSS Modules and `--msp-*` tokens ([styling.md](styling.md)); hosts can
  theme through the tokens, and a single stylesheet serves the whole tree. Leakage from host pages is handled by a reset
  scoped to `.msp-root` inside `@layer msp.reset`. Shadow DOM can be added later per element if a host needs hard
  isolation.
- Planned: **attributes only on top-level embeddables** (`src`, `format`, `theme`); everything else is a property.
- Tags are prefixed `msp-` and registered by an explicit `defineElements([...])` call (no import side effects),
  consistent with the v6 "explicit composition" rule; a tag already registered is left as it is.
- JSX typing: an element adds itself to the `MspElements` interface of `dom/jsx-types` by module augmentation
  (`interface MspElements { 'msp-x': MspElementProps<MspX> }`), which `JSX.IntrinsicElements` extends, so
  `<msp-x model={vm} />` is type-checked; JSX sets `model` as a property (§4.1).

Surviving moves: an element keeps its DOM when moved (above). Views in the mock shell are not moved as DOM: a dock or
card mounts its view through `ViewHost` (§6) and mounts it again when it is shown anew, so a view's lasting state lives
in models (§4.4).

## 6. Framework-neutral view contract

Views contributed to the shell (tabs, sections, overlays, a command's custom UI) implement (`dom/view`):

```ts
interface ViewRenderer<Ctx> {
  /** Render into `host`; return a disposer. Called when the view becomes visible, disposed when it is closed. */
  mount(host: HTMLElement, ctx: Ctx): () => void;
}
```

Helpers produce renderers for each technology:

- `domView(Component)` (`dom/view`): the built-in runtime; renders `Component({ ctx })` with `render` after what `host`
  holds, and the disposer disposes its scope and removes its nodes.
- `reactView(Component, { createRoot })` (`react/react-view`): appends a `display: contents` container to `host`,
  creates a root on it and renders `<Component ctx={ctx} />`; the disposer unmounts the root and removes the container.
  The caller passes `createRoot` (from `react-dom/client` or `preact/compat/client`), as `renderReact18` is injected
  today, so prism-ui has no React dependency. Not yet: a `ViewModelContext` provider around the root.
- `elementView('my-tag', property = 'ctx')` (`dom/view`): creates the element, sets the context on that property and
  appends it; the disposer removes it (the element disposes its rendering on disconnect, §5).

The shell mounts a renderer with `ViewHost` (`views/shell/view-host`): a `display: contents` host whose disposer is
registered with `onCleanup`, so the view is disposed with the region that shows it. Dock tabs, sections (`SectionList`)
and canvas cards use it.

State (mocks): every contributed view in the mocks is a `domView`; `reactView` and `elementView` are covered by jest
tests only (`dom/_test/element.test`, `react/_test/react-view.test`).

## 7. React and Preact on top of the view models

`react/` ships the any-ui hooks, one module each (the package exports `./*`; there is no index):

- `useObservable(source, { defaultValue?, compare? })`: the current value (a `ReadonlyState`'s at once).
- `useObservableProp(source, selector, options?)`: a slice; re-renders only when the selection changes (shallow
  comparison by default). Both are built on `useObservableSelection`.
- `useReactiveModel(model)`: mounts the model in an effect and unmounts it on cleanup; safe under Strict Mode replays,
  and it never retires the model (its creator does).

They use `useSyncExternalStore`, which React 18+ and `preact/compat` both provide, so one build serves both: a Preact
app aliases `react` to `preact/compat` as usual. A source's identity is its lifetime: a new observable on every render
(an inline `.pipe(...)`) starts over each time, and the hooks warn when the source changes on renders in a row, so
sources are memoized. Peer dependency on `react >= 18` is optional; nothing else in prism-ui imports React. The `react/`
modules are plain `.ts` using `createElement`; React test files name `@jsxImportSource react` (§4.1).

Not yet: `useViewModel()` and its context (a React view gets its `ctx` as a prop, §6). Preact is not tested; the tests
run React with `react-dom`.

The React adapter replaces today's `useBehavior` for new code; the existing `@molstar/plugin-ui` keeps its own.

Hosting prism-ui elements inside React: React 19 sets properties on custom elements natively
(`<msp-viewport model={vm} />`). From React 18 on (and with `preact/compat`), `MspElementHost` (`react/element-host`)
renders `<MspElementHost tag="msp-…" properties={{ model }} />` and sets each property on the element in a layout
effect.

## 8. Accessibility and primitives

Build on platform primitives instead of a component kit: `popover` attribute and `<dialog>` for menus and dialogs,
`inert` for covered content, `:focus-visible`, roving tabindex for trees and toolbars, combobox ARIA for the command
menu (input keeps focus; rows are `role="option"`, as in stories-untold). A small set of behaviour helpers (focus trap,
roving focus, typeahead, positioning via CSS anchor positioning with a JS fallback) lives in `prism-ui/dom/behaviors`.

State (mocks):

- Floating lists (`Select`'s options with `over`, in `components/select`; in forms they open inline under the row) are
  manual popovers (`attr:popover="manual"`) in the top layer, so no panel's clipping or fade cuts them, placed under (or
  above) their anchor in JS and following it on scroll and resize. No CSS anchor positioning yet.
- The command menu is a combobox: the input keeps focus (`role="combobox"`, `aria-controls`), arrows move the highlight,
  and rows are `role="option"` with `aria-selected`.
- Trees (`views/structure/scene-tree`, `views/state-tree`) are `role="tree"` with `treeitem`s carrying `aria-expanded`
  and `aria-selected`; no roving tabindex yet.
- Behaviours live with the components that need them, not in `dom/behaviors` (which does not exist): `installTooltips`
  (`components/tooltip`, one reused tooltip for `data-tip` elements, installed on the shell root) and `dragToReorder`
  (`components/reorder`).
- Tooltips: the rail, the canvas stack, cut-off labels (`data-tip-overflow`) and tree and form rows use the app's
  tooltip and carry no native `title` (icon buttons keep their `aria-label`); a canvas button whose menu opens on hover
  has none. `Button`'s `title` is still a native title inside panels.
- `dragToReorder` follows the document, not the handle: a release anywhere, a cancelled pointer or the window losing
  focus always ends the drag and clears its marks (`data-dragging`, `data-drop`).
- Not yet: `<dialog>`, `inert`, focus trap, typeahead.
