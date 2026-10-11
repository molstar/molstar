# Rendering: web components, JSX, React and Preact

How ui-next turns view models into DOM. The model layer ([view-models.md](view-models.md)) is renderer-free; this
document chooses the renderer for the built-in UI and defines how other renderers plug in.

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
bind directly to the text node, attribute, or child range that shows it. A runtime that does exactly that is small
(estimate 600–900 LOC with tests), has no dependency beyond rxjs, and produces real DOM nodes that any host can adopt.

The built-in UI uses this runtime (`ui-next/dom`). Contributions do **not** have to: the shell's view contract is a
mount function (§6), so a contribution can render with React, Preact, Lit, or nothing.

Gate: before committing, a Storybook spike builds three representative pieces in both A and B (Preact): the PD form
(group, mapped, conditioned, color list), the state tree with 5,000 cells and a burst of updates, and the command
command menu with keyboard navigation. Compare LOC, update cost (DOM mutations and time per burst), and authoring
friction. If A is not clearly simpler, fall back to B with the same element and contribution contracts; the model layer
is unaffected either way.

## 4. The `ui-next/dom` runtime

### 4.1 JSX

Implemented in `src/dom/` (`scope`, `render`, `jsx-runtime`, `jsx-types`, `control-flow`, `element`, `view`).

**Every `.tsx` file names its runtime with a pragma**: `/** @jsxImportSource @molstar/plugin-ui-next/dom */` (React
files: `/** @jsxImportSource react */`). The tsconfig sets the same default for editors and `tsc`, but bundlers do not
all read it (Vite 8 compiled with React's runtime; Jest's esbuild transform is global), and the package publishes `src/`
for consumers' bundlers, which apply their own JSX settings. The package exports `./dom/jsx-runtime` and
`./dom/jsx-dev-runtime`.

`jsx(type, props)` creates real nodes:

- **Strings** create elements; SVG-only tag names (`svg`, `path`, `g`, …) are created in the SVG namespace (children are
  created before their parent, so the namespace comes from the tag, not the parent).
- **Functions** are components, called **once** with props; they return nodes. There is no re-render. A component may
  return data for its parent component instead of DOM (marked with `Passthrough`, as `Match` does for `Switch`).
- **Props and children may be reactive**: an `Observable`, a `ReadonlyState`, or a writable `observableState` (its
  readonly view is used, so local state works directly: `<div hidden={open}>`). A reactive prop subscribes and writes
  the attribute/property on each value; a reactive child is a text node updated in place, or a range between two comment
  anchors whose content is rendered in a fresh scope (a function value renders lazily in that scope).
- **Properties vs attributes**: custom elements (`msp-*`) always get properties (models, arrays, callbacks); SVG and
  hyphenated names (`data-*`, `aria-*`) get attributes; other names that exist on the element are properties (with a few
  read-only exceptions). `attr:` and `prop:` force either.
- **Events**: `onClick`, `onKeyDown`, … (camel case, typed from the DOM event map) and `on:any-name` for exact names.
- **Classes and styles**: `class` and `classList` combine (the attribute is the base plus the flags that are on);
  `style` takes a string or an object of (reactive) values, camelCase or kebab-case or `--msp-*` custom properties;
  numbers are pixels except for unitless properties.
- `ref` receives the element after its props are set.

```tsx
function ComponentRow({ model }: { model: ComponentVM }) {
  return (
    <div class={styles.row} classList={{ [styles.hidden]: model.isHidden }}>
      <span class={styles.label}>{model.label}</span>
      <IconButton icon={model.isHidden.pipe(map(h => (h ? 'visibility-off' : 'visibility')))}
                  onClick={() => model.toggleVisibility()} />
    </div>
  );
}
```

`model.label` is a `ReadonlyState<string>`; the text node updates when it changes and nothing else does.

### 4.2 Ownership

Every binding is owned by the current **render scope** (a disposable group, the same primitive as `ReactiveModel`'s
mount scope). Creating a component or a list item opens a child scope; removing the node disposes it. Component bodies
can call `onCleanup(fn)` and `useScope()` for imperative work (a canvas observer, a `ResizeObserver`).

### 4.3 Control flow

- `<Show when={state} fallback={...}>{() => ...}</Show>`: mounts and disposes a branch.
- `<For each={state} key={item => item.ref}>{item => <Row model={item} />}</For>`: keyed reconciliation; an item's nodes
  and scope live as long as its key. View models are keyed by state-tree ref, so rows keep identity across updates (see
  [view-models.md](view-models.md) §4).
- `<Switch>`/`<Match>` for tagged unions (`status.kind`).
- `<Dynamic component={state} />` for slots whose content changes (the active view in a dock).

### 4.4 Rules that keep it correct

- **Read states, not values, in render.** `{model.label}` binds; `{model.label.value}` renders once. A lint rule (Biome
  plugin or a small script in CI) flags `.value` inside JSX expressions in `ui-next/**/*.tsx`.
- **Derive in the model, not in JSX.** Inline `.pipe(map(...))` in JSX is acceptable for trivial formatting because
  components run once (no resubscription per render, unlike React); anything shared or expensive is a model field.
- **No component-local mutable state that should survive moves.** Put it in a model.

## 5. Custom elements

Custom elements are used at **island** granularity: things a host could embed or a dock could move (shell, dock region,
view host, viewport, toolbar, command menu, PD form, state tree, sequence, PAE plot, snapshot stepper). Inside an
island, plain function components do the work; elements are not the re-render boundary (bindings are).

```ts
export class MspStateTree extends MspElement<StateTreeVM> {
  static tag = 'msp-state-tree';
  protected render(model: StateTreeVM) {
    return <For each={model.roots} key={n => n.ref}>{n => <TreeNode model={n} />}</For>;
  }
}
```

- `MspElement<M>` has a `model` property. `connectedCallback` renders into the element (light DOM) once a model is set;
  `disconnectedCallback` schedules disposal of the render scope on a microtask and cancels it if the element is
  reconnected in the same task, so moving an element (docking) keeps its DOM and bindings.
- An element **does not own** its model by default (the shell or the parent view model does: one owner per model). An
  element that creates its own model (`<msp-viewer src=…>` for plain pages) mounts and retires it.
- **Light DOM, no shadow root.** Styles come from CSS Modules and `--msp-*` tokens ([styling.md](styling.md)); hosts can
  theme through the tokens, and a single stylesheet serves the whole tree. Leakage from host pages is handled by a reset
  inside `@layer msp`, as `.msp-plugin` does today. Shadow DOM can be added later per element if a host needs hard
  isolation.
- **Attributes only on top-level embeddables** (`src`, `format`, `theme`); everything else is a property.
- Tags are prefixed `msp-` and registered by an explicit `defineElements([...])` call in the entry point (no import side
  effects), consistent with the v6 "explicit composition" rule.
- JSX typing: ui-next augments its own `JSX.IntrinsicElements` with each element's props, so
  `<msp-state-tree model={vm} />` is type-checked.

## 6. Framework-neutral view contract

Views contributed to the shell (panels, sections, dialogs, overlays) implement:

```ts
interface ViewRenderer<Ctx> {
  /** Render into `host`; return a disposer. Called when the view becomes visible, disposed when it is closed. */
  mount(host: HTMLElement, ctx: Ctx): () => void;
}
```

Helpers produce renderers for each technology:

- `domView(Component)`: the built-in runtime.
- `reactView(Component, { createRoot })` in `ui-next/react`: creates a root on `host`, wraps it in a `ViewModelContext`
  provider, unmounts on dispose. The caller passes `createRoot` (from `react-dom/client` or `preact/compat/client`), as
  `renderReact18` is injected today, so ui-next has no React dependency.
- `elementView('my-tag')`: any custom element; the shell sets `model`/`ctx` properties on it.

## 7. React and Preact on top of the view models

`ui-next/react` ships the any-ui hooks (`useObservable`, `useObservableProp`, `useReactiveModel`) plus `useViewModel()`
(context). They use `useSyncExternalStore`, which React 18+ and `preact/compat` both provide, so one build serves both:
a Preact app aliases `react` to `preact/compat` as usual. Peer dependency on `react >= 18` is optional; nothing else in
ui-next imports React.

The React adapter replaces today's `useBehavior` for new code; the existing `@molstar/plugin-ui` keeps its own.

Hosting ui-next elements inside React: React 19 sets properties on custom elements natively
(`<msp-viewport model={vm} />`). For React 18 the adapter exports thin wrappers that set properties through a ref.

## 8. Accessibility and primitives

Build on platform primitives instead of a component kit: `popover` attribute and `<dialog>` for menus and dialogs,
`inert` for covered content, `:focus-visible`, roving tabindex for trees and toolbars, combobox ARIA for the command
menu (input keeps focus; rows are `role="option"`, as in stories-untold). A small set of behaviour helpers (focus trap,
roving focus, typeahead, positioning via CSS anchor positioning with a JS fallback) lives in `ui-next/dom/behaviors`.
