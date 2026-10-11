# stories-untold: patterns to adopt

stories-untold (`/Users/david/projects/experiments/stories-untold`) vendors anyfold's `any-ui` and `any-ui-declarative`
and builds its frontend on them. This note records what the frontend does, what ui-next takes, and what it leaves. Paths
are relative to that repository.

## 1. Stack

React 19 function components, TSX, Vite, CSS Modules, rxjs 7. **No web components, no Preact, no custom JSX runtime**:
nothing in stories-untold validates the web-component direction; only the model layer (React-free) carries over.

Layering: `vendor/anystar/any-ui` (state, hooks, lib, commands, declarative state) → `lib/stories/stories-ui` (dumb
components) → `lib/stories/stories-frontend` (screens + models) → `apps/web` (thin shell). `lib/stories/any-molstar`
wraps Mol* as "a canvas plus rxjs streams"; Mol*'s own UI is never used.

## 2. The model layer (`any-ui/src/state`, ~650 LOC)

- `observableState<T>(initial, { compare })`: a `BehaviorSubject` with equality-suppressed `set` (shallow for plain
  records/arrays), `update(partial)` for records, and re-entrant sets queued so every subscriber sees values in order.
  Consumers get `ReadonlyState<T> = Observable<T> & { readonly value: T }`; writers stay private.
- `ReactiveModel`: `mount()`/`unmount()` (repeatable) and `retire()` (final); lifecycle methods cannot be overridden
  (`onMount`/`onUnmount`/`onRetire` hooks instead). Inside `onMount`, `this.subscribe(obs, fn)` and
  `this.registerDispose(fn)` are owned by the mount scope and disposed LIFO on unmount. One owner mounts a model; a
  second mount warns.
- `combineBehaviors({ sources, project })`: a derived `ReadonlyState` guarded against throwing projections.
- `createDisposableGroup()`: LIFO cleanup.

Conventions in `stories-frontend/src/models`:

- Public `readonly x = observableState(...)` fields; views read `model.x.asReadonly`. Logic in methods.
- Derived state is recomputed into other states in `this.subscribe(...)` blocks.
- Parents construct children in the constructor and mount/unmount them in `onMount`. Dynamic children live in an
  `observableState<Child | null>` swapped from a subscription.
- Dependencies are injected through the constructor; upward calls go through a plain `handlers` object.
- Event-like things are `Subject`s (no replay), state-like things are `BehaviorSubject`s. An outcome held as replaying
  state caused a cross-screen bug; hence "outcomes are events".
- Sources never error: an error is a value in the state.
- Stream identity is subscription lifetime: streams are created once and cached (an inline `.pipe()` in render
  resubscribes each time; `useObservable` warns).
- Models are tested without a renderer: construct, `mount()`, feed stubs.

Known rough edges the authors list: per-mount `AbortSignal` missing (hand-rolled sequence guards), module-global
singletons (`overlay`, `toasts`), models holding DOM handles, `useObservableProp` unused so some parts re-render per
streamed chunk.

## 3. Optimistic transactions (`any-ui/src/state/transaction.ts`)

`createUiTransactionCoordinator()` issues transactions (`acknowledge`, `reject`, `reconcile`);
`createOptimisticProjection(base, coordinator)` folds active contributions (`apply`, `isAbsorbed`) over the
authoritative base, drops absorbed ones on `setBase`, and reports divergences after emitting. The only consumer is
`StoryView.edit`: validate locally, `begin()`, `contribute()`, send, `acknowledge()`/`reject()`, `setBase()` on every
authoritative snapshot. Operations are idempotent by desired state, so they re-apply cleanly on a new base.

## 4. Commands (`any-ui/src/commands`, ~1,500 LOC)

- `CommandDefinition<Ctx>`: `id`, `title`, `category`, `alias` (words that start the command with arguments), `keys`
  (display) and `shortcut` (`{ key, alt?, shift?, mod? }`), `when(ctx) → true | reason` (blocked commands stay listed,
  greyed, with the reason), `form(ctx, values)`, `reshapes`, `presets(ctx)` (search-only "ways in"), `defaults(ctx)`,
  `preview(values, ctx)`, `check`, `run(values, ctx, signal)`, `assist` (LLM fills the form).
- `CommandManager<Ctx>`: `register(defs) → unregister` (screens own their commands by mount scope; unregistering cancels
  their run and closes their form), `query(text) → { recent, matches, direct }` (fuzzy subsequence scoring; alias + tail
  opens a form prefilled from the tail), `open(id, tail?, preset?)`, one run at a time with explicit cancel, `outcomes`
  as a non-replaying event stream, `handleKey` (IME/repeat/text-target hygiene, `e.code` fallback), injected recents
  store.
- `CommandInvocation`: one open form as a state machine (`editing → running → done/failed/cancelled`), field focus,
  typed text, completions, tail parsing (`tail.ts`), preview sink keyed on values + `previewDeps(ctx)`.
- **Formats** (`CommandFormats<Ctx>`): named input types with `hint`, `check(text) → value | message`, and
  `complete(text) → options`, e.g. `su.target`, `su.color`.
- Forms are built as `any-ui-declarative` definitions (zod-validated JSON) and read back as a flat field list. The
  authors' own audit: the declarative stack (~3,400 LOC incl. generated schema and the zod pin) is used only "as a value
  store and a one-action runner"; it is kept for agent-authored UI. A flat field schema was considered (~240 LOC).
- Command menu: `CommandPromptModel` (menu derived from query + context + registrations; categories in fixed order,
  recents first, disabled reasons, fallback row) + dumb `CommandMenu` (combobox ARIA; input keeps focus) + memoized
  views. The command menu is a "command mode" inside the chat prompt, not a modal.
- Commands call the same action functions as the agent's tools, so UI and agent agree.

## 5. Views

Each component calls `useObservable(model.x.asReadonly)` for exactly what it renders; granularity comes from many small
`function Part({ model })` components, `memo`, and equality-suppressed state. `useObservable` is a
`useSyncExternalStore` store that reads `.value` synchronously (first render has the value) and skips the replay. Dumb
components (`stories-ui`) take plain props and know nothing about models or Mol*. One React context (the command host).

## 6. Mol* integration

`MolstarView` (`lib/stories/any-molstar/src/view.ts`): `new PluginContext(spec)` with no Mol* UI. Behaviors and events
become cached rxjs streams (`clicked`, `hovered`, `selection`, `currentKey`, `loading`, …); mutations are methods
serialized through a promise queue; parts are addressed declaratively (MVS refs, schema rows) and resolved to loci
inside. `StageModel` mirrors app state and Mol* state both ways with explicit feedback guards. Only selection, hover,
click and snapshot key are exposed; state-tree cells, transforms and managers are not projected as view models. So the
Mol* view-model catalogue in ui-next is new design work, not a port.

## 7. Shell and styling

Fixed layout (`Split` over `react-resizable-panels` + an absolutely positioned stage with slots). No docking or panel
registry. Overlays cover a still-mounted workspace (`inert`) so the canvas and drafts survive; the render loop pauses
while covered. CSS Modules + `--su-*` design tokens; light/dark via `:root[data-theme]`; `data-ink` for text drawn over
the canvas, derived from the scene background.

## 8. Storybook

Storybook 10 with `@storybook/react-vite`, docs and a11y addons. Three kinds of stories: dumb components with inline
data, pattern mocks, and real screens over the real model graph with only the data client swapped for an in-memory one
and the canvas replaced by a slot. `play` functions test interactions (including optimistic refusal).

## 9. Take / adapt / leave

| Item                                                   | Decision                                                                                                  |
| ------------------------------------------------------ | --------------------------------------------------------------------------------------------------------- |
| `state/` (observable state, ReactiveModel, combine)    | Copy as new Mol* code (Mol* header); reformat to repo style.                                              |
| `lib/` (disposable group, shallow equality)            | Take.                                                                                                     |
| `transaction.ts`                                       | Take; use for async state-tree edits (see [../design/view-models.md](../design/view-models.md)).          |
| `hooks/` (`useObservable`, `useReactiveModel`)         | Take into the React/Preact adapter, not the core.                                                         |
| `commands/` manager, match, keys, recents              | Take; adapt `Ctx`, outcomes, keybinding scopes.                                                           |
| `commands/` form + invocation + tail                   | Adapt: forms become `ParamDefinition` params; formats become PD-aware argument types.                     |
| `declarative/` + `any-ui-declarative`                  | Leave (zod/yaml dependency, agent-authored UI is not a goal). Revisit if agent-authored forms become one. |
| `assist` (LLM fill)                                    | Keep the port in the type, no implementation.                                                             |
| Mol* facade ideas (cached streams, queue, declarative) | Adapt into view models over managers.                                                                     |
| Layout, routing, prompt-dock UX                        | Leave; Mol* needs a dockable shell and a modal command menu.                                              |
| CSS Modules + prefixed tokens + `data-theme`           | Take with `--msp-*`.                                                                                      |
| Storybook "real model, fake data" approach             | Take; mock at the view-model source boundary.                                                             |
| Vendoring process (`CHANGES.md`, in-code markers)      | Leave: the copied code is new Mol* code; stories-untold is an inspiration, not an upstream.               |
