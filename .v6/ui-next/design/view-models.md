# View models

ui-next's state layer: renderer-free classes built on the any-ui primitives (`observableState`, `ReadonlyState`,
`ReactiveModel`, transactions) that turn the plugin's managers, events, and state tree into UI-shaped, granular state.
Every view (the built-in DOM runtime, React, Preact, a test) reads the same view models.

## 1. Short answer to "do view models wrap manager behaviors?"

Yes, in two layers:

1. **Sources** turn plugin APIs into clean `ReadonlyState`s and event streams, one per plugin area, and hide the gaps
   listed in [../analysis/plugin-state.md](../analysis/plugin-state.md) §2 (void events, in-place mutation, late
   binding, silent writes, non-reactive registries). Sources know Mol*; they know nothing about panels.
2. **View models** are shaped by what a view shows: keyed trees of small objects (one per structure, component,
   representation, cell), each field its own `ReadonlyState`, plus methods that perform actions. View models depend on
   source **interfaces**, so Storybook and tests replace sources, not the plugin.

```
PluginContext ──► sources (Mol*-aware, UI-agnostic) ──► view models (UI-shaped) ──► views (dom / React / Preact)
      ▲                                                        │
      └──────────────────────── actions ◄──────────────────────┘ (also called by commands)
```

## 2. Primitives

From any-ui, unchanged in behaviour (see [architecture.md](architecture.md) §4 for the import):

- `observableState<T>(initial, { compare })`, `ReadonlyState<T>` (`Observable<T> & { readonly value: T }`).
- `ReactiveModel` with `onMount`/`onUnmount`/`onRetire`, mount-scoped `subscribe`/`registerDispose`, one owner.
- `combineBehaviors`, `createDisposableGroup`, `sameShallow`/`shallowEqual`.
- `createUiTransactionCoordinator`, `createOptimisticProjection`.

Additions ui-next needs:

- **Mount-scope `AbortSignal`** (`this.signal` inside `onMount`), aborted on unmount, instead of hand-rolled sequence
  guards; every async action started by a model takes it.
- **`keyedList(source, key, create, update)`**: maps a `ReadonlyState<readonly T[]>` to a
  `ReadonlyState<readonly VM[]>`, reusing the view model for a key, calling `update(vm, item)` on survivors, retiring
  removed ones. The emitted array is new only when membership or order changes; field changes flow through each view
  model's own states. This is what makes the DOM runtime's `<For>` touch only what changed.
- **`fromEvent(trigger, read, compare)`**: a `ReadonlyState` re-read on a void event (`selection.events.changed` →
  `selection.stats`) with equality suppression. Most sources are built with it.

Conventions:

- Public fields are `ReadonlyState`; writable states are private (`#label = observableState(...)`, exposed as
  `readonly label = this.#label.asReadonly`). Changes go through methods.
- State-like things are behaviors, event-like things are `Subject`s that never replay (command outcomes, "focus
  requested", toasts to show).
- Errors are values (`status: { kind: 'error', message }`), never stream errors.
- Streams are created once per model (fields or cached getters), never per render.
- Dependencies come through the constructor; upward calls through a `handlers` object.

## 3. Sources

One source per area, created by `createPluginSources(plugin)` after `plugin.init()`. Each is a small class (or object)
with `ReadonlyState`s and `Observable`s, disposed with the plugin UI.

| Source                 | Exposes                                                                                                                                                     | Built from                                                                                                       |
| ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| `busy`                 | `isUpdating`, `isAnimating`, `isBusy` (immediate), `isBusySettled` (debounced)                                                                              | `behaviors.state.*`, `plugin.isBusy`                                                                             |
| `stateTree(state)`     | `cell(ref): CellSource` (`label`, `status`, `flags`, `version`, `error`, `children`, `params`) cached by ref; `current`; `history` (`canUndo`, `undoLabel`) | `cell.*`/`object.*`/`changed`/`historyUpdated` events; per-ref diff of `{status, state, transform.version, obj}` |
| `structures`           | `hierarchy` (keyed by ref; stable per ref), `selected` (`trajectories`, `models`, `structures` refs)                                                        | hierarchy `behaviors.selection` + state events; reads getters after events, never `.value` alone                 |
| `volumes`, `particles` | same pattern                                                                                                                                                | their hierarchy managers                                                                                         |
| `selection`            | `stats`, `entries` (per structure: label, element count), `history`, `granularity`, `mode`                                                                  | selection/interactivity void events, `behaviors.interaction.selectionMode`                                       |
| `focus`                | `current`, `history`                                                                                                                                        | focus manager                                                                                                    |
| `measurements`         | `items` (keyed by cell ref, typed), `options`                                                                                                               | measurement `behaviors.state`, deduplicated by ref list                                                          |
| `components`           | `options` (PD values), `optionsParams`                                                                                                                      | `optionsUpdated`                                                                                                 |
| `snapshots`            | `entries`, `current`, `isPlaying`, `params`                                                                                                                 | snapshot `events.changed`, `state.snapshotParams`                                                                |
| `animation`            | `animations`, `current`, `isPlaying`, `params`                                                                                                              | animation events                                                                                                 |
| `tasks`                | `running` (list with progress, one id space, `abort()`)                                                                                                     | accumulated `events.task.*`                                                                                      |
| `log`, `toasts`        | `entries` (capped), `toasts`                                                                                                                                | `events.log`, toast `events.changed`                                                                             |
| `canvas`               | `ready`, `props`, `background`, `reprCount`, `xr`, `camera`                                                                                                 | `canvas3dInitialized`, `settingsUpdated` + `didDraw` re-read, `camera.stateChanged`                              |
| `interaction`          | `hover`, `hoverLabels`, `click` (events, loci kept inside)                                                                                                  | `behaviors.interaction.*`, `behaviors.labels.highlight`                                                          |
| `registry`             | `version`, and derived lists: representation types, color themes, presets, selection queries, formats, state actions, animations                            | wraps `plugin.register` (see below)                                                                              |

**Registry changes.** Until core emits a registry event, `createPluginSources` wraps `plugin.register` on the instance
to bump `registry.version` after each call and its undo; derived lists recompute on `version`. With the core follow-up
([../analysis/plugin-state.md](../analysis/plugin-state.md) §5) the wrapper goes away.

**Identity.** Sources key everything by state-tree ref. A hierarchy rebuild that recreates `StructureRef` objects keeps
the same `StructureSource` instance for the same ref; its fields update only when `transform.version`, `cell.state`, or
`cell.status` change.

**Batching.** Sources update after a state transaction ends (`changed` with `inTransaction: false`) and after animation
stops, as the managers already do, so views do not render intermediate trees.

## 4. View models

### 4.1 Root

`PluginUIModel` (one per plugin UI; `PluginUIContext.model`). The context holds the infrastructure (`shell`, `commands`,
`keybindings`, `contributions`, `sources`; [shell.md](shell.md) §6); the root model holds the area view models built on
it:

```ts
class PluginUIModel extends ReactiveModel {
  constructor(ctx: PluginUIContext);
  readonly scene: SceneModel;            // structures + volumes + particles as one tree
  readonly selection: SelectionModel;
  readonly tasks: TaskListModel;
  readonly notifications: NotificationsModel;  // toasts + command outcomes + log
  readonly viewport: ViewportModel;
  // created lazily by the views that need them:
  stateTree(): StateTreeModel; snapshots(): SnapshotsModel; sequence(): SequenceModel; ...
}
```

### 4.2 Catalogue

| View model                                                                                                                 | Replaces (today)                                                                              | Notes                                                                                                                                                                                                                                                                    |
| -------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `SceneModel` → `StructureNode` → `ComponentNode` → `RepresentationNode`; `VolumeNode`, `ParticlesNode`, `GenericEntryNode` | Source, Components, Volume, Particles, generic entries cards                                  | One tree for everything in the scene. Each node: `label`, `isHidden`, `isSelected` (hierarchy selection), `isFocused`, `status`, `children`, actions (`toggleVisibility`, `remove`, `focus`, `select`, `rename`). Context actions open commands.                         |
| `InspectorModel`                                                                                                           | Per-representation `UpdateTransformControl` expanders, structure-type and model-index editors | Shows the node the user picked in the scene: a `ParamsFormModel` for its transform plus contributed sections.                                                                                                                                                            |
| `ParamsFormModel`                                                                                                          | `ParameterControls` state, `TransformControlBase`                                             | Generic PD form: `values`, per-field `FieldModel` (value, visible via `hideIf`, essential/advanced, error), `isDirty`, `apply()`, `reset()`, auto-apply policy (`canAutoUpdate`, debounced), optimistic projection (§5). Used by inspector, commands, settings, dialogs. |
| `TransformEditorModel(cell)` / `ActionFormModel(cell, action)`                                                             | `UpdateTransformControl`, `ApplyActionControl`                                                | `ParamsFormModel` + the apply path (`state.updateTransform`, `State.ApplyAction`).                                                                                                                                                                                       |
| `StateTreeModel` → `CellNode`                                                                                              | `state/tree.tsx`                                                                              | Keyed by ref; `isCurrent`, `isCollapsed`, flags; `actions()` computed on demand from `actions.fromCell` (also feeds contextual commands).                                                                                                                                |
| `SelectionModel`                                                                                                           | Selection toolbar, selection stats                                                            | `stats`, `entries`, `history`, `granularity`, `mode`; actions are commands (select by query, by residue list, modify, clear, create component, apply theme).                                                                                                             |
| `FocusModel`                                                                                                               | Focus row                                                                                     | `current`, `history`, `targets` (ligands etc. per structure, computed lazily).                                                                                                                                                                                           |
| `MeasurementsModel` → `MeasurementNode`                                                                                    | Measurements card                                                                             | Typed nodes keyed by ref; adding is a command family reading `selection.history`.                                                                                                                                                                                        |
| `SnapshotsModel` → `SnapshotNode`                                                                                          | `state/snapshots.tsx`, viewport stepper                                                       | `entries`, `current`, `isPlaying`; local and remote lists separate.                                                                                                                                                                                                      |
| `AnimationModel`, `TrajectoryModel`                                                                                        | Animation and trajectory viewport controls                                                    | Plus commands.                                                                                                                                                                                                                                                           |
| `ViewportModel`                                                                                                            | Viewport controls, settings, illumination, XR, screenshot                                     | `canvas` source + `ScreenshotModel`; `canvasTheme` (light/dark chrome over the canvas, from the background: `model/canvas-theme.ts`).                                                                                                                                    |
| `SequenceModel`                                                                                                            | `sequence.tsx`                                                                                | Structure/mode/entity/chain choices as states; residue markers stay imperative inside the sequence element (thousands of residues).                                                                                                                                      |
| `TaskListModel`, `NotificationsModel`, `LogModel`                                                                          | Task overlays, toasts, log                                                                    | From `tasks`, `toasts`, `log` sources plus command outcomes.                                                                                                                                                                                                             |
| `HelpModel`                                                                                                                | Help and bindings tables                                                                      | Bindings from behavior params + trackball + command keybindings.                                                                                                                                                                                                         |

Extensions contribute their own view models (e.g. `VolsegModel` over volseg entry data, `PaePlotModel`); see
[shell.md](shell.md) §5.

### 4.3 Lifetimes

- `PluginUIModel` is created by `PluginUIContext.init()` and retired by its `dispose()`.
- Area models (`scene`, `selection`, …) are children mounted by the root's `onMount`.
- View-only models (`stateTree()`, `snapshots()`, an inspector form) are created when their view becomes visible and
  unmounted when it closes. Expensive projections therefore cost nothing while hidden.
- Node models in keyed lists are owned by their list and retired when their key disappears.

## 5. Optimistic edits

The state tree updates asynchronously (`updateTree` is a Task; reprs recompute). The transaction primitives fit two
cases where the UI should not wait:

1. **Flags** (visibility, collapse): the toggle contributes
   `{ apply: n => ({ ...n, isHidden: !h }), isAbsorbed: n => n.isHidden === !h }` to the node's projection; the source's
   next value absorbs it. A failed update rejects it and the node falls back.
2. **Parameter edits** (isovalue slider, opacity, color): `ParamsFormModel` contributes the new values; the contribution
   is absorbed when the cell's `transform.version` changes and its params equal the submitted ones, and rejected when
   the update task fails or is aborted (with a notification). Rapid edits stack and only the last one matters, because
   each contribution is a full desired state.

Everything else (adding components, loading, superposition) is not optimistic: it shows task progress instead.

## 6. Actions

Operations live in `model/actions/` as plain functions `(plugin, args, signal?) => Promise<Result>`: one implementation
per operation, called by view-model methods and by commands, so a button and its command cannot drift (stories-untold's
"commands run the same actions as the tools"). Most are thin wrappers over manager methods or `PluginCommands`; the few
operations whose logic lives in views today move here (superposition modes, lazy volume loading, quick styles, focus
targets).

## 7. React and Preact

`ui-next/react` exports `useObservable`, `useObservableProp`, `useReactiveModel`, and `useUIModel()` (context), so a
React or Preact view reads the same models:

```tsx
function StructureLabel({ node }: { node: StructureNode }) {
  const label = useObservable(node.label);
  return <span>{label}</span>;
}
```

`useBehavior` from `@molstar/plugin-ui` keeps working on raw plugin behaviors; new code should prefer view models.

## 8. Testing

View models are tested without a renderer (jest, as the rest of the repo):

- Unit: construct with fake sources, `mount()`, drive states, assert fields and calls.
- Integration: a real `PluginContext` without a canvas, loading local fixtures (`data/`), asserting the scene tree after
  builder calls. The same setup powers Storybook's integration stories ([plan.md](plan.md) §3).
