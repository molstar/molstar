# View models

ui-next's state layer: renderer-free models built on the primitives in `reactive/` (`observableState`, `ReadonlyState`,
`ReactiveModel`, transactions; any-ui's state code, [architecture.md](architecture.md) §4) that turn the plugin's
managers, events, and state tree into UI-shaped, granular state. Every view (the built-in DOM runtime, React, Preact, a
test) reads the same view models.

**State (mocks).** The view-model contracts exist as interfaces in `model/*.ts` (§4.2), with `ParamsFormModel`
(`forms/params-form`), `StatusModel` (`shell/status`) and `model/canvas-theme` as the only implementations in the
package. The Storybook mocks in `_stories/mocks/` implement the interfaces in memory (`mockScene`, `mockSelection`,
`MockAgentModel`, …) and views take them as props. The mocks skipped the plan's fake-sources level: they implement view
models directly, so sources (§3) are still a design. Not yet: sources, `PluginUIModel` (§4.1), `keyedList`, `fromEvent`
and the mount signal (§2), optimistic edits (§5), `model/actions/` (§6). All of these are Phase 4
([../plan.md](../plan.md)).

## 1. Short answer to "do view models wrap manager behaviors?"

Yes, in two layers:

1. **Sources** turn plugin APIs into clean `ReadonlyState`s and event streams, one per plugin area, and hide the gaps
   listed in [../analysis/plugin-state.md](../analysis/plugin-state.md) §2 (void events, in-place mutation, late
   binding, silent writes, non-reactive registries). Sources know Mol*; they know nothing about panels.
2. **View models** are shaped by what a view shows: keyed trees of small objects (one per structure, component,
   representation, cell), each field its own `ReadonlyState`, plus methods that perform actions. View models depend on
   source **interfaces**, so Storybook and tests replace sources, not the plugin. Until sources exist, views depend only
   on the view-model interfaces in `model/`, so a plugin-backed implementation replaces a mock without view changes.

```
PluginContext ──► sources (Mol*-aware, UI-agnostic) ──► view models (UI-shaped) ──► views (dom / React / Preact)
      ▲                                                        │
      └──────────────────────── actions ◄──────────────────────┘ (also called by commands)
```

## 2. Primitives

In `reactive/`, any-ui's code copied in as Mol\* code ([architecture.md](architecture.md) §4):

- `observableState<T>(initial, { compare })` → `ObservableState<T>` (`value`, `asReadonly`, `set`; a record state also
  has `update(partial)`); `compare` defaults to `sameShallow`. A `set` made while subscribers are being notified is
  queued and delivered after, so every subscriber sees the values in order and ends on the latest. `ReadonlyState<T>` is
  `Observable<T> & { readonly value: T }`.
- `ReactiveModel`: `mount`/`unmount`/`retire` (not overridable) with the hooks `onMount`/`onUnmount`/`onRetire`;
  mount-scoped `subscribe`/`registerDispose`, permanent `subscribePermanently`/`registerPermanentDispose`; one owner (a
  second mount while mounted warns once).
- `combineBehaviors({ sources, project })`, `createDisposableGroup`, `sameShallow`, `shallowEqual`.
- `createUiTransactionCoordinator`, `createOptimisticProjection` (§5).

Outside `reactive/`: `ParamsFormModel` (`forms/params-form`, §4.2), the form every model uses for PD params, and the DOM
runtime's keyed `<For each key>` (`dom/control-flow`), which keeps an item's nodes and scope while its key stays in the
list and renders it again when the key stays but the item is a different value.

Additions ui-next needs (not yet; Phase 4, when sources arrive):

- **Mount-scope `AbortSignal`** (`this.signal` inside `onMount`), aborted on unmount, instead of hand-rolled sequence
  guards; every async action started by a model takes it. (Commands already get a signal: `run(args, ctx, signal)`.)
- **`keyedList(source, key, create, update)`**: maps a `ReadonlyState<readonly T[]>` to a
  `ReadonlyState<readonly VM[]>`, reusing the view model for a key, calling `update(vm, item)` on survivors, retiring
  removed ones. The emitted array is new only when membership or order changes; field changes flow through each view
  model's own states, so `<For>` touches only what changed.
- **`fromEvent(trigger, read, compare)`**: a `ReadonlyState` re-read on a void event (`selection.events.changed` →
  `selection.stats`) with equality suppression. Most sources are built with it.

Conventions (settled by the mocks):

- Contracts are interfaces in `model/`; implementations (mocks now, plugin-backed later) are classes or plain objects
  behind them. Views import the interfaces; `components/` imports nothing from `model/` (its only state import is
  `ParamsFormModel` from `forms/params-form`, for the params form).
- Public fields are `ReadonlyState`; writable states are private (`#label = observableState(...)`, exposed as a getter
  or `readonly label = this.#label.asReadonly`). Changes go through methods. The one exception is
  `ImageExportModel.autoCrop`, an `ObservableState` the view's toggle button sets.
- Lists take one of two shapes:
  - **Live nodes** where rows are many and change in place (the scene tree): a `SceneNode` keeps its identity per
    state-tree ref, every field is its own state, `children` is a state of nodes.
  - **Records** where items are small and change as a whole (`StateTreeNode`, `SnapshotEntry`, `MeasurementPoint`,
    `RemoteState`, `CommandRun`, `AgentTurn`): plain objects with an `id` (`ref` for state-tree cells) in one
    `ReadonlyState<readonly T[]>`; a change replaces the record, never mutates it. Views key them by id (`<For key>`); a
    row that must not restart (a run's spinner, an open fold) follows its record by id and derives per-field streams
    with `distinctUntilChanged` (`RunRow` in `views/agent/agent-turn`).
- Lists that come from registries or are fixed for the model's life are plain arrays (`SelectionModel.actions`,
  `MeasurementsModel.kinds`, `PresetsModel.groups`, `QuickStylesModel`, `SettingsModel.groups`, `HomeModel.formats`).
- Editable params are a `ParamsFormModel` owned by the model, the same instance each time it is asked for
  (`inspect(node).form`, `entryOptions(id)`, `styleOf(node)`), so a form keeps its state while views come and go.
- Optional parts are optional members, present only when the spec or the source supports them (`HomeModel.remoteStates`,
  `RemoteStatesModel.upload`/`download`/`remove`, `MeasurementsModel.highlightPoint`, `SceneNode.menu`).
- Expensive children are lazy: `SceneNode.lazy` with `load()` computes a structure's contents (entities, then chains)
  when first expanded, so a mesoscale structure with thousands of chains costs nothing until someone looks.
- State-like things are behaviors, event-like things are `Subject`s that never replay (`CommandManager.outcomes`). A
  request a view acts on when it shows is state (`ShellModel.tabFocus`, an open overlay's `focus`), and a message to
  show is `StatusModel.post` into its `current` state; no model interface exposes a `Subject` today.
- Errors are values (`StateTreeNode.status`/`error`, a `CommandRun` with `status: 'failed'` and its `detail`), never
  stream errors.
- Streams are created once per model (fields, cached getters, `ParamsFormModel.field(name)` cached per field), never per
  render (`ParamsFormModel.visible(name)` is not cached; the params form asks once per row).
- Dependencies come through the constructor (the mock factories take a `post` for status messages, or callbacks like
  `mockScene`'s `onAction`, and what they build on: the `ShellModel` for `mockSelection`, the `MockScene` for
  `mockMeasurements` and `mockSuperposition`); upward calls through callbacks.
- **Fewer re-renders: combine, then debounce.** Where a view depends on several streams (a row's value, its visibility,
  its help; a tree row's label, description, visibility), combine them (`combineLatest`) and add a short debounce (a
  frame: `debounceTime(0)`, or `auditTime(16)` with a scheduler) where a burst of updates is expected (state-tree
  rebuilds, selection changes), so one change renders once. Not yet: today the params form's rows combine visibility
  with the filter, and a tree row binds each field to its own node (only expanded-and-has-children is combined); the
  combining and the debounce come with the sources, since the mocks change one field at a time.

## 3. Sources

**Not yet** (Phase 4: `model/sources/*` with tests against a real plugin and fixtures). The mocks implement the view
models directly; the table below is the design the plugin-backed view models are built on.

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

**Records per node.** A tree's source hands out one record per node (by ref) that changes on its own, not a tree rebuilt
on each change: views render rows that follow their record by id, and `For` renders an entry again when its record is
replaced, so a rebuilt tree costs a subtree per change ([../spike.md](../spike.md): 24 ms per relabel at 5,000 cells,
against 20 relabels in one frame with records per node). Long lists and trees are virtualized
(`components/virtual-list`: the shown nodes flattened, only the window rendered).

**Identity.** Sources key everything by state-tree ref. A hierarchy rebuild that recreates `StructureRef` objects keeps
the same `StructureSource` instance for the same ref; its fields update only when `transform.version`, `cell.state`, or
`cell.status` change.

**One change, one emission.** Code that sets several states (or one several times) for one change wraps it in `batch`
(`reactive/observable-state`): each state emits once, with its final value, when the batch ends, so views render the
change once ([../spike.md](../spike.md): a burst of ten relabels per row, 12.9 → 2.4 ms). Lists pick the current item
with `keyedSelector`, so a selection change notifies the two rows it concerns, not every row.

**Batching.** Sources update after a state transaction ends (`changed` with `inTransaction: false`) and after animation
stops, as the managers already do, so views do not render intermediate trees. What still arrives in bursts (state-tree
rebuilds, selection changes) is combined and debounced by a frame (§2, conventions).

## 4. View models

### 4.1 Root

**Not yet.** In the mocks each view takes its models as props (`StructureTab` takes `scene`, `presets`, `quickStyles`,
`options`; `ShellView` takes `status`), and stories build them with the mock factories.

`PluginUIModel` (one per plugin UI; `PluginUIContext.model`). The context holds the infrastructure (`shell`, `commands`,
`keybindings`, `contributions`, `sources`; [shell.md](shell.md) §6); the root model holds the area view models built on
it. A sketch, in the names of the interfaces of §4.2:

```ts
class PluginUIModel extends ReactiveModel {
  constructor(ctx: PluginUIContext);
  readonly scene: SceneModel;            // structures, volumes, measurements, extension rows as one tree
  readonly selection: SelectionModel;
  readonly status: StatusModel;          // messages, tasks, chips; fed from the tasks, log and toasts sources
  readonly hoverLabel: HoverLabelModel;
  readonly canvasTheme: Observable<CanvasTheme>;
  // created lazily by the views that need them:
  stateTree(): StateTreeModel; snapshots(): SnapshotsModel; measurements(): MeasurementsModel; ...
}
```

### 4.2 Catalogue

The interfaces in `model/`, their mock, and what the plugin-backed implementation is built over (Phase 4):

| Interface (`model/…`)                                          | Shape                                                                                                                                                                                                                                                                                                                                                                                                                                            | Mock                               | Built over                                                                                                                                                                          |
| -------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `SceneModel` → `SceneNode` (`scene`)                           | `roots`, `picked`, `pick(node)`, `inspect(node)` → `SceneInspection` (`form`, `slot` `inspector:<kind>` for contributed sections, read-only `facts`, `actions`). A node: `id` (the ref), `kind`, `icon`, `label`, `description`, `isHidden`, `isExpanded`, `children`, `lazy`/`load()`, `menu` (`SceneMenuItem`s, nested `items` open inline, `danger` last), `removable`; `toggleVisibility`, `toggleExpanded`, `highlight`, `focus`, `remove`. | `mockScene`                        | Structure and volume hierarchy managers, keyed by ref; the measurement manager for the Measurements row; extension rows from contributions; the cell's transform for the form (§5). |
| `PresetsModel`, `QuickStylesModel` (`scene`)                   | Presets: `groups` (as the registry groups them: ungrouped first), `current`, `apply(id)`. Quick styles: `representations`, `styles`, `applyRepresentation`, `applyStyle`.                                                                                                                                                                                                                                                                        | `mockPresets`, `mockQuickStyles`   | The structure representation preset registry; the commands `structure.preset` and `canvas.style`, so the command menu offers the same.                                              |
| `SelectionModel` (`selection`)                                 | `summary` ("12 residues", undefined when empty), `actions` (`SelectionAction`: `id`, `label`, `icon`, `run`), `clear`. What to do with the selection lives in the Selection tab's sections (slot `selection`: theme, component, Measure, Superpose), contributed through the shell.                                                                                                                                                              | `mockSelection`                    | The structure selection manager; the actions are commands of the `selection` group.                                                                                                 |
| `MeasurementsModel` (`measurements`)                           | `points` (`MeasurementPoint` records, newest first; the first four are used), `movePoint`, `movePointTo`, `removePoint`, `focusPoint`, `highlightPoint?`; `kinds` (distance, angle, dihedral, label, orientation, plane, with their point count); `add(kind)`; `options` (unit, label size and color); `added` (the Measurements row's nodes) with `styleOf(node)`.                                                                              | `mockMeasurements`                 | The measurement manager and the selection additions history.                                                                                                                        |
| `SuperpositionModel` (`superposition`)                         | `available` (two or more structures), `form` (method: chains, TM-align, UniProt, ligands, atoms; and its options), `entries` (what each structure contributes), `ready`, `hint` (what to select), `focus(id)`, `superpose()`.                                                                                                                                                                                                                    | `mockSuperposition`                | The selection and the superposition utils (logic out of the old controls, §6).                                                                                                      |
| `SnapshotsModel` (`snapshots`)                                 | `entries` (`SnapshotEntry` records: `name`, `key`, markdown `description`, `timestamp`, `image`), `current`, `saveOptions`, `entryOptions(id)`; `add`, `apply`, `update`, `replace`, `move`, `moveTo`, `remove`, `clear`; `download('session' \| 'state')`, `open(files)`.                                                                                                                                                                       | `mockSnapshotsModel`               | The snapshot manager, `PluginCommands.State`.                                                                                                                                       |
| `SnapshotStepperModel`, `TrajectoryControlsModel` (`playback`) | Stepper (the playback card): `count`, `index`, `title`, `titles`, markdown `description`, `isPlaying`, `previous`, `next`, `goTo`, `togglePlay`, `openLink(href)` (`#key` applies a snapshot). Trajectory: `label` (its structure; with several, the focused one), `count`, `index`, `isPlaying`, `seek`, `togglePlay`.                                                                                                                          | `mockSnapshots`, `mockTrajectory`  | The snapshot manager; `UpdateTrajectory`.                                                                                                                                           |
| `HoverLabelModel` (`playback`)                                 | `lines` (empty when nothing is hovered).                                                                                                                                                                                                                                                                                                                                                                                                         | Inline in stories                  | `behaviors.labels.highlight`.                                                                                                                                                       |
| `ImageExportModel` (`image-export`)                            | `preview`, `aspect`, `autoCrop`, `crop`/`setCrop` (`CropBox` in fractions, like `ViewportScreenshotHelper.relativeCrop`), `size`, `form` (resolution, format, transparency, axes, illumination), `busy`, `copy`, `download`.                                                                                                                                                                                                                     | `mockImageExport`                  | `ViewportScreenshotHelper`.                                                                                                                                                         |
| `StateTreeModel` (`state-tree`)                                | `roots` (`StateTreeNode` records: `ref`, `label`, `description`, `icon`, `transform`, folded `decorators`, `children`, `visible`, `status`, `error`), `selected`/`select`, `inspect(ref)` → `StateTreeInspection` (`updates`: the cell's transform, then its decorators; `actions`: applicable state actions; both `HomeAction`s), `toggleVisible`, `remove`.                                                                                    | `mockStateTree`                    | `plugin.state.data`; updates and actions from `actions.fromCell` (also contextual commands).                                                                                        |
| `HomeModel`, `RemoteStatesModel` (`home`)                      | `actions` (`HomeAction`: `id`, `title`, `icon`, `description`, `form`, `applyLabel`, `busy`, `apply`), `formats` (the drop hint), `openFiles(files)`, `remoteStates?` (`serverUrl`, `setServerUrl`, `states`, `loading`, `refresh`, `open`, `upload?`, `download?`, `remove?`).                                                                                                                                                                  | `mockHome`                         | The root's applicable state actions (`applyLabel` from `display.triggerLabel`), file-format detection; Remote States only when the spec turns them on.                              |
| `SettingsModel` (`settings`)                                   | `groups` (`SettingsGroup`: `id`, `title`, `icon`, `form`, `wideLabels`) by intent: Appearance, Camera, Interaction, Performance, Interface, Behaviors, Advanced. Changes apply as they are made; no resets.                                                                                                                                                                                                                                      | `mockSettings`                     | Canvas3D, Canvas3DContext, the layout and `plugin.state.behaviors` params mapped into the groups; unmapped ones land in Advanced.                                                   |
| `canvas-theme` (functions)                                     | `luminance`, `themeForLuminance`, `themeForColor`, `backgroundLuminance` (a gradient's colors averaged by its opacity; images and skyboxes count as the background color; `undefined` when transparent), `canvasThemeOf`, `canvasTheme(props$)` (distinct), `parseColor`. Dark chrome below luminance 0.179.                                                                                                                                     | Stories pass `themeForColor(…)`    | The canvas props (`canvas3d.props`, re-read on `settingsUpdated`).                                                                                                                  |
| `AgentModel`, `AISettingsModel` (`agent`)                      | `turn`, `view` (`card`, `line`, `none`), `history` (`HistoryItem`: the person's runs, agent turns), `ask`, `stop`, `approve`, `undo` (the latest turn's checkpoint only), `dismiss`, `hide`, `reopen`; `settings`: `providers`, `setup`, `setKey` (`KeyScope`), `removeKey`, `choose`, `refreshModels`, `listModels`, `addProvider`.                                                                                                             | `MockAgentModel`, `MockAISettings` | any-harness over the command registry ([agent.md](agent.md)). Opt-in: only the agent views and stories import them (`views/agent/_test/boundary.test.ts`).                          |
| `StatusModel` (`shell/status`, a `ReactiveModel` class)        | `current` (the message on show: fades, errors stay), `history`, `tasks` (`StatusTask`: `progress`, `overlay`, `abort`), `chips` (`activate`, `dismiss`); `post`, `dismiss`, `setTasks`, `setChips`.                                                                                                                                                                                                                                              | Used as is                         | The `tasks`, `log` and `toasts` sources and command outcomes; selection and focus chips.                                                                                            |

`ParamsFormModel` (`forms/params-form`) is the form inside all of them: `params`, `values`, `field(name)` (one cached
stream per field, distinct), `visible(name)` (its `hideIf` on the current values), `set`, `load(values)` (authoritative
values, not reported as a change), `reset()`, `layout()` (essentials, then advanced), and
`options.onChange(values, changed)`. The owner decides what a change does: the inspector and settings apply it as it is
made; Home actions and the state tree's updates and actions keep it until their button (`HomeAction.apply`). Not yet:
dirty state, per-field errors, an auto-apply policy (`canAutoUpdate`, debounced) and the optimistic projection (§5).

The design's other view models, and where the mocks put them:

| Planned                                                        | Replaces (today)                                                                              | State (mocks)                                                                                                                                                                                                                 |
| -------------------------------------------------------------- | --------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `InspectorModel`                                               | Per-representation `UpdateTransformControl` expanders, structure-type and model-index editors | `SceneModel.inspect(node)`: the node's form plus contributed sections; says when the item has nothing to set.                                                                                                                 |
| `TransformEditorModel(cell)` / `ActionFormModel(cell, action)` | `UpdateTransformControl`, `ApplyActionControl`                                                | `HomeAction` (a form and its apply) for Home, the state tree's updates and its actions; the apply path (`state.updateTransform`, `State.ApplyAction`) is Phase 4.                                                             |
| `TaskListModel`, `NotificationsModel`, `LogModel`              | Task overlays, toasts, log                                                                    | One `StatusModel`; the log is its `history`.                                                                                                                                                                                  |
| `FocusModel`                                                   | Focus row                                                                                     | Not yet: focusing is `SceneNode.focus()` and `MeasurementsModel.focusPoint`; `StatusModel` has chips for it, which no story sets yet. Planned: `current`, `history`, `targets` (ligands etc. per structure, computed lazily). |
| `AnimationModel`                                               | Animation viewport control                                                                    | Not yet: the Animation card is a shell overlay of sections (slot `animation`: Play over a `ParamsFormModel`, the wiggle; `_stories/mocks/animation`).                                                                         |
| `ViewportModel`                                                | Viewport controls, settings, illumination, XR, screenshot                                     | Split: canvas settings are a card of sections (`_stories/mocks/canvas-settings`), image export is `ImageExportModel`, the chrome's theme is `model/canvas-theme`.                                                             |
| `SequenceModel`                                                | `sequence.tsx`                                                                                | Not in `model/`: `views/sequence` takes `SequenceSelector`s (structure, mode, entity, chain, instance: `value`, `set`, `visible`) and a `SequenceChain`; residues stay imperative inside the strip (thousands of residues).   |
| Picking mode and level                                         | Selection mode toggle, granularity select                                                     | Not in `model/`: `views/canvas/mode-switch` props (`mode`, `setMode`, `granularity`, `granularities`, `setGranularity`), from `behaviors.interaction.selectionMode` and the picking granularity.                              |
| `HelpModel`                                                    | Help and bindings tables                                                                      | Not yet. Bindings from behavior params + trackball + command keybindings.                                                                                                                                                     |

Extensions contribute their own view models (e.g. `VolsegModel` over volseg entry data, `PaePlotModel`); see
[shell.md](shell.md) §5. In the scene tree an extension adds a row under each structure it applies to (kind `extension`:
volume streaming, assembly symmetry; `MockRowSpec` in `_stories/mocks/structure-rows`): off until its eye turns it on,
then what it shows appears under it, and picked, its options are the inspector's form. The tree knows no extension by
name. The contribution API for such rows is not in `model/` yet.

### 4.3 Lifetimes

Planned with `PluginUIModel`; the mocks live as long as their story, and cache what a view asks for again (an inspection
per node, a form per snapshot).

- `PluginUIModel` is created by `PluginUIContext.init()` and retired by its `dispose()`.
- Area models (`scene`, `selection`, …) are children mounted by the root's `onMount`.
- View-only models (`stateTree()`, `snapshots()`, an inspector form) are created when their view becomes visible and
  unmounted when it closes. Expensive projections therefore cost nothing while hidden.
- Node models in keyed lists are owned by their list and retired when their key disappears.

## 5. Optimistic edits

**Not yet.** The mocks apply every change at once (an inspector edit updates the row's label and description in the same
call), so nothing waits; the transaction primitives are in `reactive/transaction`, unused by the models so far.

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

**Writing transform params.** Forms over a cell's transform (the inspector, the state tree's updates) emit the complete
new params, each value of its PD type (Vec3 and Mat4 as arrays of their size, colors as numbers, assets kept as `Asset`
objects, object lists as new lists), never a partial object or a shared reference, and skip the update when nothing
changed (`PD.areEqual`), so `state.build().to(ref).update(params)` records exactly the new version. A form takes the
authoritative params back with `ParamsFormModel.load` after an update made elsewhere.

## 6. Actions

**Not yet**: there is no `model/actions/`. In the mocks, model methods and the stories' commands (`mockCommands`) each
report what they would do through the status line, apart from the commands that open a panel or section. Commands are
also an agent's tools ([agent.md](agent.md) §2), so one action per operation serves the button, the command menu and the
agent.

Operations live in `model/actions/` as plain functions `(plugin, args, signal?) => Promise<Result>`: one implementation
per operation, called by view-model methods and by commands, so a button and its command cannot drift (as stories-untold
has it, commands call the same action functions as the agent's tools). Most are thin wrappers over manager methods or
`PluginCommands`; the few operations whose logic lives in views today move here (superposition modes, lazy volume
loading, quick styles, focus targets).

## 7. React and Preact

The `react/` modules (`react/use-observable`, `react/react-view`, …; no barrel yet) export `useObservable`,
`useObservableProp`, `useObservableSelection`, `useReactiveModel`, `reactView` and `MspElementHost`, and later
`useUIModel()` (context; it waits for `PluginUIModel`), so a React or Preact view reads the same models:

```tsx
function StructureLabel({ node }: { node: SceneNode }) {
  const label = useObservable(node.label);
  return <span>{label}</span>;
}
```

`useBehavior` from `@molstar/plugin-ui` keeps working on raw plugin behaviors; new code should prefer view models.

## 8. Testing

View models are tested without a renderer (jest, as the rest of the repo). Now: `forms/_test/params-form.test.ts`,
`model/_test/canvas-theme.test.ts`, view tests over the mocks (`views/structure/_test/structure-tab.test.tsx` renders
`StructureTab` over `mockScene`), and `views/agent/_test/boundary.test.ts` (only the agent views and stories import the
agent model). With sources (Phase 4):

- Unit: construct with fake sources, `mount()`, drive states, assert fields and calls.
- Integration: a real `PluginContext` without a canvas, loading local fixtures (`data/`), asserting the scene tree after
  builder calls. The same setup powers Storybook's integration stories ([../plan.md](../plan.md), story levels).
