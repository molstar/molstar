# The plugin's state in the Prism viewer

What `PluginContext` holds that a person sees or changes, and where each part lives in the Prism viewer: the source that
reads it, the view model the views bind to (the mocks' contracts in `packages/prism/viewer/src/model`,
`packages/prism/ui/src/model` and `ui/src/shell`), the Storybook story that shows it, and what is missing. The mocks
cover nearly all of it at a high level; this file is the checklist that the sources and the API (README §3.3) are built
against. The plugin's observability and its gaps are in [../analysis/plugin-state.md](../analysis/plugin-state.md).

Status (2026-10-15): the view models are interfaces with in-memory mocks; no source is built over a real plugin except
the first slice's canvas, hover lines, selection summary and framing (`viewer/src/sources`).

## 1. Engine state, area by area

| Plugin state                                                                        | View model (mock contract)                                                  | Shown in (Storybook)                                    | Missing                                                                                       |
| ----------------------------------------------------------------------------------- | --------------------------------------------------------------------------- | ------------------------------------------------------- | --------------------------------------------------------------------------------------------- |
| `state.data` tree, cells, flags, `actions.fromCell`                                 | `StateTreeModel` (`roots`, `selected`; `inspect` → updates and actions)     | Shell › Viewer mode; Internals › Performance › StateTree | undo is the viewer's (§3)                                         |
| `managers.structure.hierarchy` (structures, models, components, representations)   | `SceneModel` (`SceneNode` tree, `pick`, `inspect`), `SceneInspection` slots | Views › Structure tab (all); Shell › ManyStructures      | v5's hierarchy selection ("which structures the tools act on"): Prism acts on the picked node |
| a structure's contents (entities, chains), from its model                           | `SceneNode` kinds `contents`, `entity`, `chain` (lazy children)             | Views › Structure tab; Shell › ManyStructures            | its source (read from the `Structure`, not a manager)                                          |
| `managers.volume.hierarchy`                                                         | `SceneNode` kinds `volume`, `volume-representation`                         | Views › Structure tab › EmdbMap, EmdbModel, StreamingDensity | —                                                                                             |
| `managers.particles.hierarchy`                                                      | none (`SceneNode` has no particle kinds)                                    | none                                                    | particle nodes                                                                                |
| `managers.structure.component` (options, add, presets)                              | `PresetsModel`, `QuickStylesModel`, `mockComponentOptions`, node menus      | Views › Structure tab › Representation, Structure       | —                                                                                             |
| `managers.structure.selection` (entries, stats, additions history)                  | `SelectionModel` (`summary`, `actions`, `clear`)                            | Shell › Loaded (selection bar, Selection tab)           | additions history ("undo last selection"); granularity as state                               |
| `behaviors.interaction.selectionMode`, `managers.interactivity.props` (granularity) | `ModeSwitch` props (`CanvasMode`, `Granularity`)                            | Shell › Loaded                                          | their source and view model (§3)                                                               |
| `managers.structure.focus` (current, history)                                       | the selection bar's Focus action                                            | Shell › Loaded                                          | its source and view model (§3)                              |
| `managers.structure.measurement`                                                    | `MeasurementsModel`                                                         | Shell › Loaded (Selection › Measure)                    | —                                                                                             |
| superposition (logic in v5's view)                                                  | `SuperpositionModel`                                                        | Shell › Loaded (Selection › Superpose)                  | the algorithms move out of the view into the API                                              |
| `managers.snapshot` (entries, current, options, play)                               | `SnapshotsModel`, `SnapshotStepperModel`                                    | Shell › SessionWithPlayback; Canvas › Playback          | —                                                                                             |
| remote states (`PluginConfig.State.*Server`)                                        | `RemoteStatesModel`                                                         | Shell › RemoteStates, ManyRemoteStates; Views › Home    | —                                                                                             |
| trajectory frames (model index)                                                     | `TrajectoryControlsModel`                                                   | Canvas › Playback › Trajectory                          | which trajectory, with several loaded                                                         |
| `managers.animation` (current, params, play/stop)                                   | none: the Animation card is PD sections (`registerMockAnimation`)           | Shell › Loaded (Animation card)                         | an `AnimationModel` contract                                                                  |
| `canvas3d.props` (background, lighting, camera, clipping), `settingsUpdated`        | none: Canvas settings is PD sections (`registerMockCanvasSettings`); the ink | Shell › CanvasBackground                                | a `CanvasSettingsModel` contract                 |
| `canvas3d.xr` (`isSupported`, `isPresenting`), illumination                         | toolbar `rendering` menu                                                    | Shell › Loaded                                          | XR state in a model                                                                           |
| `managers.camera`, `camera.stateChanged`                                            | toolbar `reset` menu; `paddedFocus` (framing source)                        | Shell › Loaded; RealPlugin                              | —                                                                                             |
| `helpers.viewportScreenshot`                                                        | `ImageExportModel`                                                          | Shell › ExportAtSection                                 | —                                                                                             |
| `behaviors.labels.highlight`                                                        | `HoverLabelModel` (`hoverLines` source, built)                              | Shell › RealPlugin                                      | —                                                                                             |
| `state.behaviors` (behavior params, custom properties)                              | `SettingsModel` (group `behaviors`)                                         | Shell (Settings tab)                                    | —                                                                                             |
| `config` (servers, providers, performance)                                          | `SettingsModel` (groups `performance`, `advanced`); the viewer's options    | Shell (Settings tab)                                    | config has no change event; options are read at start                                        |
| `managers.toast`, `events.log`                                                      | `StatusModel` (`current`, `history`); Log tab                               | Shell (status column, Log)                              | —                                                                                             |
| `managers.task` (progress, finished)                                                | `StatusModel.tasks`; the busy overlay                                       | Shell › Loaded (density, surface tasks)                 | Cancel in the status line                                     |
| `managers.dragAndDrop`, open files                                                  | `HomeModel.openFiles`; `ShellView`'s `onDropFiles`                          | Views › Home tab                                        | —                                                                                             |
| registry (actions, formats)                                                         | `HomeModel.actions` (`HomeAction` with `form`, `applyLabel`), `formats`     | Views › Home tab                                        | a behavior keeping the latest registered set (§3)                |
| registry (presets, themes, selection queries)                                       | `PresetsModel`, `QuickStylesModel`; generated commands `preset.*`, `styles.*` | Shell › CommandMenuOpen                                 | same                                                                                          |
| `managers.markdownExtensions`                                                       | `SnapshotStepperModel.description`, `openLink`                              | Shell › SessionWithPlayback                             | —                                                                                             |
| sequence (v5's sequence panel)                                                      | `SequenceViewers` (`viewer/src/sequence`, built)                            | Shell › RealPlugin, RealPluginPredicted                 | ([../../seq/prism-sequence.md](../../seq/prism-sequence.md))                                 |
| bindings (behavior `bindings` params, trackball)                                    | Help tab sections (`bindings` of the viewer)                                | Shell › Help                                            | —                                                                                             |
| `layout`, `leftPanelTabName`                                                        | not used: the shell is the layout                                           | —                                                       | —                                                                                             |
| `customStructureControls`, `customImportControls`, `genericRepresentationControls`  | not used: extensions place their UI in the viewer's slots (README §3.1)     | `testing/structure-rows.ts` (volume streaming, symmetry) | —                                                                                             |

## 2. The viewer's and the UI's state

| State                                                    | Where                                                          | Shown in                                     |
| -------------------------------------------------------- | -------------------------------------------------------------- | -------------------------------------------- |
| docks, tabs, sizes, pinned, compact (`viewState`)        | `ShellModel` (`layout`, `dock(id)`, `setSize`, `setPinned`, …) | Shell › Pinned, Compact                      |
| open tab, card, section, page (`navigation`)             | `ShellModel` (`show`, `openOverlay`, `openPage`, `sectionOpen`) | Shell › ExportAtSection; Shell › Assistant    |
| theme, fullscreen                                        | `ThemeModel`, `FullscreenModel`                                | Overview › Foundations › Themes; ExpandByContext |
| the canvas ink                                           | `canvasInk` (`ui/src/model/canvas-theme`)                      | Shell › CanvasBackground                     |
| the sequence's style, colors, structure followed         | `SequenceViewer`                                               | Shell › RealPlugin                           |
| the agent's turn, history, settings                      | `AgentModel`, `AISettingsModel`                                | Shell › Assistant (all)                      |

`ShellModel` holds both rows today; README §3.3 splits it into the `navigation` and `viewState` managers.

## 3. What the mocks do not cover

- **Undo**: not the plugin's `canUndo` (its history is off) nor a manager's: the viewer snapshots both engines before
  an agent's turn and restores them on undo ([sequence.md](sequence.md) §7). The selection's additions history stays
  the plugin's.
- **Animation and canvas settings** are PD sections registered by the mocks, with no model: they need `AnimationModel`
  and `CanvasSettingsModel` contracts, so the API, commands and tools have something to call. Canvas props are set
  through one action; a path that changes them without `settingsUpdated` is a bug to fix in core.
- **Selection mode, granularity and focus** are the plugin's state (`behaviors.interaction.selectionMode`,
  `managers.interactivity.props`, `managers.structure.focus`); what is missing is their source and view model (the
  mocks pass them as view props), which the agent's account and the sequence read.
- **Particles** have no scene nodes; they get nodes, and Storybook gets stories over every data type (structures,
  models and trajectories, assemblies, maps and streamed density, segmentations, particles, shapes and meshes).
- **Tasks** in the status line cannot be cancelled; only a task that holds the plugin (the busy overlay) has Cancel.
  The running tasks are read from `events.task` (progress, finished), as the classic UI does; no core change.
- **Registry changes**: a behavior keeps the latest registered set as a state (David), which the menus, the generated
  commands and lazy extensions follow.

## 4. Reading the plugin: one shape for every source

The plugin's observability is uneven (David, 2026-10-15: "state being a variable and having a `changed` event only in
many of the places"), but changing it would break its API, so it stays. The viewer's `sources/` turn every pattern into
one shape, `ReadonlyState<T>` (a synchronous `value`, emits only on change, equality-suppressed), so views, view models
and the agent's account never see the plugin's patterns:

| The plugin's pattern                                    | Examples                                                                      | The source                                                                                |
| ------------------------------------------------------- | ----------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| a `BehaviorSubject`                                     | `behaviors.state.*`, `interaction.selectionMode`, `labels.highlight`, `focus` | taken as is, with equality                                                                |
| a mutable `state` field and a void `changed` event      | snapshots, animation, toasts, component options, interactivity props          | `fromChanged(read, event, compare)`: re-read on the event, emit when the value differs   |
| a getter and no event                                   | `canvas3d.props` (outside `settingsUpdated`), config                          | set through the viewer's own action, which updates the state; a silent path is a core bug |
| cells mutated in place, a coarse event                  | `state.data` cells, the structure hierarchy                                   | per-node records, rebuilt only for the cells that changed (view-models.md §3)             |
| events to accumulate                                    | `events.task` (progress, finished), `events.log`                              | a `scan` into a list state                                                                |

One helper per row lives in the viewer's `sources/`, tested against a real plugin; nothing else reads the plugin. A later
core cleanup (states instead of fields plus events) would change these helpers only.
