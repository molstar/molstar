# The plugin's state in the Prism viewer

Where the viewer reads what `PluginContext` holds, and which of the plugin's managers it uses. The code is the truth for
how.

## 1. Engine state, area by area

Numbered rows (code cites `§1, row 2`); files under `src/sources/` unless a path is given.

| #   | Plugin state                                               | Read by → shown through                                           |
| --- | ---------------------------------------------------------- | ----------------------------------------------------------------- |
| 1   | `state.data`                                               | `features/state-tree/source.ts` → `StateTreeModel`                |
| 2   | `managers.structure.hierarchy`                             | `scene.ts` (a record per node) → `SceneModel`                     |
| 3   | a structure's entities and chains                          | `scene.ts`, on first expand → contents nodes                      |
| 4   | `managers.volume`, volume streaming                        | `scene.ts` → volume nodes, the streaming row                      |
| 5   | `managers.particles`                                       | `scene.ts` → particle nodes                                       |
| 5a  | shape providers in the state tree (no manager)             | `shapes.ts` → `scene.ts`'s Shapes group, `shape.*`                |
| 6   | `managers.structure.component`                             | `structure-tab.ts` → `PresetsModel`; `look.ts` → `SceneLookModel` |
| 7   | `managers.structure.selection`                             | `selection.ts`, `interaction.ts` → `SelectionModel`               |
| 8   | selection mode, `managers.interactivity.props`             | `interaction.ts` → `ModeSwitch`                                   |
| 9   | `managers.structure.focus`                                 | `interaction.ts` → `focus.*`                                      |
| 10  | `managers.structure.measurement`                           | `measurements.ts`, `features/measurements/source.ts`              |
| 11  | superposition (v5: in the view)                            | `features/superposition/operations.ts` → `SuperpositionModel`     |
| 12  | style layers (overpaint, transparency, emissive, clipping) | `style-layers.ts` → `representation.layers`                       |
| 13  | `managers.snapshot`, `.markdownExtensions`                 | `features/snapshots/source.ts` → `SnapshotsModel`                 |
| 14  | remote states                                              | `features/remote-states/` → `RemoteStatesModel`                   |
| 15  | trajectory frames                                          | `trajectory.ts` → the scene's `SceneTrajectory`, the canvas chip  |
| 16  | `managers.animation`, procedural animation                 | `features/animation/source.ts` → `AnimationModel`                 |
| 17  | `canvas3d.props`                                           | `viewer-sources.ts`, `settings.ts` → Canvas card, Settings        |
| 18  | `canvas3d.xr`                                              | `xr.ts` → toolbar XR item                                         |
| 19  | `managers.camera`                                          | `framing.ts`, `features/export-image/source.ts`                   |
| 20  | `helpers.viewportScreenshot`                               | `features/export-image/operations.ts` → `ImageExportModel`        |
| 21  | `behaviors.labels.highlight`                               | `pointer.ts` → `HoverLabelModel`                                  |
| 22  | `state.behaviors`, `config`                                | `settings.ts` → `SettingsModel`; config read at start             |
| 23  | `managers.toast`, `events.log`                             | `toasts.ts` → status column, History's Log                        |
| 24  | `managers.task`                                            | `viewer-sources.ts` `tasks` → status pill (Cancel), busy overlay  |
| 25  | registry: actions, formats                                 | `home.ts` (follows registrations) → `HomeModel.actions`           |
| 26  | registry: presets, themes                                  | `structure-tab.ts` → `PresetsModel`, `QuickStylesModel`           |
| 27  | sequence                                                   | `sequence.ts` → `SequenceViewers` ([sequence.md](sequence.md))    |
| 28  | bindings                                                   | `ui/bindings.ts` → Help tab                                       |

## 3. View state and undo

- **Undo is the viewer's history**, not the plugin's: checkpoints per step ([sequence.md](sequence.md) §7.2–7.3).
- **The selection's additions history is view state**: the ordered picks Measure reads and reorders; never recorded.
- **Selection mode, granularity and focus are the plugin's own state** (`interaction.ts`), so headless has them; view
  state, no step.
- **Canvas props are set through the plugin's emitting path** (`settingsUpdated`); a silent path is a core bug. XR is
  read by `xr.ts` once the canvas exists.

## 4. Reading the plugin: one shape for every source

Every plugin pattern becomes `ReadonlyState<T>` (synchronous `value`, emits only on change) held by a `SourceScope`
(`sources/state.ts`):

| #   | Pattern                               | Helper                                       |
| --- | ------------------------------------- | -------------------------------------------- |
| 1   | a `BehaviorSubject`                   | `fromBehavior`, with equality                |
| 2   | a mutable field and a void `changed`  | `fromChanged(read, event, compare)`          |
| 3   | a getter, no event (`canvas3d.props`) | set through the viewer's action, which emits |
| 4   | cells mutated in place, coarse event  | per-node records, rebuilt per changed cell   |
| 5   | events to accumulate (tasks, log)     | `scanList`                                   |

- **Row 4 is incremental per changed part** (packages/prism/AGENTS.md, bulk operations): every Mol\* change rebuilds the
  hierarchy, so a source that rebuilt everything at each would be N² over N loads. The scene's sync (`scene.ts`
  `syncRows`) keeps a structure's rows while its signature is the same (its hierarchy refs' cells, transforms, objects
  and hidden flags, its name, its hidden contents, its emptied components), and its clean-up looks only at the ids let
  go and the new ones; a measurement is described once per cell object and transform (`measurements.ts`). Guarded in
  `test:scaling` (`scene-sync.scaling.ts`, `measurements.scaling.ts`).
- **Only `api/`, `sources/`, `sequence/`, `testing/` and a feature's `source.ts`/`operations.ts` import the plugin**
  (`src/_test/boundary.test.ts`). Feature view models read sources and act through the API.

## 5. Which of the plugin's managers the viewer uses

Selective (README decision 24): engine work stays the plugin's; UI-era policy and side effects are the viewer's.

- **As they are**: `state`, `builders`, registries, `dataFormats`.
- **Read through `sources/`, written through `api/`**: the structure, volume, particle, interactivity, camera,
  animation, task, toast and markdown managers, `lociLabels`, `asset`.
- **`managers.snapshot`**: snapshots, session files, remote states; checkpoints use `plugin.state` directly
  (`sources/checkpoint.ts`).
- **`PluginCommands`**: not for UI work (a few state commands inside operations). **`layout`**: only its root, for the
  canvas host.
- **Not used**: `dragAndDrop`, `leftPanelTabName`, the custom controls and param editors (extension slots, rows and
  commands: [viewer-spec.md](viewer-spec.md) §3).
