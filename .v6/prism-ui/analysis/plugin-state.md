# Plugin state and command surface

What the plugin runtime exposes for a view-model layer to project, and where it falls short. `C` =
`packages/plugin/core/src`; `K` = `packages/core/src`.

## 1. What is observable today

| Area                | Stream(s)                                                                                                                                                                        | Kind                                                               | Read on change from                                                            |
| ------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------ | ------------------------------------------------------------------------------ |
| Busy                | `behaviors.state.{isUpdating, isAnimating, isBusy}`                                                                                                                              | Behavior (`isBusy` debounced 750 ms; `plugin.isBusy` is immediate) | —                                                                              |
| Interaction         | `behaviors.interaction.{hover, click, drag, key, keyReleased, selectionMode}`                                                                                                    | Behavior (loci payloads)                                           | —                                                                              |
| Hover labels        | `behaviors.labels.highlight`                                                                                                                                                     | Behavior                                                           | —                                                                              |
| Left tab            | `behaviors.layout.leftPanelTabName`                                                                                                                                              | Behavior                                                           | —                                                                              |
| Data/behavior state | `state.data.events.cell.{created, removed, stateUpdated}`, `object.{created, removed, updated}`, `changed`, `historyUpdated`; `state.data.behaviors.{currentObject, isUpdating}` | Events + Behavior                                                  | `cell.status`, `cell.state`, `transform.version`, `canUndo`, `latestUndoLabel` |
| Structure hierarchy | `managers.structure.hierarchy.behaviors.selection`                                                                                                                               | Behavior (whole hierarchy)                                         | `hierarchy.current`, `selection`                                               |
| Volume / particles  | `managers.{volume,particles}.hierarchy.behaviors.selection`                                                                                                                      | Behavior                                                           | same pattern                                                                   |
| Selection           | `managers.structure.selection.events.{changed, additionsHistoryUpdated, loci.*}`                                                                                                 | Void events                                                        | `stats`, `entries`, `additionsHistory`                                         |
| Focus               | `managers.structure.focus.behaviors.current`, `events.historyUpdated`                                                                                                            | Behavior + void                                                    | `history`                                                                      |
| Measurements        | `managers.structure.measurement.behaviors.state`                                                                                                                                 | Behavior                                                           | —                                                                              |
| Component options   | `managers.structure.component.events.optionsUpdated`                                                                                                                             | Void                                                               | `state.options`                                                                |
| Interactivity       | `managers.interactivity.events.propsUpdated`                                                                                                                                     | Void                                                               | `props` (fresh copy per read)                                                  |
| Snapshots           | `managers.snapshot.events.{changed, opened}`                                                                                                                                     | Void                                                               | `state` (immutable `List` entries)                                             |
| Animation           | `managers.animation.events.{updated, applied}`                                                                                                                                   | Void                                                               | `state`, `current`, `getParams()`                                              |
| Toasts              | `managers.toast.events.changed`                                                                                                                                                  | Void                                                               | `state.entries` (immutable `OrderedMap`)                                       |
| Tasks               | `events.task.{progress, finished}`                                                                                                                                               | Events                                                             | no list; accumulate                                                            |
| Log                 | `events.log`                                                                                                                                                                     | Events                                                             | `log.entries` (reassigned)                                                     |
| Layout              | `layout.events.updated`                                                                                                                                                          | Void                                                               | `layout.state`                                                                 |
| Canvas settings     | `events.canvas3d.settingsUpdated`                                                                                                                                                | Void, incomplete                                                   | `canvas3d.props`                                                               |
| Canvas              | `canvas3d.{didDraw, commited, commitQueueSize, reprCount, resized}`, `camera.stateChanged`, `xr.{isSupported, isPresenting}`                                                     | Behavior                                                           | —                                                                              |
| Screenshot          | `helpers.viewportScreenshot.behaviors.{values, cropParams, relativeCrop}`                                                                                                        | Behavior                                                           | —                                                                              |
| Actions             | `state.data.actions.events.{added, removed}`                                                                                                                                     | Void                                                               | `actions.fromCell(cell, plugin)` (pull)                                        |

Registries (representations, themes, presets, selection queries, data formats, loci labels, drag-and-drop, markdown
extensions, custom properties) and `config` have **no change notification**; `query.structure.registry` has only a
`version` counter. `PluginContext.register(entry)` is the single entry point that changes them.

## 2. Gaps that matter for granular view models

1. **Coarse hierarchy updates.** The hierarchy behaviors emit one composite object per change; all `*Ref` objects are
   recreated on every rebuild while `cell` identities stay stable. Updates are suppressed during transactions and
   animation and flushed afterwards; the `current`/`selection` getters call `sync(false)` and can run ahead of
   `behaviors.selection.value`.
2. **Cells mutate in place.** `status`, `state` flags, `obj`, `params`, `errorText` change on the same object;
   `cell.stateUpdated` does not say which. `transform.version` changes with params/tags/parent, not with cell flags.
3. **Void events and shared mutable collections** (selection entries Map, history arrays, focus history) mean a view
   model must re-read and diff.
4. **Non-reactive registries and config.** Parameter option lists (`getAddParams`, theme lists) are computed when the
   schema function runs and go stale silently.
5. **Silent write paths.** `layout.setProps` (no event, no expand handling), direct `canvas3d.setProps` and
   `canvas3dContext.setProps` (no `settingsUpdated`).
6. **Late binding.** `managers.interactivity`, `managers.lociLabels`, `builders.structure` exist only after `init()`;
   `canvas3d` only after mount.
7. **Measurement churn.** `StructureMeasurementManager.sync` builds new arrays on almost every `changed`, re-emitting
   `behaviors.state`.
8. **Tasks.** No list of running tasks; `progress.id` (run id) differs from the id `requestAbort` takes (task id).

## 3. Commands and input

- `PluginCommands` (`C/commands.ts`): 32 typed message-bus commands (State, State.Snapshots, Interactivity, Layout,
  Toast, Camera, Canvas3D) with random UUID ids, no titles, no enumeration API, `Promise<void>` results, subscribers run
  serially. `Interactivity.Structure.Highlight/Select` have no subscribers (dead).
- Operations the UI performs **without** a command: undo, selection modify/clear, focus set/clear, animation play/stop,
  snapshot next/previous, selection-mode toggle, left-panel tab, screenshot, component add/update, every
  hierarchy/component/measurement/volume manager method. A command menu needs its own registry over these operations.
- Bindings (`K/util/binding.ts`) live as hidden `bindings` params on behaviors (hover/click select, focus, camera keys
  `I`/`O`/`T`/`G`/Shift+Space, snapshot Ctrl+arrows) and in Canvas3D trackball attribs (W/A/S/D/R/F/Q/E, drag/scroll).
  `InputObserver` forwards keys only while the pointer is over the canvas and the target is the body or canvas. There is
  no global keyboard service, no when-clauses, no conflict detection. The help view's lookup of `clickSelectToggle` is
  wrong (the key is `clickToggle`).

## 4. State tree and actions

- `State` exposes an immutable `tree`, mutable `cells`, the events above, `updateTree` (Task), `transaction`,
  `applyAction`, `updateCellState`, `undo`/`canUndo`/`latestUndoLabel`.
- Cell flags: `isGhost` (hide from UI), `isLocked`, `isHidden`, `isCollapsed`.
- `StateTransformer.definition`: `display { name, description }`, `params(a, plugin)`, `isApplicable(a, plugin)`,
  `canAutoUpdate`, `isDecorator`.
- `StateAction.definition`: `display`, `params(a, plugin)`, `isApplicable(a, transform, plugin)`, `run`.
- `state.data.actions.fromCell(cell, plugin)` is the applicable-action list for a cell (pull only).
- `PluginSpec.Action.customControl`/`autoUpdate` are dropped by `registry-entry.ts`.

## 5. Core follow-ups worth proposing

None of these blocks Prism (the source layer works around them), but each removes a workaround. The additive ones are
done in v6 (Phase 4); the two that break compatibility wait until the classic UI is retired and are recorded in
[../../plans/breaking-v6-changes.md](../../plans/breaking-v6-changes.md).

| Change                                                                                            | Removes workaround                       | Compatible                                  |
| ------------------------------------------------------------------------------------------------- | ---------------------------------------- | ------------------------------------------- |
| `plugin.events.registryChanged` (emitted by `register()` and its undo)                            | Re-deriving menus on a timer or per open | Yes: a new event                            |
| `cell.stateUpdated` payload says what changed (`status` / `flags`)                                | Diffing cell snapshots                   | Yes: optional fields                        |
| `layout.setProps` emits; `canvas3d.setProps` path emits `settingsUpdated`                         | Polling settings                         | Yes: existing events, emitted in more cases |
| `TaskManager` exposes running tasks as a behavior, one id space                                   | Accumulating progress events             | Yes: a new behavior                         |
| Stable string ids for `PluginCommands` (`'state.remove-object'`)                                  | Mapping UUIDs for telemetry/scripting    | Yes: ids next to the existing ones          |
| Fix help binding key                                                                              | —                                        | Yes                                         |
| Remove dead `Interactivity.Structure.*` commands and `customControl`                              | —                                        | No: deferred to the classic UI's retirement |
| Move `customStructureControls`/`customImportControls`/`genericRepresentationControls` out of core | React types in a non-UI package          | No: deferred to the classic UI's retirement |
