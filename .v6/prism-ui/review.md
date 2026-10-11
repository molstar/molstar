# Prism spec review (2026-10-11)

## Verdict

The contracts are not ready to freeze. One blocker remains, and a few majors sit underneath it:

- **The blocker:** `PrismContext` is meant to be the headless core, but it names viewer-mode models that live in an
  opt-in package. Nothing yet says how modes and modules attach their own models or how long those models live.
- **Structural majors:** the headless context has no command registry or contributions it can use, and the plugin import
  boundary cannot hold as written.

Most other majors are additive fixes to documents or contracts (sources, run attribution, a reactive `when`, the layout
spec). Two are local runtime bugs (the For region bug and the delegated event walk). One documentation pass covering B1,
M1–M5 and the Phase 4 plan rewrite is enough to freeze the contracts Phase 4 builds against. The command, shell and
agent contracts can then be frozen per phase.

## Must fix before the build

### B1. PrismContext names viewer models it cannot import, and model ownership and lifetime are undefined (blocker)

_Merges ctx-1, plan-1, the placement part of vm-1, the eager/lazy part of vm-7, and the lifetime part of plan-11._

**Problem.** `PrismContext` is planned for `@molstar/prism` and has these fields:

- `scene: SceneModel`
- `selection`
- `snapshots`
- "… measurements, superposition, playback, settings"

All of these contracts live in `@molstar/prism-viewer`, an opt-in mode. modes.md §2 says a story-only app ships none of
the viewer's code.

The documents describe three different shapes:

- **shell.md §6:** all fields eager.
- **view-models.md §4.1:** some eager, plus lazy `stateTree()`, `snapshots()` and `measurements()`.
- **architecture.md §2 and plan Phase 4:** "the model/ contracts move into core". But `ui/src/model` holds only
  canvas-theme, playback and settings.

`CommandContext` (commands.md §2.1) names `StructureNode`, `CellNode` and `FocusEntry`, which exist in no file. There is
no slot through which a mode, the stories module or the agent attaches typed models. There is also no lifetime rule: the
mock in `viewer/src/testing/viewer-mode.tsx:339-400` creates the models again on every activation and never retires
them.

**Where:**

- shell.md §6 lines 566-577 and 620-622
- view-models.md §4.1 lines 160-176 and §4.3 lines 233-237
- architecture.md §2 lines 32-34
- modes.md §1 line 17, §2 lines 29-34, §3 line 46
- plan.md Phase 4 lines 121-124
- commands.md §2.1 lines 130-140
- `packages/prism/viewer/src/model/*.ts`

**Change:**

1. Limit `PrismContext` to frame-wide infrastructure: plugin, sources, commands, status, history, settings, and the
   hover and canvas sources.
2. Add a typed module slot, for example `PrismModule { id; models(ctx): M; dispose() }` with `ctx.models.get(key)`.
3. Modes, stories and the agent attach their models through that slot. The viewer contracts and their plugin-backed
   implementations stay in `@molstar/prism-viewer`.
4. Give the viewer a headless entry (`@molstar/prism-viewer/headless`) so `createPrism` tests can include it.
5. Make `CommandContext` composable: a core part (busy, activeView, activeMode) plus parts each mode supplies (see M1).
6. State one lifetime rule (see Q1). Then make shell.md §6, view-models.md §4.1/§4.3, architecture.md §2 and plan Phase
   4 say the same thing.
7. List in Phase 4 exactly which model files move to core.

**Phase:** 4 (blocks the core carve-out and the `createPrism` test). This also shapes 5 and 6.

### M1. The command registry has no single owner, context type or package

_Merges ctx-3 and cmd-1._

**Problem.** The documents disagree with each other and with the code:

- shell.md §6 types `commands: CommandRegistry<PrismContext>` and says "commands get the model they work with, not the
  context".
- commands.md §2.1 gives `CommandContext` a `prism: PrismContext` field plus viewer-only fields (structures,
  sceneTarget, currentCell, selectionMode), and says "UI commands get PrismUIContext".
- commands.md §1 says the command core may move to `@molstar/plugin`. architecture.md §3 says it moves to
  `@molstar/prism`.
- commands.md §3 says built-ins are "Registered by `PluginUIContext`".
- No `CommandRegistry` type exists. The code has `CommandManager<Ctx>` over `ReadonlyState<Ctx>`.
- The mocks build one manager per mode (`viewer/src/testing/commands.ts:163`, `stories/src/testing/index.tsx:126`). The
  frame registers no commands, and the `mode.*` commands are registered into each mode's manager
  (`ui/src/testing/shell.tsx:316-329`).

The per-mode managers are a mock artefact (modes.md §1 already describes one registry). The document contradictions are
real.

**Change:**

1. Freeze one `CommandManager<CommandContext>` owned by `PrismContext`, living in `@molstar/prism`.
2. `CommandContext` is a narrow derived state: a frame base (busy, activeView, activeMode) plus mode slices looked up by
   key, for example `ctx.get(SceneCommandContext)`. It has no `prism` field.
3. Commands close over the models their contribution gives them.
4. The frame registers its commands once (mode.*, theme, command-menu.open). Modes register into the shared manager on
   activation.
5. Fix commands.md §1 (the core is `@molstar/prism`) and §3 (`PrismContext`, not `PluginUIContext`).
6. Replace `CommandRegistry<PrismContext>` in shell.md with the real `CommandManager<Ctx>` shape.

**Phase:** 4 (types move to core) and 5.

### M2. Headless and UI are not split for contributions and commands

_Merges ctx-2, plan-4 and the navigation-port part of cmd-1._

**Problem.** Commands, command groups, formats, `models?(ctx)` and `structureRows` exist only as fields of
`UIContribution`. Only `PrismUIContext` activates them, from `ui.contributions`. `PrismSpec` is only `{ plugin }`, so
`createPrism({ spec })` has an empty registry. That contradicts Decision 2 (the registry is headless and serves the
agent's tools) and plan Phase 5 ("Command tests run headless").

Phase 5 also builds pieces that need the UI and are only built in Phase 6:

- `opens` and `panel.*` (need `ShellModel`)
- "open as form" (needs a dialog host that no phase builds)
- the KeybindingService (a `PrismUIContext` field)

Nothing marks a command as UI-only, so the generated tool list would offer commands that fail headless.

**Where:** shell.md §5 lines 499-516 and 548; shell.md §6 lines 591-592, 620-621, 632-634; commands.md §2 lines 61-70
and §2.1 line 131; plan.md Phase 5 lines 162-167 and Phase 6 line 177.

**Change:**

1. Add `PrismContribution { id; commands; commandGroups; formats; models(ctx); structureRows }`, listed in
   `PrismSpec.contributions` and activated by `PrismContext.init()`.
2. Make `UIContribution extends PrismContribution` and add the view fields (tabs, sections, overlays, canvas parts,
   paramEditors, settings).
3. Ship `DefaultPrismContributions()` in core next to `DefaultUIContributions()`.
4. `opens` and panel commands reach the UI through a `navigate` port that `PrismUIContext` supplies. With no UI, or when
   the target is not registered, the command's `when` is false and its reason says why.
5. Mark UI-only commands so tool generation skips them.
6. Split Phase 5 into:
   - **5a (headless):** core commands over actions, generated and contextual commands, formats, `execute`, tests in
     Node.
   - **5b (after a `PrismUIContext` skeleton):** UI commands, the KeybindingService, the dialog host (M27).

**Phase:** 4-5.

### M3. The plugin import boundary cannot hold for extensions, mode sources and commands

_Merges ctx-9 and the boundary part of plan-2._

**Problem.** shell.md §6 and architecture.md §3 allow `@molstar/plugin` imports only in `@molstar/prism`'s
`model/sources` and `model/actions`, and the planned test checks exactly that. Several documented paths break the rule:

- commands.md §1 says a command's run "may dispatch a PluginCommand", and §2 puts `prism` (with its plugin) on
  `CommandContext`.
- Mode-specific sources (the state tree, MVS story rendering) would have to live in core. That contradicts modes.md §2.
- Every Phase 7 port wraps a controller that takes a `PluginContext` (GeometryControls, Mp4Controls, volseg, PAE).
- `models?(ctx: PrismContext)` gives extensions `ctx.plugin`.

"Not used, ever" is not quite true for PluginLayout: `layout.root` is the container the canvas is sized against. Prism
avoids the body mutation only because it leaves `spec.layout` unset.

**Where:** shell.md §6 lines 605-616 and §5 lines 515 and 551; architecture.md §3 lines 141-148; commands.md §1 line 11
and §2 line 131; plan.md Phase 7 line 201; `packages/plugin/core/src/context.ts:330-333, 390`.

**Change:**

1. Define the boundary by folder role, not by package. Any Prism package (core, a mode, an extension's prism entry) may
   have `sources/` and `actions/` folders that import the plugin.
2. Turn the boundary test into a reusable checker that runs on every package with a Prism module.
3. Commands call actions only. Remove "dispatch a PluginCommand" from commands.md §1.
4. Restate the PluginLayout rule: Prism never sets `plugin.layout`, and an action calls `plugin.handleResize()` when the
   shell's canvas inset changes.

**Phase:** 4 (the boundary test), 7 (extension ports).

### M4. The Phase 4 work list contradicts other decisions, and the core carve-out has a package cycle

_Merges plan-3, the IconName part of vm-1, and the cycle part of cmd-1._

**Problem:**

- The Phase 4 view-model order includes Sequence, which plan.md:192-193 says not to wire.
- It includes Help, which has no view (see Q9).
- It includes Tasks, Notifications and Log, which view-models.md:214 already folds into StatusModel.
- `commands/types.ts:10`, `model/settings.ts:9` and five viewer contracts import `IconName` from
  `components/icon/icon.ts`. That module is the DOM module; `IconName` is derived from the generated Lucide table.
  Moving `commands/` and `model/` into core is therefore not "a file move": it creates a core→prism-ui cycle. The
  architecture.md §3 line 112 rule ("icon names are data", `import type` allowed) needs revising for the carve-out, not
  enforcing.
- keyedList, fromEvent and the mount-scope AbortSignal are assigned to Phase 4 in view-models.md:12-14 but are missing
  from the plan text.

**Change:**

1. Rewrite the Phase 4 list in this order: the PluginSources interfaces and fake sources; the reactive additions
   (keyedList, fromEvent, AbortSignal); a core icon key (a string-literal type in core, which prism-ui maps to glyphs);
   the HistoryModel contract (M5); then the view models.
2. Drop Sequence and Help from the order.
3. Scene exposes a rows hook that the Phase 6 `structureRows` contributions feed.
4. Revise architecture.md line 112.

**Phase:** 4.

### M5. History and run attribution are defined only in the opt-in agent package

_Merges ctx-11 (history part), vm-9 and cmd-3._

**Problem.** shell.md §6 and view-models.md §4.1 give core `PrismContext` a `history: HistoryModel`. That interface
exists nowhere. `RunSource`, `CommandRun` and `HistoryItem` live in `packages/prism/agent/src/model.ts:12-71`, which
core must not import.

`CommandOutcome` (`ui/src/commands/types.ts:37-47`) carries no run id, source, arguments or timing, and `open()` takes
no source.

When an assistant answers inside a field, `invocation.ts:484-488` emits a `done` outcome. That enters recents and
History, and in the mock it releases the agent checkpoint (`mock-agent.ts:268-269`) although nothing ran.

Two things are unspecified:

- whether a run the person confirms after the assistant proposed it counts as the person's or the assistant's
- what an app without the agent shows

**Change:**

1. Define these in `@molstar/prism`:
   - `RunSource`
   - a core `RunRecord { id, commandId, title, args, source, turn?, status, startedAt, durationMs, result }`
   - `HistoryModel`, holding runs, with an item-kind extension point
2. The agent's `CommandRun` extends `RunRecord` with `read`, `approval` and `progress`, and the agent adds turns.
3. Every outcome carries its `RunRecord`.
4. Assistant answers get a separate `answered` outcome that is never recorded as a run.
5. A run the person confirms is `{ kind: 'person', via?: 'assistant' }`.
6. Without the agent, the app shows Log over StatusModel. With the agent, the agent contribution replaces Log with
   History (see M18).

**Phase:** 4 (contract), 5.

### M6. The source catalogue does not cover what view models are built over, and it still lists toasts

_vm-8, plus the toasts part of ctx-9._

**Problem.** Under the one-layer rule, view models see only sources and actions. Several entries in view-models.md
§4.2's "Built over" column rely on plugin APIs that no source in §3 provides:

- `state.data.actions.fromCell`
- PD definitions of transformers and actions (`definition.params(a, plugin)`)
- ViewportScreenshotHelper
- writable Canvas3D, Canvas3DContext and behavior params
- data formats and format detection
- trajectory frame count and index
- loci for entity and chain rows
- remote states

`PresetsModel.current` (`viewer/src/model/scene.ts:112`) has no source at all, because core does not record which preset
was applied.

§3, §4.1 line 170 and the §4.2 StatusModel row feed status from a `toasts` source. shell.md §6 says toasts are "Not
used, ever". MVS also logs its errors (`mvs/runtime/src/load.ts:210`), so the log covers them.

**Change:**

1. Add these sources to §3: `actions` (per ref, versioned), `params` (resolve a PD definition and current values for a
   ref), `screenshot`, `settings` (read and write, with change events), `formats`, `trajectory`, `remoteStates`, and
   loci helpers for derived scene rows.
2. Either define `PresetsModel.current` per focused structure, from a tag the preset action sets, or drop it.
3. Remove the `toasts` source and state that core toasts are covered by the log.

**Phase:** 4.

### M7. SceneNode identity does not work for rows without a cell, and lazy children go stale

_vm-2._

**Problem.** `viewer/src/model/scene.ts:30-34` keys every node by "the state-tree ref". The `contents`, `entity` and
`chain` rows have no cell, and an `extension` row exists before its cell does. The mock builds ids by string
concatenation (`testing/scene.ts:473-500`), which can collide with real refs.

`load()` computes children "once" (lines 49-50). A model-index or assembly change, or a preset re-creating cells, leaves
those children stale. Nothing defines what `picked` or expansion become when their node is retired.

**Change:**

1. Specify a scene-key scheme: `ref` for cells; `${structureRef}#contents`, `#entity:${id}`,
   `#chain:${entityId}:${labelAsymId}` and `#ext:${contributionId}` for derived rows.
2. Add `ref?: string` for the backing cell.
3. Recompute lazy children when the owning structure cell's `transform.version` or `obj` changes.
4. `picked` becomes undefined when its key disappears.
5. State whether expansion is UI state keyed by scene key or the cell's `isCollapsed` flag. (view-models.md §4.3 already
   retires nodes whose key disappears.)

**Phase:** 4.

### M8. StateTreeNode nests its children, which forces the tree rebuilds the spike ruled out

_vm-3._

**Problem.** `viewer/src/model/state-tree.ts:23` embeds `children?: readonly StateTreeNode[]` inside immutable records,
and `roots` holds the whole tree. A change to one cell replaces that record, every ancestor and `roots`. view-models.md
§3 "Records per node" rules this out, on the spike's evidence of 24 ms per relabel at 5,000 cells.

The status union also drops core's `processing` state (`core/src/state/object.ts:105`).

**Change:**

1. Use `roots: ReadonlyState<readonly string[]>` and `node(ref): ReadonlyState<StateTreeNode> | undefined`, with
   `children: readonly string[]` holding refs.
2. Map all four core statuses, or document that `processing` folds into `pending`.

**Phase:** 4.

### M9. Contracts have no results, busy states or node errors for async work

_vm-4._

**Problem.** These methods return `void`:

- `PresetsModel.apply` and the QuickStyles apply methods
- `SuperpositionModel.superpose`
- `SnapshotsModel.add`, `apply`, `replace`, `download`, `open`
- `MeasurementsModel.add`
- `ImageExport` copy and download
- `HomeModel.openFiles`
- `RemoteStates` refresh and open
- `SelectionAction.run`

Only `HomeAction.busy`, `ImageExportModel.busy` and `RemoteStates.loading` exist as busy states.

view-models.md §6 says actions return `Promise<Result>`, but the model methods drop the result. `SceneNode` has no
status or error, while core's `doError` leaves a failed cell in `error` with its objects removed and the update
resolving normally (`core/src/state/state.ts:1040-1060`).

**Change:** Freeze one rule (see Q4). Recommended:

1. Async model methods return `Promise<ActionResult>` (`{ ok: true } | { ok: false; message }`), never reject, and also
   post to StatusModel.
2. Each long operation exposes `busy: ReadonlyState<boolean>`.
3. Add `status` and `error` states to `SceneNode`, fed from the cell.

**Phase:** 4.

### M10. Registry-derived lists and form definitions are fixed at construction

_vm-6._

**Problem.** view-models.md §2 lines 79-80 make registry-derived lists plain arrays: `PresetsModel.groups`, QuickStyles,
`SelectionModel.actions`, `SettingsModel.groups` and `HomeModel.formats`. §3 lines 130-132 introduce `registry.version`
precisely because registries change, and contributions register into them after init.

`ParamsFormModel` takes a fixed `readonly params: P` and computes defaults once (`forms/params-form.ts:28-34`). Core
re-evaluates `definition.params(src, plugin)` per parent object (`state.ts:1175-1176`), so inspector option lists go
stale after a model or assembly change.

**Change:**

1. Make these lists `ReadonlyState<readonly …[]>`, recomputed on `registry.version` or on the relevant source.
2. Define what `SettingsModel` shows before the canvas is ready and in a headless context.
3. Give `ParamsFormModel` a `reset(params, values)` (or `params: ReadonlyState<P>`) that keeps the cached field streams
   and `load`s normalized values.
4. Require forms over cells to re-derive when the parent object changes.

**Phase:** 4.

### M11. Lazy models, inspections and shared nodes have no ownership in the contract

_vm-7, minus the eager/lazy part already in B1._

**Problem.** No interface in `viewer/src/model` extends ReactiveModel or exposes `dispose`, so "mounted by the context"
and "unmounted when the view closes" have no handle in the contract.

The documents also conflict on form lifetime. view-models.md §2 lines 81-82 says a form is "the same instance each time"
so it keeps its state. §4.3 lines 235-236 says it is unmounted when its view closes.

Two further problems:

- With one owner per ReactiveModel, a model shown in two places is torn down by whichever view closes first.
- `inspect(node)` caches inspections whose forms subscribe to cells, and nothing retires them. The mock's `Map` leaks
  (`testing/scene.ts:745-764`).

`MeasurementsModel` cannot be built over `SceneModel` either. The mock uses `scene.measurements` and
`scene.measurementStyle`, which exist only on MockScene (`testing/measurements.tsx:82-83`).

**Change:**

1. Lazily created view models are obtained as ref-counted leases: `{ model, release() }`, with an optional grace period
   so form state survives a close and reopen.
2. A SceneInspection belongs to its node and is retired with the node's key.
3. Add `measurements: ReadonlyState<readonly SceneNode[]>` to `SceneModel`, and define `styleOf(node)` as
   `scene.inspect(node).form`.

**Phase:** 4.

### M12. The optimistic-edit rules do not match State.updateTree

_vm-5._

**Problem.** view-models.md §5 lines 247-253 makes two assumptions that core does not meet.

- **Absorb on equal params.** A contribution is absorbed when "params equal the submitted ones". But `resolveParams`
  rewrites `transform.params` through `PD.normalizeParams` (`state.ts:1174-1181`), so equality may never hold and the
  edit stays pending forever.
- **Reject on task failure.** A transform error does not reject the update unless `rethrowErrors` is passed
  (`state.ts:259, 290`). And `revertOnError` re-runs the old tree, which matches neither rule (`state.ts:393-401`).

The flag projection `n => ({ ...n, isHidden })` is written for a record, but `SceneNode.isHidden` is a per-field state
on a live node.

**Change:**

1. A contribution settles on the first `stateUpdated` or `transform.version` change for its ref after submission.
2. `ok` absorbs the edit and `load`s the normalized params. An error, removal or abort rejects it and posts to status.
   No equality check.
3. Specify the projection per field.
4. Decide whether flags need optimism at all (see Q5).
5. Fold in rt-12: the code that called `begin()` owns resolving the transaction, and disposing a projection detaches it
   without finishing it.

**Phase:** 4.

### M13. For loses a row whose only node is a region placeholder

_rt-1._

**Problem.** `insertRegion` (`dom/render.ts:156-157, 173-183`) starts as one empty Text node. On the first value that is
not text, it wraps that node in comment anchors and removes it.

For (`dom/control-flow.ts:255-256`) treats a one-node row as `first === last`. A row whose only child is a Show, Dynamic
or stream that starts empty or as text therefore holds a detached node. Later moves, removals and `scope.dispose` act on
that detached node and leave the real content in place.

**Change:**

1. Keep the placeholder as the `end` anchor, or use the one-node shortcut only for Elements.
2. Add a `control-flow.test.tsx` case: a row that is a `<Show>` starting falsy, then reordered and removed.

**Phase:** all (scene tree, section and toolbar rows).

### M14. The delegated event walk follows the live parentNode

_rt-2._

**Problem.** `dispatchDelegated` (`dom/render.ts:279-301`) walks `node.parentNode` while handlers run, instead of the
path fixed at dispatch.

- A handler that detaches its own subtree (for example, closing the Show it sits in) cuts the walk short, so the
  ancestors' handlers never run.
- When an inner root is detached, the outer root continues from `parentNode === null` and stops.
- A ShadowRoot's `parentNode` is null, so the walk never crosses a shadow boundary.

**Change:**

1. Walk `e.composedPath()` up to the root's index, and record that index for outer roots.
2. Document in rendering.md §4.1 that a native listener calling `stopPropagation` above delegated content blocks the
   delegated handlers.
3. Add tests for a handler that removes its own subtree, for a removed inner root, and for a native `stopPropagation`.

**Phase:** all.

### M15. Focus, Escape and modality have no contract

_Merges ux-4, ux-5 and rt-7._

**Problem.** Escape order ("innermost first", shell.md §2 lines 122-130) depends on the order in which listeners were
registered on `window` and `document`:

- `shell-view.tsx:241-251`
- `canvas-chrome.tsx:213-219` and `428-447`
- the capture-phase listener in `agent-turn.tsx:138-151`
- `mode-switch.tsx:145-153`
- the ⌘K listener in `testing/shell.tsx:349-356`

RailMenu has no Escape handler. commands.md §6 puts the KeybindingService on the plugin root, so with several shells on
one page, one Escape or ⌘K acts on all of them.

Nothing says where focus goes when something opens or returns when it closes:

- Closing the panel or a card unmounts the focused element, and focus drops to `<body>`. Deleting a tree row or closing
  a tab does the same, because `clearNodes` and For's remove do not check `activeElement`.
- BusyOverlay claims `aria-modal`, but the rest of the shell stays operable and Cancel never receives focus.
- A Show that is closing stays interactive and focusable for about 220 ms (`control-flow.ts:72-86`).

**Change:**

1. Add a "Focus and dismissal" subsection to shell.md §2.
2. Add a layer stack on `ShellModel` or `PrismUIContext`: `pushLayer({ close, modal, returnFocus }) → undo`. A single
   Escape handler closes the top layer, attached by the KeybindingService on the root element.
3. Opening a layer focuses its heading or first control. Closing returns focus to `OpenOverlay.anchor`, the tab button
   or the canvas.
4. Modal layers (the command menu, the busy overlay, dialogs) set `inert` on the rest of the shell and trap Tab. The
   busy card focuses Cancel.
5. At runtime level, give Show and For a `focusReturn` target used when the removed range contains `activeElement`, and
   set `inert` on content that is closing.
6. Update commands.md §6 to match.

**Phase:** 6 (shell). The runtime hook affects every phase.

### M16. Contribution `when` is not reactive and comes in three forms, and render cannot reach the contribution's models

_Merges ctx-6, ux-10 and the render-context part of plan-2._

**Problem.** Visibility has three forms:

- a tab's `when(CommandContext)`, which does not exist on `ShellTab`
- a section's `when(ShellViewContext): boolean` (`shell.ts:116`)
- a toolbar item's `visible: Observable<boolean>` (`shell.ts:140`)

SectionList evaluates `when` only when the slot's registrations change (`section-list.tsx:36-39`). The Superpose
section, assembly symmetry and the volseg tab therefore never follow the data. `superposition.tsx:119` reads a snapshot.

`ShellViewContext` (`shell.ts:23-31`) has no models, and `ContributionModels` is undefined. Factory closures work, but
the declarative `models?(ctx)` field cannot be reached from `render`.

**Change:**

1. Use one reactive form for tabs, sections, overlays and toolbar items: `when?: (ctx) => Observable<boolean>`,
   subscribed once per item on activation.
2. Make the render context `ViewRenderer<ContributionViewContext<M>>` with `{ shell, models: M, focus, subject, dock }`,
   where M is the return type of `models()`.
3. Models are retired when their contribution is deactivated.

**Phase:** 6-7.

### M17. Same-id registration has no consistent rule, and undoing a replacement loses the original

_Merges ctx-5 and cmd-8._

**Problem.** `registerTabs` (`shell.ts:344-371`) and `registerOverlays` (373-382) overwrite an existing id. Undoing the
second registration deletes the id outright, so the first is never restored. `registerSections` and `registerToolbar`
never check ids, so duplicates collide in `<For key>`. CommandManager (`manager.ts:147-166`) replaces and never
restores.

commands.md §4 says "replaces" and §9 says "error". §9 is the target state, already noted as "Not yet". agent.md §6
relies on replacing Log; that works in the mock only because `withAddOnTabs` filters the list first.

**Change:**

1. Use one rule for every registry: a duplicate id throws at activation, naming both owners, unless the contribution
   declares `replaces: '<id>'`.
2. A declared replacement keeps a shadow stack, and undo restores the previous entry.
3. Apply `ui.commands.hide` and keybinding overrides in the manager, so the menu, keys, `execute` and the tool list all
   respect them.
4. Test: register A, register B, undo B; A must be restored.

**Phase:** 5-6.

### M18. Modes have no contract or spec field, and contributions cannot be scoped to a mode

_ctx-4._

**Problem.** The mocks need these things per mode, and the frozen contract has none of them:

- **Mode identity** (`MockMode`): id, title, icon, description.
- **Mode parts** (`MockModeParts`, `ui/src/testing/shell.tsx:56-86`): start tab, prefixes, canvas, overlays, playback,
  context, canvasLead, canvasAside, empty, drop handling.

`CanvasItem` (shell.md §5 line 507) is undefined; the code has `CanvasToolbarItem`. `PrismUISpec` has no `modes` or
initial mode, and `ui.contributions` is one flat list.

An extension command such as `export.geometry` opens the viewer's Export overlay. With story mode active, `openOverlay`
silently does nothing (`shell.ts:448-451`), but the command is still listed.

**Change:**

1. Add
   `PrismMode { id; title; icon; description?; start?; contributions; canvas?: { content; overlays; playback; context; lead; aside }; empty?; drop? }`.
2. Add `ui.modes` and `ui.initialMode`, and keep frame contributions separate from mode contributions.
3. Contributions may declare `modes?: readonly string[]`.
4. A command whose `opens` target is not registered has a false `when`.
5. Replace `CanvasItem` with `CanvasToolbarItem` plus named canvas-part slots.

**Phase:** 6.

### M19. Mode switching over one shared plugin is unspecified

_Merges ctx-10 and plan-11, with cmd-6 noted._

**Problem.** Both modes run on the same `PrismContext.plugin`. Switching from a loaded viewer session to story mode
leaves plugin state that the document does not hold, which breaks modes.md §3 ("nothing in story mode writes plugin
state the document does not hold"). `addScene` would capture that state.

The story player (a SnapshotStepperModel) and the viewer's Snapshots tab may share the plugin's snapshot manager.

The mock hides all this: it keeps `empty` per mode (`viewer-mode.tsx:~345`, `shell.tsx:218-221`). Status and history
entries that refer to retired nodes outlive their mode. The remembered mode choice has no store.

The agent checkpoint is fixed to a plugin snapshot (agent.md §4 lines 83-85), which is wrong in story mode.

**Change:**

1. Add `enter(ctx, from)` and `leave(ctx)` to `PrismMode`, and pick one plugin-state policy (see Q2).
2. The story player uses its document, not the shared snapshot manager.
3. Derive `empty` from the shared sources.
4. Status and history entries for retired nodes become inert.
5. Take the agent checkpoint through a `Checkpointer` port that the mode supplies, and release it on a mode switch.
6. Name the preference store (see M21).

**Phase:** 6. It does not block Phase 4.

### M20. createPrismUI has no attach phase, canvas host or dispose owner

_ctx-8, with the verifier's corrections._

**Problem.** Lifetimes have only `init()` (shell.md §6 lines 620-622), and `createPrismUI` "builds, initialises and
mounts".

Canvas3D exists only after `plugin.mountAsync` or `initViewerAsync` with a container
(`plugin/core/src/context.ts:279-322`). Nothing says:

- who calls it
- which element hosts the canvas when several elements share one context (§7)
- what `createPrismUI` returns
- who calls `dispose()`

Elements do not own their models. FullscreenModel and ThemeModel already bind late, and `canvas.ready` already covers a
late canvas at the source level.

**Change:**

1. Specify `init()` → `attach(root, canvasHost)` (calls `plugin.mountAsync`; the canvas sources start) → `detach()` →
   `dispose()`.
2. `createPrismUI` returns `{ ctx, dispose }`, and the caller is the sole owner.
3. Exactly one element per context hosts the canvas.

**Phase:** 6.

### M21. The layout spec and persistence do not match ShellModel

_Merges ctx-7 and ux-9._

**Problem:**

1. **Forgotten placement.** Unregistering a tab removes its id from `dock.tabs` (`shell.ts:360-369`), so every mode
   switch forgets where the user moved that tab.
2. **Expanded but empty dock.** `active: state.active ?? tab.id` (line 358) keeps a stale stored id, which
   `DockModel.update` (line 212) resolves to undefined: the dock is expanded with no tab.
3. **Spec types.** `layout?: Partial<ShellLayout>` forces full DockStates for all four docks, unlike
   `ShellOptions.docks`. Tab placement is specified twice. `tabs.hidden` has no API. `canvas.toolbar` (command ids)
   overlaps `layout.buttons` (toolbar item ids). `ThemeTokens` is undefined, and ThemeModel takes only theme ids.
4. **No save or restore path.** `layout` is a plain getter with no change stream, and there is no `restore`,
   `applyLayout` or `reset`. `setSize` clamps every dock to 180 px. A restore that goes through it would turn the top
   dock's 88 px into 180; the constructor does not clamp.
5. **No store.** `persistLayout: { key }` names no store, and the mode and theme choices are "remembered" with no field.

**Where:** shell.md §3 lines 277-285, §6 lines 636-649; modes.md §4 line 52; `shell.ts:81-99, 336-371, 422-424`.

**Change:**

1. Type `ui.layout` as `ShellOptions` plus `tabs?: Record<id, { dock?, order?, hidden? }>`.
2. Remove `canvas.toolbar`, or key it by toolbar item ids alongside `buttons`.
3. Drop `ThemeTokens` until ThemeModel supports custom tokens.
4. Add `ui.persist?: false | { store: PrismStore; key }`, persisting `{ layout, mode, theme }`.
5. Add `layoutChanges` (debounced), `restore(layout)` with merge rules, `reset()` and `hideTab`.
6. Keep placement for unregistered tabs, resolve a missing `active` to the first registered tab, and use a minimum size
   per axis.

**Phase:** 6.

### M22. Programmatic execution has no API, a single-run lock, and no argument contract

_Merges cmd-2 and cmd-4._

**Problem.** CommandManager holds one in-flight `#run` and refuses `open()` while it exists (`manager.ts:105, 318-327`).
`cancel()` takes no run id. `execute(id, values)` is "Not yet" (commands.md §4 lines 178-187).

agent.md expects reads to run in parallel and the person's commands to keep running during a turn (§2 line 40, §5 lines
137-138). CommandMenu closes on any `done` outcome (`command-menu.tsx:127-130`), so an agent run finishing would close a
menu the person is typing in.

`run` already accepts full `CommandArgs`. But:

- nothing specifies how programmatic arguments are validated
- `params(ctx, values)` depends on context and on `reshapes`, so one static JSON Schema is wrong or goes stale
- the tool list may change within a turn

**Change:**

1. Add `execute(id, args, { source, signal, turn? }): RunHandle { id; done: Promise<CommandResult>; cancel() }` as the
   single run path. `open`, `confirm` and keybindings all go through it.
2. Concurrency: write runs are queued, reads run concurrently. Expose `running` and per-run `cancel(runId)`.
3. The menu closes only on its own run's outcome.
4. Validate arguments by rebuilding `params(ctx, args)` until the `reshapes` fields stop changing, then run PD
   validation, format checks and `check`.
5. Schema rules: a `reshapes` field becomes a `oneOf`; context-dependent selects become strings with a run-time
   `format`, not enums.
6. The tool list is fixed per turn, and `when` is answered with its reason when a tool is called.

**Phase:** 5a.

### M23. Undo has no contract, and checkpoint release depends on command outcomes

_cmd-5, with cmd-6 noted._

**Problem.** `undoable?: boolean` is a bare flag (`types.ts:28`). The person's Undo is `state.data.undo()`, which covers
state-tree transactions only. The agent's undo is a separate snapshot, and the two interact in no defined way.

The checkpoint is released only on command outcomes (`mock-agent.ts:268-270, 317-327`). Inspector edits, 3D picking and
toolbar closures do not go through commands (commands.md §7), so "Undo turn" can discard them, which agent.md §1
forbids. Meanwhile navigation (`panel.home`) and answers release the checkpoint even though they change nothing.

**Change:**

1. Define the release trigger in core HistoryModel: any state, selection or component change not tagged with the active
   turn. Tag the transactions an agent run commits. Navigation, reads and answers do not release.
2. Decide whether `state.undo` while a checkpoint is held restores the turn or releases it, and document it.
3. An `undoable` command must commit exactly one transaction, enforced by `execute`.
4. Take checkpoints through the mode's Checkpointer (M19).

**Phase:** 5.

### M24. The command menu's assistant port is synchronous and lossy

_cmd-10._

**Problem.** `CommandMenuAssistant` (`command-menu.tsx:44-52`) has `ready(): boolean`, read once at render (line 404),
and a synchronous `propose(text)` that returns a tail string. `manager.open` re-parses that tail (lines 74-85).

A real model call is asynchronous and cancellable. Values the tail cannot express are lost, and proposed fields are not
marked as proposed. The core `CommandAssistant.propose(request, signal): Promise<AssistResult>` (`types.ts:177-179`) has
a different shape. README lines 114-117 record both ports as settled.

**Change:**

1. Make `ready` a `ReadonlyState<boolean>`.
2. Make `propose(text, signal)` async, returning `{ command, values, note? } | undefined`.
3. Add a "resolving" row state that is aborted on input or close.
4. `manager.open(id, { tail?, values?, proposed? })`.
5. Share `AssistResult` between both ports, and update the README entry.

**Phase:** 5.

### M25. The agent client is UI-bound, and the agent's ports are not contribution fields

_Merges cmd-11 and the ports part of ctx-11._

**Problem.** `MolstarUIClient` (agent.md §4 lines 59-81) mixes headless operations (tools, invoke, state, turns, events)
with UI-only focus (the open panel and section, the hovered residue). It is planned only as a UI contribution. That
conflicts with Decision 2 and with "the agent's tools in Node".

The agent reaches the UI through:

- a CommandMenuAssistant
- ShellView's `agent` slot
- `commandDone`
- replacing Log

None of these is a `UIContribution` field. The mock uses `MockShellAddOnParts` (`ui/src/testing/shell.tsx:92-99`).
`ai: { enabled }` has no place in the spec.

The boundary test is cited as `agent-boundary.test.ts` in plan.md:110, agent.md:153 and view-models.md:295. The file is
`opt-in-boundary.test.ts`.

**Change:**

1. Put `PrismClient` in core over `PrismContext`: `tools`, `invoke` via `execute`, `state`, turns through the
   Checkpointer, and `events`.
2. Build `focus()` from focus sources. `PrismUIContext` adds the panel, section and hover parts.
3. Add `commandMenuAssistant?`, `statusSlot?` and `outcomes?` to `UIContribution`.
4. Log → History uses `replaces` (M17).
5. Agent options become contribution options, not a core spec key.
6. Provide `attachAgent(prism)` for headless use.
7. Fix the test file name in all three places.

**Phase:** 5 (contracts only; see M26).

### M26. The plan builds agent wiring and prism-headless that the decisions defer

_plan-5._

**Problem:**

- Phase 5 (plan.md:169-173) lists the PD→JSON Schema generator, the state and focus accounts, the client, checkpoints,
  any-harness wiring, the MCP server and a real key store. That contradicts "agent features are mocks only".
- Phase 4 builds `@molstar/prism-headless` (plan.md:126-127), which architecture.md:37 marks "later".
- Phase 6's apps/prism includes the story document and the agent. The story document over MVS is "its own slice" and is
  unscheduled (modes.md §5). plan.md:184-186 implies an early milestone but does not state it.

**Change:**

1. Move the agent wiring to a named slice after Phase 8. Keep only the contracts core needs: M5, M22, and optionally the
   agent command fields (cmd-7).
2. Drop prism-headless from Phase 4, or define its deliverable.
3. State apps/prism's Phase 6 milestone as viewer mode on a real plugin.
4. Schedule the "story document over MVS" slice.

**Phase:** 4-6 scope.

### M27. No phase builds the dialog host, and Zenodo has three homes

_plan-6._

**Problem.** Zenodo is described three ways:

- commands.md §2.1 line 113: a panel reached by `opens`
- shell.md §4.2 and §5 and plan.md:199: a dialog
- extensions.md:29: "inline in the command menu", which commands.md §2 forbids

`CommandTarget` has only tab and overlay. `UIContribution` has no dialog kind, and `openDialog` is "Not yet"
(shell.md:285). Phase 5's "open as form" and Phase 7's Zenodo both need a dialog host.

**Change:**

1. Add a dialog host to Phase 5b or 6: `ShellModel.openDialog` and `closeDialog`, a `dialogs` contribution kind, a
   `{ dialog: id }` target, and its place in the layer stack (M15) and the keybinding scopes.
2. Make Zenodo a dialog, and fix commands.md and extensions.md to match.

**Phase:** 5-7.

### M28. Embeddable elements and embedded mode are promised but unscheduled

_Merges ux-11 and plan-7._

**Problem.** Phase 6 builds only `<msp-plugin>`. Later phases need elements nobody builds:

- `<msp-pae-plot>` in Phase 7
- `<msp-viewport>`, `<msp-snapshot-stepper>` and `<msp-params-form>` in Phase 8

`embedded` has no layout (shell.md:275), and ShellView cannot drop the rail: the card and stack insets are derived from
`RailWidth` (`shell-view.tsx:143-146`). Nothing says where rail items and cards go in embedded mode.

shell.md §7 says each element takes "a `ui` or `model`" property; MspElement has only `model`. Context sharing and
ownership across elements are undefined (see M20).

The most immediate problem: `.shell` uses `container-type: size` (`shell.module.css:9-24`), so `<msp-plugin>` inside a
container with auto height renders at 0 px.

Markdown extensions (`#key` links, command links, swatches), which mvs-stories needs, are not ported, and no phase
schedules them (shell.md §4.3).

**Change:**

1. **Phase 6:** document that the host must give the element a definite height, or set a `min-height` on it.
2. Add a work item before Phase 8 covering:
   - the embedded layout as ShellView options
   - each §7 element's property type
   - an ownership rule: the host creates the context, elements never dispose it, one element owns the canvas
3. Port the markdown extensions before mvs-stories.
4. Mark the rest of §7 as open, not frozen.

**Phase:** 6 (sizing), 7-8.

### M29. Phase 7 ports disagree with the designs, and the plan omits some extension hooks

_plan-9._

**Problem:**

- **Assembly symmetry** is a structure row in shell.md §4.1a, a Structure-tab section in extensions.md:30, and an
  inspector section in plan.md:197.
- **Model-export** is commands only in extensions.md:27, but an Export › Models section in shell.md §4.2.
- **Overlap:** `sceneEntries` and `structureRows` (shell.md §5 lines 510-513) have no stated difference.
- **Missing from Phase 7:** anvil's Membrane Orientation, the `leftPanelTabName` writers (sb-ncbr, volseg, g3d), and the
  volume-streaming editor.

**Change:**

1. Pick one form per extension and update extensions.md, shell.md and Phase 7 together: assembly symmetry as a structure
   row; model-export as an Export section plus a command.
2. Merge `sceneEntries` into `structureRows`, or define the difference.
3. Add the missing ports.
4. Specify how an action marks its task as an overlay task (`StatusTask.overlay` exists; the source-side flag does not).

**Phase:** 7. The contribution overlap is part of the freeze.

## Should fix

**rt-3. Re-entrant regions and For.** If a plain BehaviorSubject re-emits synchronously while its region renders,
`insertRegion` (`render.ts:184-189`) inserts stale content built in a disposed scope. For's reconcile has no guard.
Insert only if `scope === next`, add a `reconciling` flag to For, and state the guarantee in rendering.md §4.2. Low risk
while sources are built on `observableState`, which queues re-entrant sets.

**rt-4. JSX emitted from streams.** `map(v => <Comp/>)` builds its bindings outside the region's scope. On later
emissions, `requireScope` throws inside the subscriber. This only affects JSX that creates bindings. Add a §4.4 rule and
lint ("a stream child emits text or `() => JSX`"), a dev warning, and a clearer `requireScope` message.

**rt-5. MspElement model lifetime.** A model the element creates is never retired (`element.ts:50-56, 84-111`).
`elementView` sets `ctx` while MspElement reads `model` (`view.ts:27-34`). In rendering.md §5, retire the owned model
when it is replaced or the element is finally disposed, and re-create it on reconnect. Default `elementView` to `model`
for MspElements.

**rt-6. `batch` across states.** combineBehaviors projects from combineLatest's raw values
(`combine-behaviors.ts:48-51`), so two states set in one batch emit an intermediate pair. There is no visible frame, but
views render twice. Read `.value` in the projection, and correct view-models.md §3 and spike.md to say that batch
dedupes per state.

**rt-8. For fallback during time slicing.** The fallback is shown when `order.length === 0` (`control-flow.ts:354-358`),
so a non-empty list whose entries are all deferred flashes its empty state. Decide from `items.length` and add a test
with slicing on. (The duplicate-key variant does not happen.)

**rt-9. Switch edge cases.** `combineLatest([])` never emits, so a Switch with no cases never shows its fallback. Nested
arrays of cases are dropped (`control-flow.ts:127-140`). Flatten children deeply and use `of({ index: -1 })` when there
are no cases.

**rt-10. Shadow-root embedding.** Styles are adopted into `ownerDocument`, the `.msp-root` check stops at shadow
boundaries and the class is never removed, and delegation uses the global `document` (`element.ts:65-66`,
`render.ts:310`). State in rendering.md §5 whether shadow-hosted embedding is supported (Q7). If it is, use
`getRootNode()`.

**rt-11. Time slice length.** `sliceEnd` is reset by a microtask (`control-flow.ts:177-183`), so promise chains within
one task each get a new slice. Reset at a task boundary, as `flushWaiting` does, and re-measure with sources driven by
async updates in Phase 4.

**vm-10. Label format.** HoverLabelModel and MeasurementPoint promise plain text, but Mol* labels are HTML or non-string
values (`loci-label.ts:14-15`), and they include text from files. Freeze the conversion in the source contract: strip
tags, or structured spans from an allow-list. Drop labels that are not strings. This also covers plan-8's hover-label
item.

**vm-11. Object URLs and timestamps.** Snapshot images and previews need object URLs that nobody revokes. Timestamps are
stored as pre-formatted "2 min ago" strings (`snapshots.ts`, `image-export.ts:26-27`, `home.ts:37-38`). The model owns
and revokes URLs per asset or preview generation. Timestamps are epoch milliseconds, and views format them.

**vm-12. Unspecified Phase 4 panels.** Component options and add-component, picking mode, canvas settings, animation and
focus have no model interfaces. List them in plan.md as not frozen, or add minimal interfaces: ComponentsModel,
PickingModel, AnimationModel, FocusModel.

**cmd-7. Agent command fields.** agent.md marks `read` and `approval` as planned. As booleans they cannot carry an
approval reason that depends on the arguments ("Fetches data from www.ebi.ac.uk"). Plan `effect: 'read' | 'write'`,
`approval?(args, ctx): string | undefined` and `audience?`, and require JSON-serializable `CommandResult.data`.

**cmd-9. Keybinding scope.** `CommandDefinition` has `keys?` and `shortcut?` but no scope (`types.ts:72-73`). The spec
override allows one chord per command, while several are planned. The turn card's capture-phase Esc listener pre-empts
everything. Replace both fields with `keybindings?: { chord; scope? }[]`, make the override `chords: KeyChord[] | null`,
and make the KeybindingService the only root listener (see M15).

**cmd-12. Command menu and turn card accessibility.** The menu has:

- no `aria-activedescendant` and no option ids
- a constant `aria-expanded="true"`
- no `aria-modal`
- no focus restore when it closes (`command-menu.tsx`)

The turn card's whole section is a polite live region, so streamed text is announced. A pending approval does not move
focus. Specify these in commands.md §5 and agent.md §5:

- the menu: option ids, active descendant, `aria-modal`, focus restore
- the turn card: announce state changes and the final answer only; announce a pending approval assertively, focus it,
  and make Deny the default

**ux-1. A card hides the canvas stack when the panel is collapsed.** `rightInset` counts only `right.active`
(`shell-view.tsx:119-125`), so a card opened over a collapsed panel covers the stack and the end of the selection bar.
Use `f && (active || overlay)` in `rightInset` only, and add a test.

**ux-2. Theme step tokens leak into nested themes.** The light block never sets `--msp-step-*`, and `auto` sets nothing
outside its dark media query (`theme.css`). The light canvas chrome under a dark UI therefore gets dark-theme shades.
Every theme block sets every base token, and `auto` gets the light values outside the query.

**ux-3. Theming and overrides.**

1. `color-scheme: light` on `:root` changes the host page. Move it off `:root`.
2. Make `element.ts:26` and styling.md §3 agree on where the tokens are set.
3. Reword "override any token on any element" to "base tokens on an element with `data-msp-theme`".
4. State now that generated class names are private.

`ThemeTokens` and the `data-msp-part` names stay deferred, as plan.md:113 says.

**ux-6. Canvas stack menus.** ToolButton opens on focus and closes 160 ms after blur, and has no `aria-expanded`.
ToolMenu and RailMenu claim `role=menu` without the menu keyboard pattern (`canvas-chrome.tsx:112-353`). The features
themselves are reachable through commands. Add `aria-expanded` and `aria-controls`, open and focus the menu on Enter or
ArrowLeft, move between items with arrow keys, and return focus on Escape. Otherwise drop `role=menu`.

**ux-7. Rail and splitter semantics.** The rail has three `role=tablist` navs for one panel. A tab collapses the panel
when clicked again. The tabpanel is not labelled by its tab. The splitter is not focusable and handles only horizontal
resizing, with a hard-coded 300/340 reset that is wrong for the top and bottom docks (`dock-view.tsx`). Pick the rail
semantics (Q8). Make the splitter focusable with `aria-valuenow` and arrow keys, take the reset from `DefaultSizes`, and
take the axis from the dock id.

**ux-8. Dock states the view ignores.** `setFloating` and `setPinned` accept any dock, but the view renders floating
only for right and top, and pinned only for right. The compact grid has no `left` area. `embedded` has no layout. The
compact theme switch renders its `face`, not the static icon §2 describes. Narrow the model to what is built, specify
compact behaviour for every dock, and align the theme-switch sentence with the code.

**ux-12. Status live region.** The whole status column is `aria-live=polite`, including the task pills and the agent
slot. Errors get no `role=alert` (`canvas-chrome.tsx:462`). Put only the message text in the live region, use
`role=alert` for errors, and keep tasks and the agent card outside it. Name the ✕ buttons after their task or chip.

**cmd-6. Agent checkpoint in story mode.** This is handled under M19 and M23 through a Checkpointer port. Phase 5 should
not hard-code a PluginState snapshot.

**plan-8. Parity gaps.** SceneNodeKind has no particles kind, although a particles source is designed. The classic Help
tab's mouse and behavior binding tables have no home, and the Help view that commands.md §6 promises is unscheduled. Add
particles to SceneNodeKind and its mock, and schedule a Help view with the KeybindingService in Phase 6 (Q9).

**plan-10. Testing and acceptance.**

1. Add a parity matrix (current-ui §5 → Prism home → test) as Phase 6's exit criterion.
2. Add shared contract tests that run against both the mocks and the real models.
3. Add a perf bench over `createPrism`.
4. Make the a11y criteria executable checks.
5. Add `styles:check` and `icons:check` to CI.
6. Name a host package for the integration stories.
7. Resolve plan.md:77 (Storybook `play` tests) against architecture.md §9 (Jest, not `play`); see Q10.

**plan-12. API ledger.** `display.triggerLabel` is in core (`core/src/state/action.ts:58-61, 105`;
`transformer.ts:167-170, 282`) but not in changelog.md. shell.md §6 says removal candidates are recorded in the ledger,
but they are not. Add the changelog entry. Make a changelog entry part of each compatible follow-up's definition of
done. Reword shell.md §6 to "will be recorded in Phase 4".

## Open questions for the lead

1. **Mode model lifetime (B1, M19).** Should a mode's models live from first activation for the life of the context, or
   be retired when the mode is left? **Recommendation:** retire on leave, the same lifetime as the mode's contributions.
   Their state lives in the plugin, and keeping them risks leaks into the other mode. Persist only preferences.
2. **Plugin state on a mode switch (M19).** The options are: confirm and clear (offering to save a snapshot); snapshot
   and restore per mode; or import viewer state as a story scene. **Recommendation:** confirm and clear with a snapshot
   offer, and load story mode over a cleared plugin. It is the only option that keeps stories reproducible.
3. **Same-id rule (M17).** Should a duplicate throw, or should the last registration win with a stack?
   **Recommendation:** throw by default; an explicit `replaces` keeps a shadow stack so Log → History works and is
   reversible.
4. **Async model methods (M9).** `Promise<ActionResult>`, or `void` plus status and busy states? **Recommendation:**
   `Promise<ActionResult>` that never rejects, also posted to status, plus `busy` states for long operations.
5. **Optimism for flags (M12).** `updateCellState` is synchronous. Do visibility toggles need optimistic projection?
   **Recommendation:** no. Use optimism only for parameter edits that run a transform.
6. **Agent and headless scope (M26).** Should the agent wiring leave Phase 5, and prism-headless leave Phase 4?
   **Recommendation:** yes to both. Phase 5 freezes only the run, history, execute and checkpoint-port contracts. The
   `createPrism` Node test lives in core.
7. **Shadow-root embedding (rt-10).** Supported in v6? **Recommendation:** not supported. Document light DOM only.
   Revisit when a wrapper library needs it.
8. **Rail semantics (ux-7).** Toggle buttons or a single tablist? **Recommendation:** toggle buttons with
   `aria-expanded` and `aria-controls`. That matches "click again collapses" and "a card covers the panel".
9. **Help (plan-8).** Where do mouse and behavior bindings and keybinding conflicts live? **Recommendation:** a Help
   overlay contribution, built in Phase 6 with the KeybindingService and also reachable from the command menu's `?`.
10. **Interaction tests (plan-10).** Jest or Storybook `play`? **Recommendation:** Jest is the gate for behaviour.
    `play` tests stay as smoke tests, run by the Storybook test runner in CI.
11. **Unbuilt Phase 4 models (vm-12).** Freeze minimal interfaces now, or list them as not frozen? **Recommendation:**
    list them as not frozen and add interfaces when each panel is built.

## Checked and fine

These held up against the code:

- **The settled decisions** stand as written.
- **observableState** queues re-entrant sets correctly.
- **FullscreenModel and ThemeModel** bind late and do not need the root at construction.
- **The `canvas` source** already models a late canvas through `ready`.
- **view-models.md §4.3** already retires keyed-list nodes when their key disappears.
- **`run` signature:** it accepts full `CommandArgs`, so programmatic execution needs a contract but no signature
  change.
- **StatusTask** already has an `overlay` flag.
- **In-card widgets** that handle Escape on their own element (Select) already run before the document listeners.
- **ShellModel's constructor** does not clamp spec sizes.
- **PluginLayout's body mutation** is avoided as long as Prism leaves `spec.layout` unset.
- **MVS load errors** are logged as well as toasted, so a log-fed status does not lose them.
- **The opt-in boundary test** exists and covers the agent package. Only the file name cited in the documents is wrong.
