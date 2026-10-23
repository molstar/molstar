# The viewer's API

`viewer.api`: operations by area, the one way state changes, called alike by commands, tools, scripts, the URL and view
models. The plugin's `managers` and `builders` done again: one run path, targets, results, rxjs states. Every operation
with its flags and params: [operations.md](operations.md) (generated). The code is the truth for how (`src/api/`).

## 1. An operation

`OperationDefinition` (`api/area.ts`): `id` `<area>.<verb>`, `title`, `params` (PD, or a function of the context for
registry-driven ones), flags `read`, `history: false` (view state), `undoable: false`, `quiet`, `record: false`,
`approval` (`network` | `upload` | `destructive`; destructive wins), `addresses(values)` (what a network call reaches;
absent: the agent asks), `modules`, `when` (true or the reason), `exists` (identity), `project` (optimistic edit),
`surface` (§10), `structureMenu` (an entry in a structure row's "…" menu that runs it for that structure where `applies`
holds: model-archive's PAE fetch for an AlphaFold DB entry), `aliases`, `run(values, ctx)`.

- **Params are PD.** One definition gives the command's form, the tool's flat input and validation (`checkValues`,
  `api/run.ts`). Targets are text fields with the `target` format ([targets.md](targets.md)).
- **Desired state.** An operation says what should be; the same call twice changes nothing the second time
  (`ApiResult.already`: `unchanged`, no history step), so a retried step is safe.
- **Absent versus given.** A value the operation can keep is optional (`OptionalParam`, `api/params.ts`): absent is
  kept, given is set. No sentinel means "keep". A call that gives nothing to change fails `invalid`, naming what it
  takes.
- **Switches take the state wanted**, never a toggle (commands flip it: `toggle` in `commands/operations.ts`). With
  settings: a union `{ name: 'off' } | { name: 'on', params }` (`SwitchParam`), settings partial; without: boolean `on`;
  visibility is `visible`. A choice with its own settings is a union too (`export.image`'s `size`).
- **Inputs: one loading vocabulary.** Content-reading operations (`load.data`, `load.volume`, `load.mvs`, `load.svs`,
  `session.open`) take exactly one of `url`, `data` (text, `'base64,…'`, bytes) and `file`; none or several fail
  `invalid`. `load.url` takes a url only, `load.files` several files. `format: 'auto'` asks Mol\*'s format registry by
  extension, else by content; one it cannot tell fails naming the formats. `binary` follows the extension, else the
  format. Modules are the input's, per call; a url asks as `load.url` does (`api/input.ts`, `api/formats.ts`).
- **Nodes**: operations on state tree cells take `node`, a CSS-like selector ([nodes.md](nodes.md)); acting on several
  is one step, acting on one fails `ambiguous` with the matches.
- **Words**: `list` reads a collection, `get` one thing or a state; set modes are `add`, `remove`, `intersect` (and
  `set`); a read names a structure `{ ref, label }`; an area's subscribable state is a noun (`structure.loaded`). One
  state, one operation (`trajectory.frame`, `canvas.camera`). A renamed operation keeps its old id in `aliases` until
  the release (listed nowhere). An empty `remove` acts only when there is one, else `ambiguous`; `all: true` clears.
- **Identity for operations that create things** (README decision 30): `exists` says what makes it the same
  (`load.entry`: source and id; `load.url`: the URL; `component.add`: structure and target; `measure.*`: kind and
  targets in order; `snapshot.add`: none). If it exists the call does nothing and says so, unless `copy: true`.
- **Failure is all or nothing**: a failed or aborted run restores its step's checkpoint ([sequence.md](sequence.md)
  §7.3). **Never throws**: failures are results `{ ok: false, code, message, hint? }`.
- **One run path** (`api/run.ts`): `api.run(id, values, { signal, copy })` checks params, modules, `when` and `exists`,
  opens a history step for an undoable operation (none when `unchanged`), applies `project`, runs, records the run (id,
  values, origin, result, time, step or group). A call through `ctx.api` is part of the caller's step.
- **Origins** (`RunOrigin`): `person`, `agent` (with its turn), `script`, `url`, `link` (the sequence,
  [sequence.md](sequence.md) §4).
- **Plugin boundary**: only `api/`, `sources/` and features' `operations.ts` and `source.ts` import the plugin as values
  (`import type` anywhere; `_test/boundary.test.ts`).

## 2. Areas

One object per area; its name is the area's id (`viewer.api.component.add` is `component.add`): one name in the API, the
command menu, the tools, the run log and the docs.

| `viewer.api.…`                 | Over (the plugin)                                                   | Brought by                    |
| ------------------------------ | ------------------------------------------------------------------- | ----------------------------- |
| `load`                         | `builders.data`, the loaders, state actions                         | core; `loaders`               |
| `session`                      | sessions, the remote-states client                                  | core; `remote-states`         |
| `structure`, `contents`        | `managers.structure.hierarchy`, hierarchy presets; entities, chains | `structure`                   |
| `component`, `representation`  | `managers.structure.component`, themes, the style layers            | `structure`                   |
| `styles`, `annotations`        | representation presets and looks; what extensions can apply         | `structure`; `annotations`    |
| `volume`                       | `managers.volume`                                                   | `volume`                      |
| `particles`                    | `managers.particles`                                                | `particles` (lazy)            |
| `trajectory`                   | the model-index transforms; a simulation box's display (`periodic`) | `trajectory`                  |
| `select`, `focus`, `highlight` | `managers.structure.selection`, `.focus`, `managers.interactivity`  | core                          |
| `measure`                      | `managers.structure.measurement`                                    | `measurements`                |
| `superpose`                    | the superposition algorithms                                        | `superposition` (lazy)        |
| `snapshot`, `animation`        | `managers.snapshot` (its list only), `managers.animation`           | `snapshots`; `animation`      |
| `camera`, `canvas`             | `managers.camera`, `canvas3d` props, illumination, XR               | core; `canvas-settings`       |
| `export`                       | `helpers.viewportScreenshot`; extensions' exports                   | `export-image`; extensions    |
| `sequence`                     | the Seq\* viewer ([sequence.md](sequence.md))                       | core (with a sequence viewer) |
| `state`                        | `state.data` (the tree, applicable actions), by `node`              | `state-tree` (lazy)           |
| `settings`, `history`, `view`  | config, behavior params, options; prism-ui's history, navigation    | core                          |

- **Core areas** (`CoreApi`) are in every viewer. Built-in features' areas and their operations in core areas
  (`load.entry`, `canvas.lighting`, `session.save`, …) are typed as present (`BuiltInApis`, `BuiltInOperationFeatures`);
  with the feature left out, reading one throws naming it. Code that works either way asks `api.area(id)` or
  `api.describe(id)`. Lazy and gated areas and extensions' are optional (`Partial<ExtensionApis>`).
- **Hiding contents** (README decision 18): `contents.visibility` makes each component hold what it held minus the
  hidden entities and chains; the hidden set is custom state (`prism-viewer.contents`), and `structure.preset` and
  `component.add` apply it in their step. A component left empty is a null cell, restored when shown again.
- **Registries are arguments, not operations**: one `structure.preset` takes a preset id, `representation.color` a theme
  name; completions and the tool's enum come from the registry.
- **Events** for embedding hosts (v5's `structureInteractivity`): `viewer.events.hover`, `click`, with the target under
  the pointer.

## 3. The shape every area shares

- `ApiArea` (`api/area.ts`): `id`, made once per viewer from an `ApiAreaContext`; operations declared as fields with
  `this.operation(def)` (the method and its definition are one: `fn.definition`); states are `ReadonlyState`s from the
  sources; `contribute(name, def)` for built-in features.
- **Made once, bound by facades**: `viewer.api`, `api.as(origin)`, a group's API and `ctx.api` rebind the same area's
  operations (§5); its states are the same objects through every facade.
- **Values are partial**: a missing value is the param's default.
- **Reads** are the states (views, scripts) and `read: true` operations, which return text and data for commands and
  tools.

## 4. The types around it

In `api/area.ts`: `ApiAreaContext` (plugin, sources, targets, modules, managers, operations, extensions, journal; no
origin), `ApiContext` (adds `api`, `origin`, `signal`, `copy`), `ApiErrorCode` (`invalid`, `blocked`, `unavailable`,
`not-found`, `ambiguous`, `confirm`, `aborted`, `failed`), `ApiModuleKey` (`canvas`, `files`, `network`, `storage`,
`page`; what the host offers, `ApiModules.has`), `Projection`, `ViewerApi` (the areas plus `run`, `as`, `group`, `area`,
`describe`, `definitions`, `running`, `events.started`/`finished`).

## 5. Origins, grouping, scripts

- **Bound by origin**: `viewer.api` runs as `script`; view models and commands use `api.as({ kind: 'person' })`, tools
  `{ kind: 'agent', turn }`, the loaders' URL calls `url`, the sequence `link`. The history and the record read it.
- **Part of a run**: calls through `ctx.api` join the caller's step (no own step, no wait on the history's queue). A
  call from another origin waits and is its own step.
- **`api.group(label, work, { signal?, onFailure? })`**: one history step and one record group; fails as a whole,
  restoring its checkpoint, when the work throws, returns a failure, or a call inside failed (`onFailure: 'continue'`
  undoes just that call). The result carries `runs` and the union of `changed`. In an agent's turn it joins the turn.
- **Scripts import loaders as needed**: v5's `load*` are plain functions over the plugin
  (`@molstar/plugin-extension/loaders`); there is no `load` facade pulling them all in.
  `viewer.run(label, (plugin, api) => …)` runs raw plugin work as one opaque step; its `api` calls join it. A plugin
  change outside the API still shows but is no step, and the steps before it can no longer be undone
  ([sequence.md](sequence.md) §7.3).
- **Declarative scenes** are MolViewSpec, loaded with `load.mvs` in one step; no second builder language.

## 6. Extensions and what is public

- **Public surface**: `@molstar/prism-viewer` (`createPrismViewer`, `DefaultPrismViewerSpec`, the types, `ApiResult`,
  `node`); `/api` a curated list for scripts and extensions (`src/api.ts`: the area base class, result and operation
  types, `checkValues`, node selectors; the runner stays inside); other subpaths one by one in `package.json` exports.
- **`viewer.sources` is read-only**: every source a state; the marks' writers are the API's (`highlight.mark`,
  `unmark`).
- **An extension's areas** come from `ViewerExtension.api` (factories by area id, each made once): operation ids start
  with the area's id (checked). Only built-in features (`builtIn: true`) add operations to core areas (`contributes`),
  so core areas do not grow with Mol\*'s extensions ([viewer-spec.md](viewer-spec.md) §3).
- **Typed when imported**: an extension augments `ExtensionApis` (`declare module '@molstar/prism-viewer/api'`).
- **Lazy**: stub areas from `stub.operations`; a call loads the chunk, the real area replaces the stub in place, the
  call is replayed (a failed load is the result).
- **Left out**: the area is absent; commands and tools over it say why they cannot run.

## 7. Code: the agent and scripts writing against the API

Decided, not built: a `viewer-code` skill (`run_viewer_code { code }`) runs TypeScript against `ViewerApi` in a sandbox
(an `async function executor(api)`; limits on time, calls, output), in a group (§5) bound to the agent's turn, all or
nothing; the typings the model reads are generated from the API's definitions, and also serve as the people's reference
(`apps/docs`).

## 8. What is promised

- **Core areas are the package's public contract**: ids, params and results follow semver once released (operation ids
  are also what the run log, kept records and skills name). Extensions' areas carry no promise (README decision 13).
- **`viewer.plugin` stays reachable** and is not part of the promise.
- **v5's `Viewer` methods** map onto areas (`load*` → `load`, `structureInteractivity` → `highlight`, `select`,
  `events`; `handleResize` → nothing, the shell sizes the canvas).

## 9. Commands over the operations

- **One command per operation by default** (`commands/operations.ts`, `operationCommand`): its id, title and params;
  form from the params, inline sentence from the flat fields ([commands.md](../../ui/_spec/commands.md)).
- **Generated commands** where a list reads better: one per preset, quick style, selection query, applicable state
  action, each the same operation with that argument fixed. The Select group (`commands/select.ts`, its queries from
  `sources/select-queries.ts`): `select.set` per query of the plugin's registry (by category) and per entity, chain,
  non-standard residue and element of the loaded structures (made again when they change), with the group's shared mode;
  a highlighted row lights its atoms (`lociPreviewSink`) and is counted (ui's commands.md §5).
- **Opening commands** (`<feature>.open`, `panel.<tab>`) run `view.open`: the way into a custom UI, whose actions are
  the area's operations.
- **Keybindings and toolbar buttons name commands**, never operations, so hiding a command hides them.
- **Own forms** (ui's commands.md §5.1; `views/command-forms/`, `sources/command-forms.ts`, given to the menu by the
  shell): Show as… (`representation.set`) is the types the component supports as icon rows (those it shows marked, one
  it does not picked), the picked type's size, quality and opacity under them (its hidden `params`); Color by theme…
  (`representation.color-theme`) is the themes that apply, by the registry's categories, painted with their legends on
  the structure, and the uniform theme's color (its hidden `color`). A click on a theme or a swatch applies it (nothing
  more to set); a click on a type picks it (its options follow), Enter or a double click adds it. The inspector's Add
  representation… and Set coloring… open them filled in. Why: a long registry list reads better as what it draws than as
  words, and the sentence keeps the keys.
- **The query language shows only for a query** (`operationCommand`): a command with `target` and `language` shows the
  language in its sentence only while the target text reads as a query (not a keyword or an address list), optional
  (empty: the default); it follows the target (`reshapes`). Why: "ligand · PyMOL" read as if the coloring were PyMOL's.
- **Own forms read as sentences**: the subject, `in` its structure, then the choice (`as` a type, `by` a theme), by the
  command's `order` and `fields`.
- **A select whose default is the empty choice is asked for** (`formDefaults`): its "None" would fill the field and the
  form open on the next one (the theme's form opened on the query language).

## 10. Coverage

Every operation has a command, a tool, or both, or `surface: 'ui-only' | 'agent-only'` with a reason
(`src/_test/coverage.test.ts`). Each area's file lists its operations in its header comment.

## 11. One API under commands and tools

- **The API is the only way state changes.** Commands (people: forms, previews, keys) and tools (the agent: flat inputs,
  in skills, [agent.md](agent.md)) are two surfaces over it (README decision 8). What a view shows comes from read-only
  states fed by the sources, never from a call's return value.
- **Managers and views**: `PrismViewer.managers` (`commands`, `transactions`, `history`, `navigation`, `viewState`,
  `notifications`; generic, DOM-free, prism-ui's) and `views` (root view models: `shell`, `scene`, …), beside `plugin`.
  `transactions` are optimistic edits: shown at once, absorbed when the source reports them, else rolled back.
- **Undo is the viewer's, one history**: every undoable step (a person's call, an agent's whole turn) is one entry with
  a checkpoint taken before it (`plugin.state.getSnapshot` with chosen flags plus custom state, never the snapshot
  manager); Undo restores it ([sequence.md](sequence.md) §7). Mol\*'s state history is off.
- **Three layers of state**: the engines' (Mol\*'s tree, camera, canvas, selection, focus; Seq\*'s document), the
  viewer's own (what the sequence follows, the agent's turn), the UI's (open panel, docks). Ephemeral UI (hover, scroll,
  an unsaved form) is not state.
- **Headless**: `@molstar/prism-viewer/headless`, the same API, tools and account without a DOM; an operation whose
  module the host lacks is listed and says why it cannot run. prism-ui's `reactive/`, `commands/` and `forms/` model are
  DOM-free, so no separate core package.
- **Evals** run an agent over the headless viewer and check the state account ([agent.md](agent.md) §6).
- **Not Stories**: no closed document language; the "document" is the plugin's state tree plus the sequence's
  SeqViewSpec. Operations compile to state updates and sequence patches; undo restores snapshots, not inverses.

## 12. The camera as a view, and viewers that follow each other

`viewer.camera` (`sources/camera.ts`; not an operation: a pointer-rate stream and a cut belong neither in the run path
nor in the record or the history, and `camera.get`/`camera.set` are rounded, a result, a transition: for people and the
agent):

- `snapshot()`: the view now, full precision: `{ mode, fov (rad), position, target, up, radius }`; none without a
  canvas.
- `changed`: an `Observable` of the view after a frame drew a different one (the person dragging, a call, a transition's
  frames). Not replaying. It speaks only when the view differs from the one last announced or set, by value
  (`sameView`), so a view set `silent` is never announced.
- `set(snapshot, { silent? })`: cuts the camera to the view (no transition: one under way ends there; the radius is
  clamped to the scene's) and asks for a draw; false without a canvas. `silent` marks the view the camera took (after
  the clamp) the last announced, so a host that mirrors one camera onto another does not hear its own setting back: no
  echo loop, however many viewers follow each other, and a follower's clamp never comes back to its leader.
- `syncCameras(viewers, from?)` (`@molstar/prism-viewer`): each viewer's `changed` sets the others silently; `from`
  applies its view to the rest at once; returns the unlink. The radius is part of the view so the clipping and fog
  match; a viewer whose scene is smaller clamps it (silently: the leader keeps its own).
