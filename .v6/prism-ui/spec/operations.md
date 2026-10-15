# The viewer's API: operations, and the commands over them

The one way the viewer's state changes (README §3.3): operations, area by area, with the commands people use over them.
The agent's tools over the same operations are [agent.md](agent.md); what each area reads from the plugin is
[state.md](state.md); what it covers of v5 is [parity.md](parity.md). Status: **accepted** in shape (2026-10-15); the
operation lists are a first pass, filled in as each area is built.

## 1. An operation

```ts
interface Operation<P extends PD.Params = PD.Params, R = unknown> {
  id: string;                                    // '<area>.<verb>': 'component.add', 'canvas.background'
  title: string;                                 // for people: the command's title, the history entry
  params: P | ((ctx: ApiAreaContext) => P);      // one definition: the command's form and the tool's schema derive from it
  read?: true;                                   // changes nothing: no history step, may run alongside others
  history?: false;                               // changes view state, not content: no history step (select, focus, camera, highlight, view)
  undoable?: false;                              // a step the history cannot restore (an upload); default true
  approval?: 'network' | 'upload' | 'destructive';   // asks before it runs, from the agent (agent.md §4)
  modules?: readonly ApiModuleKey[];             // what it needs of the host (canvas, files, network); absent: says why
  when?(ctx: ApiAreaContext): true | string;     // can it run now; the string is the reason (no origin: commands read it)
  exists?(values, ctx: ApiContext): string | undefined;  // identity: what it would create exists, said as the result (decision 30)
  project?(values: PD.ValuesFor<P>): Projection; // an optimistic edit (managers.transactions), when the result is known
  surface?: { kind: 'ui-only' | 'agent-only'; reason: string };   // §4
  run(values: PD.ValuesFor<P>, ctx: ApiContext): Promise<ApiResult<R>>;   // ctx.signal, ctx.api (part of this run)
}

type ApiResult<R = unknown> =
  | { ok: true; did: string; changed: Changed; value?: R; warnings?: string[]; next?: string; unchanged?: true }
  | { ok: false; code: ApiErrorCode; message: string; hint?: string };   // api.md §4

interface Changed { nodes?: string[]; selection?: boolean; view?: boolean; sequence?: boolean }   // what a reader re-reads
```

- **Params are PD.** The same definition gives the command's form (prism-ui's `ParamsFormModel`), the tool's input
  (`paramsShape`, flat fields for small models) and validation. Targets are text fields with the `target` format
  ([targets.md](targets.md)).
- **Desired state.** An operation says what should be, not what to do: `representation.set` with the same values twice
  changes nothing the second time (`ApiResult.already("… already so.")`: `unchanged`, no history step), so a retried
  step is safe.
- **Identity for operations that create things** (decision 30): each says what makes the thing it creates the same
  (`load.entry`: source and id; `load.url`: the URL; `component.add`: structure and target; `measure.*`: the targets in
  order; `snapshot.add`: none, every call adds). When it exists, the call does nothing and says so ("1TQN is already
  loaded."), unless `copy: true` asks for another. Where the person's intent is unclear, the agent asks (`ask`,
  agent.md §3).
- **Failure is all or nothing**: a failed or aborted run restores the checkpoint its step took (sequence.md §7.3), so a
  loader that failed after its download leaves nothing behind.
- **Never throws.** Failures are results with a code, a message for the person or the model, and a hint naming what to
  do instead.
- **One run path**: `api.run(id, values, { origin, signal })`. It checks the params, the modules and `when`, says
  "already so" for what `exists` names, opens a history step for an undoable operation (the checkpoint before it,
  sequence.md §7; none when the result is `unchanged`), applies `project` through `managers.transactions`, runs, and
  records the run (id, values, origin, result, time, its step or group). It never throws. A call an operation makes
  through `ctx.api` is part of its step. Commands, tools, the JS API, URL actions and the
  sequence's actions all call it; a view model calls the operation's typed function, which calls it.
- **Origins**: `person`, `agent` (with its turn), `script`, `url`, `link` (the sequence, sequence.md §4).
- **Where**: `packages/prism/viewer/src/api/<area>.ts`, and a feature's `features/<name>/operations.ts` (API code by
  role: what it adds to core areas or its own); with `sources/`, the only code with value imports of the plugin
  (`import type` is allowed anywhere; an extension's `plugin` part imports its own behaviors and entries;
  `_test/boundary.test.ts` checks it).
- **The types around it** (`ApiContext`, `ApiModuleKey`, `Projection`, `ViewerApi`) are defined in [api.md](api.md),
  which turns these areas into the viewer's API objects.

## 2. Areas

| Area           | Operations (first pass)                                                                                                     | Reads                         | From                |
| -------------- | --------------------------------------------------------------------------------------------------------------------------- | ----------------------------- | ------------------- |
| `load`         | core: `files`, `action` (a registry state action by id: Home's blocks); the loaders feature: `entry` (`pdb`, `pdb-ihm`, `afdb`, `emdb`, `model-archive`, id), `url` (url, format, binary), `mvs`, `svs`, `snapshot` (url or remote id), `trajectory`, `volume` | — | core; loaders |
| `session`      | `open` (`.molx`, `.molj`), `save`, `clear`; `remote.list`, `remote.open`, `remote.upload`, `remote.remove`                  | `remote.list`                 | core; remote feature |
| `structure`    | `preset` (hierarchy and representation presets by id), `model`, `assembly`, `visibility`, `remove` (the whole entry)        | `get` (the scene outline)     | core                |
| `contents`     | `visibility` (entity or chain: left out of every component, below)                                                          | `get` (entities, chains)      | core                |
| `component`    | `add` (target, label), `modify` (target, mode `union`/`subtract`/`intersect`), `remove`, `rename`, `visibility`, `options` (hydrogens, quality, lighting, …) | — | core |
| `representation` | `set` (on a component: type and params), `update` (a representation's params, the inspector), `visibility` (a row's eye), `remove`, `color` (target, color or theme), `size`, `transparency` (target, alpha: the overlays), `emissive`, `material`, `clipping` (target: excluded from the clip objects, Mol*'s clipping layer), `reset` (a layer) | — | core |
| `annotations`  | `apply` (what the structure carries: quality scores, validation, charges; extensions register what they can apply) | `list` | core; extensions |
| `styles`       | `apply` (quick style by id: `default`, `cartoon`, `spacefill`, `surface` (representation presets); `plain`, `illustrative` (looks)) | —                     | core                |
| `volume`       | `isovalue`, `representation`, `update` (ui-only), `visibility`, `remove`, `streaming`                                       | `get`                         | core                |
| `particles`    | `representation`, `update` (ui-only), `visibility`, `remove`                                                                | `get`                         | core                |
| `select`       | `set` (target, mode `set`/`add`/`remove`/`intersect`), `all`, `clear`, `invert`, `around` (radius), `expand` (to residues, chains, …); `mode`, `granularity` | `get` (what is selected) | core           |
| `focus`        | `set` (target), `clear`                                                                                                     | `get`                         | core                |
| `highlight`    | `set` (target), `clear` (for embedding hosts and the link; no history step)                                                 | —                             | core                |
| `measure`      | `distance`, `angle`, `dihedral`, `label`, `orientation`, `plane` (targets), `remove`, `options`                             | `get`                         | measurements        |
| `superpose`    | `chains`, `tm-align`, `ligands`, `atoms`, `db` (targets)                                                                    | —                             | superposition       |
| `snapshot`     | `add`, `apply`, `update`, `replace`, `move`, `remove`, `next`, `previous`, `play`, `stop`                                   | `list`                        | snapshots           |
| `animation`    | `play` (by name and params), `stop`, `wiggle`                                                                               | —                             | animation           |
| `trajectory`   | `frame`, `play`, `stop`                                                                                                     | `get`                         | core                |
| `camera`       | `reset`, `zoom`, `orient`, `frame` (target)                                                                                 | `get`                         | core                |
| `canvas`       | `background`; canvas settings: `settings` (ui-only), `lighting`, `illumination`, `camera`, `clipping`; `xr`                 | `get`                         | core; canvas settings |
| `export`       | `image` (size, crop, transparent, format), `image-copy`, `preview` (ui-only); extensions add `models`, `geometry`, `animation` (a session is `session.save`) | —           | export-image; extensions |
| `sequence`     | `style` (text, tracks), `colors`, `follow` (structure), `legend`                                                            | `get`                         | core                |
| `state`        | `remove`, `visibility`, `update` (a cell's params), `action` (an applicable state action by id)                             | `get` (the tree)              | state tree          |
| `settings`     | `set` (a group's values)                                                                                                    | `get`                         | core                |
| `history`      | `undo`, `undo-turn` (redo later, sequence.md §7.3)                                                                          | `list`                        | core                |
| `view`         | `open` (tab, card, section, page; `reveal` a node), `expand`, `fullscreen`, `theme` (UI only, no history step)               | —                             | core                |

- **Hiding contents** (decision 18): `contents.visibility` makes each component of the structure hold "what it held,
  except the hidden entities and chains" (a MolScript `exceptBy` over the component's own selection); the components'
  own params and the hidden set are the viewer's state (`registerCustomState('prism-viewer.contents')`, so snapshots
  and checkpoints carry them), and the operations that make components (`structure.preset`, `component.add`) apply it
  inside their step. A component left empty is a null cell, restored when the part is shown again.
- **Events** for embedding hosts (v5's `structureInteractivity`): `viewer.events.hover` and `viewer.events.click`,
  with the target under the pointer (README §3.6).
- **From** says which feature brings the area (viewer-spec.md §3.1): `core` is always there; a feature's areas come and
  go with it, and so do their commands and tools.
- **Registries are arguments, not operations**: one `structure.preset` takes a preset id, one `representation.color`
  takes a theme name, one `select.set` takes a selection query's name as its target; the registry lists them (as
  completions, as an enum in the tool's schema, refreshed by the registry-state behavior, state.md §3). Commands may
  still show one entry per preset (§3).
- **Extensions add areas** under their id (`geo-export.export`, `assembly-symmetry.enable`); only the viewer's built-in
  features add operations to core areas (`contributes`: `canvas.lighting`, `load.entry`; viewer-spec.md §3).

## 3. Commands over the operations

- **One command per operation by default**, with the operation's id, title and params; its form is the params' form, its
  inline sentence the flat fields (commands.md §2.1).
- **Generated commands** where a list reads better in the menu: one `preset.<id>` per preset, `styles.<id>` per quick
  style, `select.query.<name>` per selection query, one per applicable state action; each runs the same operation with
  that argument fixed.
- **Opening commands** (`measure.open`, `superpose.open`, `snapshots.show`, `panel.<tab>`) run `view.open`: they are
  the way into a custom UI, whose own actions are the area's operations.
- **Keybindings and toolbar buttons name commands**, never operations directly, so hiding a command hides them too.
- The menu's groups and prefixes (commands.md §5, §10) stay; the `select` group's shared `mode` is `select.set`'s
  `mode`.

## 4. Coverage

A test enumerates the operations and fails for one with neither a command nor a tool, unless it says why
(`surface: 'ui-only' | 'agent-only'` with a reason). The parity matrix's API column names the areas; each area's file
lists its operations in its header comment, citing this section.
