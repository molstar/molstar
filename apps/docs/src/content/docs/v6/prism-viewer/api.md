---
title: The Prism viewer API
description: "Operations by area, results, origins, undo and groups, targets as strings, describing operations, and the events a script can watch."
---

`viewer.api` is what scripts, embedding pages, the viewer's own commands and the assistant all call. It has one object
per area, and each area has operations as methods:

```ts
const viewer = await createPrismViewer({ target, spec });
const { api } = viewer;

await api.run('load.entry', { source: 'pdb', id: '1tqn' });
await api.styles.apply({ style: 'cartoon' });
await api.representation.color({ target: 'ligand', color: 'red' });
await api.select.set({ target: 'A:120-130' });
await api.camera.frame({ target: 'selection' });
```

Every operation takes one object of values and returns a promise of a result. Values you leave out take their defaults,
so `api.component.add({ target: 'resn HEM' })` is a whole call. Operations a feature adds to a core area (the loaders'
`load.entry`) are optional members (`api.load.entry?.(…)`); `api.run('load.entry', …)` reaches them by id. A second
optional argument is `{ signal?: AbortSignal; copy?: boolean }`.

Operations never throw. A failure is a result.

## Results

```ts
const r = await api.component.add({ target: 'resn HEM', label: 'Heme' });
if (r.ok) {
  console.log(r.did);        // one line, for a person
} else {
  console.warn(r.code, r.message, r.hint);
}
```

| Field       | When   | Meaning                                                                      |
| ----------- | ------ | ---------------------------------------------------------------------------- |
| `ok`        | always | `true` or `false`.                                                           |
| `did`       | `ok`   | What happened, in one line.                                                  |
| `changed`   | `ok`   | What a reader should read again: `{ nodes?, selection?, view?, sequence? }`. |
| `unchanged` | `ok`   | `true` when it was already so. Nothing changed and no history step was made. |
| `value`     | `ok`   | The data of a read, or what the operation made.                              |
| `warnings`  | `ok`   | Things worth knowing.                                                        |
| `next`      | `ok`   | An obvious next step, when there is one.                                     |
| `code`      | `!ok`  | Why it failed (below).                                                       |
| `message`   | `!ok`  | What went wrong, for a person.                                               |
| `hint`      | `!ok`  | What to do instead.                                                          |

An operation states a desired state: running `representation.color` with the same values twice changes nothing the
second time, and the second result is `unchanged`. An operation that creates something checks whether it exists:
`load.entry` for an entry that is loaded, or `component.add` for a component that is there, does nothing and says so.
Pass `{ copy: true }` as the second argument to get another.

Failure codes:

| Code          | Meaning                                                              |
| ------------- | -------------------------------------------------------------------- |
| `invalid`     | The values are wrong.                                                |
| `blocked`     | Not possible now, for example Undo with nothing to undo.             |
| `unavailable` | The host lacks what the operation needs, such as a canvas headless.  |
| `not-found`   | Something it names does not exist, or a target selected nothing.     |
| `ambiguous`   | More than one thing matches.                                         |
| `confirm`     | Ask the person first, then call again.                               |
| `aborted`     | The signal was aborted.                                              |
| `failed`      | It ran and failed. A failed or aborted run restores what it changed. |

## Areas

The operation `component.add` is `api.component.add(...)`, and it is the same operation the command menu, the run log
and the assistant call. A hyphen in an operation id is camel case on the area: `history.undo-turn` is
`api.history.undoTurn`, `superpose.tm-align` is `api.superpose.tmAlign`. These are the areas the viewer has today. The
core's (`load`, `session`, `select`, `focus`, `highlight`, `camera`, `canvas`, `settings`, `history`, `view`) are always
there. The built-in features' areas (`structure`, `contents`, `component`, `representation`, `styles`, `annotations`,
`trajectory`, `volume`, `particles`, `measure`, `snapshot`, `animation`, `export`) and their operations in core areas
(`load.entry`, `load.url`, `session.save`, `canvas.lighting`, …) come with the default spec and are typed as present, so
`await api.load.entry({ id: '1tqn' })` needs no `!`: reading one in a viewer whose spec left its feature out throws,
naming the feature. The lazy ones (`superpose`, `state`) and remote states are optional in the type, so check for them
(`api.superpose?.chains(...)`) when you do not control the spec; `api.area(id)` and `api.describe(id)` ask for any.

| Area                  | Operations                                                                                                                                                       |
| --------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `load`                | `files`, `action`; with the `loaders` feature: `entry`, `url`, `data`, `mvs`, `svs`, `snapshot`, `trajectory`, `volume`, `viewport`.                             |
| `session`             | `open`, `save`, `clear`; with remote states: `remote.list`, `remote.open`, `remote.upload`, `remote.remove`.                                                     |
| `structure`           | `preset`, `model`, `assembly`, `symmetry`, `visibility`, `remove`, `get`, `transform`, `transform-reset`, `unit-cell`.                                           |
| `contents`            | `get`, `visibility`: hide an entity or a chain from every component of a structure.                                                                              |
| `component`           | `add`, `modify`, `remove`, `rename`, `visibility`, `options`, `list`.                                                                                            |
| `representation`      | `set`, `update`, `visibility`, `remove`, `color`, `size`, `transparency`, `emissive`, `material`, `clipping`, `reset`, `layers`, `layer-update`, `layer-remove`. |
| `styles`              | `apply`: a quick style, `default`, `cartoon`, `spacefill`, `surface`, `plain` or `illustrative`.                                                                 |
| `annotations`         | `list`, `apply`: show what a structure carries (B-factors, secondary structure, and what extensions add).                                                        |
| `select`              | `set`, `all`, `clear`, `invert`, `around`, `expand`, `mode`, `granularity`, `get`.                                                                               |
| `focus`               | `set`, `clear`, `get`, `again`, `history`.                                                                                                                       |
| `highlight`           | `set`, `clear`, `mark`, `unmark`.                                                                                                                                |
| `camera`              | `reset`, `reset-axes`, `zoom`, `orient`, `frame`, `set`, `get`.                                                                                                  |
| `canvas`              | `background`, `palette`, `xr`, `get`; with `canvas-settings`: `settings`, `lighting` (its switches, illumination among them), `camera`, `clipping`.              |
| `trajectory`          | `frame`, `play`, `stop`, `get`.                                                                                                                                  |
| `volume`, `particles` | `get`, `representation`, `update`, `visibility`, `remove`, and more for volumes: `isovalue`, `streaming`, `streaming-view`.                                      |
| `sequence`            | `style`, `colors`, `follow`, `get`. Present when the spec has a sequence.                                                                                        |
| `settings`            | `get`, `set`.                                                                                                                                                    |
| `history`             | `undo`, `redo`, `undo-turn`, `list`.                                                                                                                             |
| `view`                | `open`, `close`, `expand`, `fullscreen`: open a tab, a card or a section, reveal a node, or expand or fullscreen the viewer. No history step.                    |
| `measure`             | `distance`, `angle`, `dihedral`, `label`, `orientation`, `plane`, `get`, `remove`, `visibility`, `options`.                                                      |
| `superpose`           | `chains`, `tm-align`, `ligands`, `atoms`, `db`, `list`, `remove`.                                                                                                |
| `snapshot`            | `add`, `apply`, `update`, `replace`, `move`, `remove`, `clear`, `next`, `previous`, `play`, `stop`, `list`, `link`, `seek`, `transition`.                        |
| `animation`           | `play`, `stop`, `wiggle`, `get`.                                                                                                                                 |
| `export`              | `image`, `image-copy`, `preview`, `viewport`.                                                                                                                    |
| `state`               | `get`, `select`, `visibility`, `update`, `apply`, `remove`, `action`: the raw state tree, by node selector. Lazy.                                                |

`packages/prism/viewer/_spec/operations.md` lists every operation with its parameters and flags, written from the code.

Extensions add areas of their own, under their id. A script that imports an extension gets its area typed. See
[Extensions](/docs/v6/prism-viewer/extensions/).

### Examples

```ts
// Load, then show a ligand.
await api.run('load.entry', { source: 'afdb', id: 'Q8W3K0' });
await api.run('load.url', { url: 'https://example.org/model.cif', format: 'mmcif' });
await api.run('load.mvs', { url: 'https://example.org/story.mvsj' });

// Select and focus. Neither is a history step.
await api.select.set({ target: 'A:248-256' });
await api.select.set({ target: 'resn HOH', mode: 'remove' });   // mode: set, add, remove or intersect
await api.focus.set({ target: 'ligand' });

// Color: a CSS color, or a theme by name.
await api.representation.color({ target: 'polymer', theme: 'chain-id' });
await api.representation.color({ target: 'selection', color: '#ff8800' });

// Make a component and give it a representation.
await api.component.add({ target: 'A:120-130', label: 'Loop', representation: 'ball-and-stick' });

// A quick style, a preset by registry id ('auto' is the default), the canvas.
await api.styles.apply({ style: 'illustrative' });
await api.structure.preset({ preset: 'auto' });
await api.canvas.background({ color: 'white' });
```

The values of each operation are described by `api.describe(id)` (below). Where a value takes a registry entry, such as
a preset id or a color theme, the failure message and hint list the choices.

### Reading

A read is an operation with a `value` and no history step. `api.select.get()`, `api.component.list()` and
`api.structure.get()` return text in `did` and data in `value`. A few areas also expose a `ReadonlyState` to subscribe
to: `structure.list`, `canvas.state`, `history.steps`, `annotations.providers` and `animation.state`. For the selection
and the focus as states, use `viewer.sources.interaction`.

```ts
const sel = await api.select.get();
if (sel.ok) console.log(sel.did, sel.value);

api.history.steps.subscribe((steps) => console.log(steps.length, 'undoable steps'));
```

## Targets

Operations that act on part of a structure take a `target`. A target is a string, and the same string means the same
thing in the API, the commands and the assistant's tools.

| Form             | Examples                                                                               |
| ---------------- | -------------------------------------------------------------------------------------- |
| A keyword        | `all`, `polymer`, `protein`, `nucleic`, `branched`, `ligand`, `ion`, `water`, `coarse` |
| The person's own | `selection` and `focus`: what is selected and focused now.                             |
| Addresses        | `A`, `B:120`, `B:200-210`, `B:120A`, `C:122@OG1`, `A:10-20, C`. Author numbering.      |
| Query text       | `resn HEM`, `resn HIS and chain A`, `within 5 of resn STI`: PyMOL by default.          |

A string is read as a keyword, then as addresses, then as query text. Reading is lenient about spaces, case and a
trailing comma. What cannot be read fails with the part that did not read and a hint.

Set `language: 'vmd' | 'jmol' | 'molscript'` to read the query text in another language. A target can also be schema
items, for example `{ auth_asym_id: 'A', beg_auth_seq_id: 248, end_auth_seq_id: 256 }`, or `{ molql: expression }` with
a MolScript expression as JSON. Whatever the form, everything compiles to MolQL, and that is what the history and the
session record store.

When several structures are loaded, say which with `structure: '1TQN'`, a label or a state ref. Without it a target
applies to the structure you are working with: the focused one, else the selected one, else the one in view, else every
loaded structure. The result says which it was.

## Nodes

Operations on cells of the state tree (`state.*`; a representation's `update`, `visibility` and `remove`; `volume.*` and
`particles.*`) take a `node`: a selector in the shape of CSS, compiled to Mol\*'s `StateSelection`.

```ts
await api.representation.visibility({ node: 'component[type=ligand] > representation', visible: false });
await api.state.select({ node: 'structure[label="1CBS"] representation:hidden' });
await api.state.apply({ node: 'component[type=ligand]', transformer: 'structure-representation-3d' });
```

| Part        | Examples                                                                                        |
| ----------- | ----------------------------------------------------------------------------------------------- |
| Kinds       | `data`, `trajectory`, `model`, `structure`, `component`, `representation`, `volume`, `*`        |
| Attributes  | `[label="1CBS"]`, `[type=cartoon]`, `[tag=mine]` or `.mine`, `[transformer$=component]`, `#ref` |
| Operators   | `=`, `*=` (contains), `^=` (starts with), `$=` (ends with); text compares without case          |
| States      | `:hidden`, `:visible`, `:error`, `:not(...)`                                                    |
| Combinators | a space (anything below), `>` (a child), `,` (either)                                           |

A ref alone works too. An operation acts on every match of its kind as one step; one that needs a single cell
(`state.get`, `state.action`) lists the matches when there are several. In scripts, `node()` from
`@molstar/prism-viewer/api` writes the same text: `node('component', { type: 'ligand' }).child('representation')`.

## Origins

Every call has an origin that says who made it. `viewer.api` runs as `script`. `api.as(origin)` gives the same API for
another:

```ts
const asPerson = api.as({ kind: 'person' });
await asPerson.styles.apply({ style: 'cartoon' });
```

| Kind                       | Used by                                  |
| -------------------------- | ---------------------------------------- |
| `person`                   | The commands and forms in the UI.        |
| `script`                   | `viewer.api`, and `viewer.run`.          |
| `agent` (with a `turn` id) | The assistant's tools.                   |
| `url`                      | The URL's loading parameters.            |
| `link`                     | The sequence view acting on the 3D view. |

History and the session record keep the origin of each step. The status column reports the results of a person's and the
URL's calls, not a script's.

## History and undo

The viewer has one undo history. Every undoable operation run is one step, and a checkpoint of the scene is taken before
it. Undo restores the latest checkpoint.

- Operations that change what you see and not the scene are not steps: `select`, `focus`, `highlight`, `camera`, `view`,
  trajectory frames, and reads. Undo does not restore the camera, the selection, the focus or the open panels.
- A result that is `unchanged` makes no step. An operation that cannot be restored, such as an upload, makes no step
  either and is only recorded.
- A run that fails or is aborted restores its checkpoint, so a load that fails half way leaves nothing behind.
- Undo is not available while a run is in flight. Redo puts back what Undo took back, until a new step is made.
- Changing the plugin outside the API (a script that calls a loader on `viewer.plugin` without `run`) still shows, but
  makes the steps before it unsafe: History marks them, and they can no longer be undone.
- For a very large scene, the viewer turns per-step undo off and says so once.

```ts
await api.history.undo();            // undoes the latest step; { ok: false, code: 'blocked' } when it cannot
await api.history.redo();            // puts back what Undo took back
const list = await api.history.list();   // list.value: the undoable steps, newest first
```

### Groups

`api.group` makes the calls inside it one step, and one group in the record:

```ts
const g = await api.group('Prepare the figure', async (api) => {
  await api.canvas.background({ color: 'white' });
  await api.styles.apply({ style: 'illustrative' });
});
// g.ok, g.did === 'Prepare the figure', g.runs: the calls, g.changed: all they changed
```

The group fails as a whole if the work throws, returns a failed result, or a call inside fails: the checkpoint is
restored and the result says so. `{ onFailure: 'continue' }` undoes only the failed call and goes on. Calls an operation
makes through its own `ctx.api` are part of its step in the same way.

For work that no operation covers yet, `viewer.run(label, (plugin, api) => ...)` is one step: it is recorded as one
opaque run, and the calls it makes through `api` join it.

### An assistant's turn

An assistant's turn is one step, a group with its `agent` origin. `api.history.undoTurn({ turn })` undoes it. If the
person made steps after the turn, it fails with the code `confirm`; call again with `discardLater: true` to undo those
too.

## Describing operations

```ts
const d = api.describe('representation.color');
if (d) {
  d.params;      // the values, as parameter definitions (a color theme select lists the current themes)
  d.available;   // true, or the reason it cannot run now
  d.supported;   // true, or the reason this viewer cannot run it at all
  d.definition;  // id, title, and the flags: read, history, undoable, approval, modules
}

api.definitions.subscribe((all) => console.log([...all.keys()]));   // changes as extensions come and go
```

`api.run(id, values?, options?)` runs an operation by its id, which is how commands, tools and the run log name them. An
id that does not exist gives a `not-found` result.

## Events and state

```ts
api.events.started.subscribe((run) => {});
api.events.finished.subscribe((run) => console.log(run.operation, run.origin.kind, run.result));
api.running.subscribe((runs) => {});                    // the runs in progress
```

A run record has `id`, `operation`, `title`, `values`, `origin`, `started`, `result`, `ms`, `step`, `group` and
`nested`.

The viewer's other observables are `viewer.state.ready`, `.busy` and `.empty`, `viewer.events.hover`, `.click` and
`.error`, and `viewer.managers.history` (`steps`, `canUndo`, `undoBlocked`, `running`, `checkpointing`).

States are `ReadonlyState`: an rxjs `Observable` that replays its current value on subscribe and also has a synchronous
`.value`. Events are plain observables that replay nothing.

## What is stable

The core areas and the built-in features' are the public contract: their ids, values and results follow semver once the
package is released. Areas that extensions add, and `viewer.plugin`, carry no such promise.
