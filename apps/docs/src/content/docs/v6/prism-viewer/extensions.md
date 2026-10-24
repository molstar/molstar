---
title: Writing a Prism viewer extension
description: "The ViewerExtension fields, part ids and slots, lazy extensions, overrides, sequence tracks and behaviors, and the model quality extension as a worked example."
---

An extension adds to the viewer: engine pieces (color themes, behaviors), operations with their commands, annotations,
tracks in the sequence, skills for the assistant, and views. It is a plain object that you put in the spec's
`extensions` list. The extensions are code in the Mol\* repository. The types are the viewer's own, and there is no
stability promise beyond the repository's.

<!-- prettier-ignore-start -->
:::caution[Work in progress]
`apps/prism-viewer/src/extensions.ts` lists the extensions the app ships: model quality, model and geometry export, video export, Zenodo, assembly symmetry, membrane orientation, Kinemage and the data-only ones.
:::
<!-- prettier-ignore-end -->

An extension should start with operations, commands and skills, and add views only where a view is needed.

## ViewerExtension

```ts
import type { ViewerExtension } from '@molstar/prism-viewer/spec';
```

| Field         | What it is                                                                                                                                    |
| ------------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| `id`          | Lowercase words joined by dashes. Every part id starts with it. It cannot be the id of a core API area.                                       |
| `title`       | For people.                                                                                                                                   |
| `plugin`      | The engine side: `{ registry?, behaviors?, config? }`. Registry entries (color themes, representations, formats), behaviors and config items. |
| `api`         | Its API areas, by area id: `{ [areaId]: (ctx) => new MyArea(ctx) }`. Each is made once. The operation ids start with the area id.             |
| `commands`    | `(api) => CommandDefinition[]`. Called once its operations exist.                                                                             |
| `annotations` | `AnnotationProvider[]`: what it can show of a structure. They appear in `annotations.list` and `annotations.apply`.                           |
| `sequence`    | `SequenceBehavior[]`: tracks and menus in the sequence view. See below.                                                                       |
| `skills`      | `ViewerSkill[]`: tools for the assistant over its operations. A skill, never a core tool.                                                     |
| `ui`          | `(ctx) => ViewerExtensionUI`: its views. A function, so the views can be a chunk of their own.                                                |
| `settings`    | Reserved for its group in Settings. Not shown yet.                                                                                            |

`builtIn` and `contributes` are for the viewer's own features: only they may add operations to core areas such as `load`
and `canvas`. An extension adds an area of its own.

An extension is headless-safe when `plugin`, `api`, `commands`, `annotations`, `sequence` and `skills` import no DOM.
Keep the DOM in `ui`.

### Parts and ids

Anything with an id that the viewer registers is a part: a tab, a section, a toolbar item, an overlay, a sequence
behavior. A part id is `<owner>.<kind>.<name>`, where the owner is the extension id:

```text
model-archive.sequence.plddt
geo-export.section.export
snapshots.tab.snapshots
core.toolbar.reset              (the core's own parts have the owner core)
```

The kinds are `tab`, `section`, `row`, `selection-action`, `node-action`, `toolbar`, `overlay`, `command`, `tool`,
`skill`, `sequence`, `api` and `settings`. Operation ids and command ids are their own: `<area>.<verb>`. A duplicate id
is an error at creation, and the message names both owners.

### Views

`ui` returns any of these. Each id is a part id with the matching kind.

| Field              | What it is                                                                                                                                                                                                                                          |
| ------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sections`         | A section `{ id, slot, title, icon?, order?, keywords?, render }` in one of the viewer's slots.                                                                                                                                                     |
| `tabs`             | A tab `{ id, title, icon, dock?, order?, render }` for a large view of its own.                                                                                                                                                                     |
| `overlays`         | A card over the canvas `{ id, title, render }`, opened by a toolbar item or a command.                                                                                                                                                              |
| `toolbar`          | An item `{ id, icon, title, command?, overlay?, group?, order?, place?, visible? }` that runs a command or opens one of its overlays.                                                                                                               |
| `rows`             | Rows in the Structure tab's tree `{ id, when(node), whenChanged?, render }`, under the nodes `when` holds for (a structure's assembly symmetry).                                                                                                    |
| `selectionActions` | Rows of the selection bar's menu `{ id, title, hint?, icon?, order?, command?, run? }`: what to do with the selection, run as a command or by `run`.                                                                                                |
| `nodeActions`      | Rows of the scene tree's node menus, which the inspector shows too, `{ id, title, icon?, order?, kinds, when?(node), run(node) }`: on the node kinds in `kinds` that pass `when`, after the node's own rows and before Remove. `run` gets the node. |

The slots are a closed list that the viewer owns. An extension cannot add rail buttons or cards of its own beyond these.

| Slot               | Where a section goes                                                                  |
| ------------------ | ------------------------------------------------------------------------------------- |
| `home.import`      | Home, under the loading actions.                                                      |
| `selection`        | The Selection tab.                                                                    |
| `export`           | The export card.                                                                      |
| `canvas-settings`  | The canvas settings card.                                                             |
| `animation`        | The animation card.                                                                   |
| `settings`         | Settings, after the viewer's groups.                                                  |
| `inspector:<kind>` | The inspector of a node kind, such as `inspector:structure` or `inspector:component`. |

`render` is a Prism view, `domView(...)` from `@molstar/prism-ui`. The context `ui` receives has `api` (the API as the
person), `sources` (read-only states over the plugin) and `managers`. It has no plugin, so a view changes the scene
through operations.

### Overrides

The composition changes any part by its id, whoever registered it:

```ts
const spec = {
  ...FullPrismViewerSpec(),
  overrides: {
    'model-archive.sequence.plddt': null,                        // no pLDDT track
    'geo-export.section.export': { title: 'Export geometry' },    // a new title
    'model-archive.color-qmean': null,                           // a command, by its command id
  },
};
```

An id that nothing registers is an error at creation. To drop a whole extension, leave it out of the list.

## An operation, a command

An area is a class with its operations as fields. The area's `id` is its key in `api`. This one frames the ligands of a
structure by composing two core operations:

```ts
import { ApiArea, ApiResult } from '@molstar/prism-viewer/api';
import { operationCommand } from '@molstar/prism-viewer/commands/operations';
import type { ViewerExtension } from '@molstar/prism-viewer/spec';

class LigandViewApi extends ApiArea {
  readonly id = 'ligand-view';

  readonly show = this.operation({
    id: 'ligand-view.show',
    title: 'Show the ligands',
    params: {},
    history: false,                                  // a view change: not an undo step
    async run(_, ctx) {
      const selected = await ctx.api.select.set({ target: 'ligand' });
      if (!selected.ok) return selected;
      const framed = await ctx.api.camera.frame({ target: 'selection' });
      if (!framed.ok) return framed;
      return ApiResult.ok('Showing the ligands.', { selection: true, view: true });
    },
  });
}

declare module '@molstar/prism-viewer/api' {
  interface ExtensionApis { 'ligand-view': LigandViewApi }     // types `viewer.api['ligand-view']` for scripts that import this
}

export const LigandViewExtension: ViewerExtension = {
  id: 'ligand-view',
  title: 'Ligand view',
  api: { 'ligand-view': (ctx) => new LigandViewApi(ctx) },
  commands: (api) => {
    const area = api['ligand-view'];
    return area ? [operationCommand(api, area.show, { category: 'Structure', icon: 'eye' })] : [];
  },
};
```

The `params` of an operation are one parameter definition (`ParamDefinition` from
`@molstar/core/util/param-definition`). The command's form, the assistant's tool fields and the validation all derive
from it. Besides `history: false`, the definition can say:

| Field                 | Meaning                                                                                                                                                            |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `read: true`          | It changes nothing: no history step.                                                                                                                               |
| `undoable: false`     | The history cannot restore it, for example an upload: no checkpoint. The assistant asks every time.                                                                |
| `approval`            | `'network'`, `'upload'` or `'destructive'`: the assistant asks before running it. See [Assistant](/docs/v6/prism-viewer/assistant/#what-it-may-do-without-asking). |
| `addresses(values)`   | For `network`: the addresses a call reaches, so the assistant can check them.                                                                                      |
| `modules`             | What it needs of the host: `canvas`, `files`, `network`, `storage`, `page`. Without it, it fails with `unavailable`.                                               |
| `when(ctx)`           | `true` or the reason it cannot run now.                                                                                                                            |
| `exists(values, ctx)` | For operations that create something: a message when it exists already, and the run does nothing.                                                                  |

The `declare module` block types the area for scripts that import the extension. An operation's own calls go through
`ctx.api`, which is part of its step: if `show` were an undoable operation, the select and the camera calls would be
inside its one undo step.

## Lazy extensions

A heavy extension can load on first use. `LazyViewerExtension` has the same `id` and `title`, a `load` function, and a
`stub` that says what it will bring:

```ts
import type { LazyViewerExtension } from '@molstar/prism-viewer/spec';

export const LazyLigandView: LazyViewerExtension = {
  id: 'ligand-view',
  title: 'Ligand view',
  stub: {
    parts: [],
    operations: [{ id: 'ligand-view.show', title: 'Show the ligands' }],
    commands: [{ id: 'ligand-view.show', title: 'Show the ligands', category: 'Structure' }],
  },
  load: async () => (await import('./ligand-view.js')).LigandViewExtension,
};
```

The stub is what exists before the load: operations, command titles, placeholders for views, and skills for the
assistant's catalog. Calling a stub operation, running a stub command or disclosing its skill loads the chunk as a task
in the status column, swaps the real extension in, and replays what asked for it.

- A stub must declare every part, operation and command the extension brings. A loaded extension that brings what the
  stub did not declare, or lacks an operation it declared, or declares different flags (`read`, `undoable`, `approval`,
  `modules`), fails the load. Activation is all or nothing: the stub is put back and a later call tries again.
- `stub.owns` lists transformer ids, formats, file extensions, representations and themes that the extension answers
  for. Opening a file with one of its extensions, or a `.molj` session that needs one of its transformers, loads the
  extension first.
- `preload: true` loads it during creation, before `ready`, still as a chunk of its own.
- Overrides apply to the stub's parts too, so they are checked when the spec is created.

`@molstar/prism-viewer/features/lazy` shows the pattern: the heavier built-in features as stubs whose `load` imports the
feature's module. Put the stub in a module that does not import the rest of the extension, or the chunk is loaded with
the stub.

## Sequence tracks

The sequence band shows one structure at a time; the structure tree has a Sequence node under every structure, whose
children are its tracks. A **track provider** is how an extension adds sequence data: it says what its track is and when
it applies, and the viewer puts it in and takes it out, per structure. The person adds it (Add track), removes and hides
it; the agent does the same through `sequence.track.add`, `.remove`, `.hide` and `.show`. Each change is a step.

```ts
import type { SequenceTrackProvider } from '@molstar/prism-viewer/spec';
```

| Member                 | What it is                                                                                                                                                    |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `id`, `label`, `order` | The id is `<extension id>.sequence.<name>`; the label names it in the tree and under Add track; `order` places its track (ruler 0, letters 10, secondary 20). |
| `track`                | Its track's ref in the document (start it with the extension's id).                                                                                           |
| `applies(structure)`   | The structure carries its data. Read only the structure and what the viewer holds: a provider fetches nothing.                                                |
| `shown`                | In by default; `false` for a track only offered under Add track.                                                                                              |
| `build(ctx)`           | Its data nodes and its one track (Seq\* document nodes) over `ctx.structure`, `ctx.complex` (the complex ref a computed node names) and `ctx.letters`.        |
| `menu`, `legend`       | The track's menu (colorings of the 3D view by theme; the first is its icon in the band's gutter) and its legend blocks, shown while it is in.                 |
| `whenChanged`          | A state told when `applies` may answer otherwise without a new structure (a property computed later).                                                         |

List providers in the extension's `sequence.tracks`. The viewer's own (ruler, letters, secondary structure) are
providers too, in `spec.sequence.tracks`, and `overrides` replaces or removes any of them by id.

## Sequence behaviors

For what is not a track, a rule that follows the scene, a **sequence behavior** adds to the document by patch. It reads
what the sequence shows and proposes changes. List them in the extension's `sequence.behaviors`.

```ts
import type { SequenceBehavior } from '@molstar/prism-viewer/spec';
```

| Member                  | What it is                                                                                                                                                                                                                                                                        |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `id`, `title`           | The id is `<extension id>.sequence.<name>`; the host checks it.                                                                                                                                                                                                                   |
| `attach(ctx)`           | Starts the behavior and returns its dispose function.                                                                                                                                                                                                                             |
| `ctx.shown`             | A `ReadonlyState` of what the document shows. A new value for each document (another structure, another letter theme), `undefined` for nothing. It carries the `structure`, the `complex` ref, the `data` node, the `section` and the track to put yours `after` (none: the end). |
| `ctx.patch(shown, ops)` | Applies Seq\* operations to that document, all or nothing. It returns false and does nothing if the document was replaced meanwhile.                                                                                                                                              |
| `ctx.offer(offer)`      | Gives the track menus (`menus`, by track ref) and the legend blocks (`legend`). `undefined` withdraws them.                                                                                                                                                                       |

A new document has none of the behavior's nodes, so the behavior adds them again whenever `shown` changes. Give its
nodes refs that start with its id, and make `patch` insert them.

## An example: model quality

`packages/extensions/model-archive/src/prism-viewer.ts` is the `model-archive` extension, shown in the app. It makes the
quality scores of predicted models visible: pLDDT and QMEAN.

```ts
export const ModelArchiveExtension: ViewerExtension = {
  id: 'model-archive',
  title: 'Model quality',
  plugin: {
    registry: [{ structure: { themes: { color: [PLDDTConfidenceColorThemeProvider, QmeanScoreColorThemeProvider] } } }],
    behaviors: [],
  },
  annotations: QualityAnnotations,
  sequence: { tracks: [PlddtTrackProvider] },
  commands: (api) => [
    operationCommand(api, api.annotations.apply, {
      id: 'model-archive.color-plddt',
      title: 'Color by pLDDT',
      category: 'Structure',
      icon: 'paintbrush',
      alias: ['plddt', 'confidence', 'alphafold'],
      values: { annotation: 'plddt' },
    }),
    // ... the same for QMEAN
  ],
};
```

It uses four of the fields and no views:

1. **`plugin.registry`** registers the two color themes with the engine, so
   `representation.color({ theme: 'plddt-confidence' })` works for every caller.
2. **`annotations`** are two providers, `plddt` and `qmean`, each with `isApplicable(structure)` and an `apply` that
   calls `ctx.api.representation.color(...)`. A provider that applies through other operations is part of the step of
   `annotations.apply`. `api.annotations.list()` tells which a structure carries, and
   `api.annotations.apply({ annotation: 'plddt' })` colors it.
3. **`commands`** make "Color by pLDDT" and "Color by QMEAN" from `annotations.apply` with the `annotation` value fixed.
   `operationCommand` gives the command its form, its availability and its undo entry.
4. **`sequence`** has a track provider: a pLDDT track under the viewer's tracks, in by default for a structure that
   carries pLDDT.

The provider, reduced:

```ts
export const PlddtTrackProvider: SequenceTrackProvider = {
  id: 'model-archive.sequence.plddt',
  label: 'pLDDT',
  order: 30,
  track: PlddtTrack,
  applies: (structure) => QualityAssessment.isApplicable(structure.model, 'pLDDT'),
  build: (ctx) => ({
    data: [{ kind: 'computed', ref: PlddtData, params: { provider: 'molstar:plddt', options: { of: ctx.complex } } }],
    track: { kind: 'track', ref: PlddtTrack, params: { label: '', height: 2 },
      children: [ /* a layer drawing the values as a line in the four bands */ ] },
  }),
  menu: [{ id: 'color-plddt', label: 'Color by pLDDT', icon: 'paintbrush', theme: PlddtTheme }],
  legend: [{ title: 'Model confidence', items: Bands.map((b) => ({ label: b.label, color: b.color })) }],
};
```

The track's menu item colors the 3D view with the theme. The legend blocks appear in the sequence legend card. The
viewer itself knows nothing of pLDDT: the whole feature is in this one package, and the app turns it on by listing it in
`MolstarViewerExtensions`.

The extension has no skill and no view yet, and the pairwise (PAE) plot waits for a view of its own. The other ported
extensions (`apps/prism-viewer/src/extensions.ts`) follow the same shape: model and geometry export with Export
sections, video export and Zenodo as lazy extensions with stubs (`prism-viewer-stub.ts`), assembly symmetry, ANVIL's
membrane and Kinemage as operations, commands and skills.
