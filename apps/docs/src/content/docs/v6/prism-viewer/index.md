---
title: Prism viewer
description: "The Mol* v6 viewer on the Prism UI: what it is, how to run it, its URL parameters and the main areas of its interface."
---

The Prism viewer is the v6 viewer for molecular structures: Mol\*'s own tools (loading, the scene, selections, the
sequence, snapshots, export) on the new Prism UI. It does what the v5 Viewer does, with each panel designed for its task
rather than ported.

Everything the viewer changes goes through one API. A person using a form or the command menu, a script on the page, the
URL and the optional [assistant](/docs/v6/prism-viewer/assistant/) all call the same operations, and every operation
returns the same kind of result. See [the API](/docs/v6/prism-viewer/api/).

<!-- prettier-ignore-start -->
:::caution[Work in progress]
The viewer is not published to npm yet. The package names below are the repository's workspace packages, and the API can still change. Volumes and segmentations are not ported yet.
:::
<!-- prettier-ignore-end -->

## What there is

| Piece                   | What it is                                                                                                                       |
| ----------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| `apps/prism-viewer`     | The app: a page, the URL parameters, Mol\*'s extensions, the assistant. Built as ES modules with chunks.                         |
| `@molstar/prism-viewer` | The library the app is made from: `createPrismViewer`, the API, the extension types, a headless entry.                           |
| `prism-viewer.js`       | The app's single-file build: one ES module a page imports with no build step. See [Embedding](/docs/v6/prism-viewer/embedding/). |

## Run the app

From the repository root:

```sh
pnpm install
pnpm dev:apps -- prism-viewer   # watch mode, prints the page URL
pnpm build:apps                 # builds every app, this one included, into apps/prism-viewer/build
```

The build holds `index.html`, `index.js` with its `chunks/` and `prism-viewer.js`.

## Page globals

Once the viewer is made, the app sets `window.molstar` for the browser console and for scripts on the page:

| Name        | What                                                                                                                                                                                                                                                                     |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `viewer`    | The viewer (`viewer.plugin` is its Mol\* plugin).                                                                                                                                                                                                                        |
| `api`       | `viewer.api`. See [API](/docs/v6/prism-viewer/api/).                                                                                                                                                                                                                     |
| `lib`       | Mol\*'s library without a build step, as v5's Viewer has it: `structure`, `volume`, `shape`, `loci`, `math`, `plugin` (`PluginContext`, `PluginSpec`, `DefaultPluginSpec`, `StateActions`, …) and `extensions` (MolViewSpec, wwPDB struct-conn, Model Archive, loaders). |
| `lib.prism` | The Prism viewer's own: `createPrismViewer`, `DefaultPrismViewerSpec`, `FullPrismViewerSpec`, `node` (writes a [node selector](/docs/v6/prism-viewer/api/#nodes)) and `MS` (the MolScript builder).                                                                      |

```js
const { api, lib } = window.molstar;
const ligand = lib.prism.node('component', { type: 'ligand' }).child('representation');
await api.representation.visibility({ node: ligand, visible: false });
const all = lib.prism.MS.struct.generator.all();
```

`lib` leaves out v5's classic UI (`PluginUIContext`, `DefaultPluginUISpec`, the view models) and `StateTransforms`:
import a transformer from its module. The app loads it beside the viewer, in a chunk of its own; `prism-viewer.js` and
`@molstar/prism-viewer` do not include it.

## URL parameters

The app reads its parameters from the query string, with v5's names. Names are case-sensitive (lower case), the first of
a name wins, and an empty value counts as absent. A parameter the app does not know is ignored, and the console logs the
ignored names once.

### Loading

The loading parameters run once the viewer is ready, as operations of the
[`load` area](/docs/v6/prism-viewer/api/#areas) with the `url` origin, in the order of this table. A load that fails is
reported and does not stop the ones after it.

| Parameter                                                                 | Loads                                                                          |
| ------------------------------------------------------------------------- | ------------------------------------------------------------------------------ |
| `snapshot-id`                                                             | A state from the state server (`load.snapshot`, which is `session.open`).      |
| `snapshot-url`, `snapshot-url-type`                                       | A `.molx` or `.molj` state. Without the type, the extension, else the content. |
| `structure-url`, `structure-url-format`, `structure-url-is-binary`        | A structure file (`load.url`). The format defaults to `mmcif`.                 |
| `mvs-url` or `mvs-data`, `mvs-format`                                     | A MolViewSpec file or inline data (`load.mvs`). Without a format, it is told.  |
| `pdb`, `pdb-ihm` (`pdb-dev` is the same), `emdb`, `afdb`, `model-archive` | A database entry by id (`load.entry`).                                         |
| `url`, `url-format`, `url-is-binary`                                      | Any file by address (`load.url`). Without a format, the address decides.       |

The `*-is-binary` parameters are true for `1` or `true`; left out, the format and the address decide (a `.bcif` address
is binary).

### Options

These are read when the viewer is created.

| Parameter                                                              | Sets                                                                                           |
| ---------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------- |
| `pixel-scale`, `pick-scale`                                            | The resolution of the canvas and of picking.                                                   |
| `transparency`                                                         | `auto`, `blended`, `wboit` or `dpoit`.                                                         |
| `resolution-mode`                                                      | `auto`, `scaled` or `native`.                                                                  |
| `power-preference`                                                     | `default`, `high-performance` or `low-power`.                                                  |
| `graphics-backend`                                                     | `auto` (default: WebGPU where available, else WebGL), `webgpu`, `webgl`, `webgl2` or `webgl1`. |
| `allow-major-performance-caveat`, `disable-antialiasing`               | Flags.                                                                                         |
| `illumination`                                                         | Flag: path-traced global illumination.                                                         |
| `bg-color` (or `viewport-background-color`)                            | The canvas background, a CSS color.                                                            |
| `viewport-focus-behavior`                                              | `default`, `secondary-zoom` or `disabled`.                                                     |
| `pdb-provider`                                                         | `pdbe`, `rcsb` or `pdbj`.                                                                      |
| `emdb-provider`                                                        | `pdbe` or `rcsb`.                                                                              |
| `volume-streaming-server`, `volume-streaming-disabled`, `state-server` | The servers the viewer talks to.                                                               |
| `collapse-left-panel` (or `hide-controls`)                             | Flag: start with no panel open.                                                                |

A flag is true for an empty value, `1` or `true`.

A few v5 parameters are accepted without a warning but do nothing yet: `map-provider`, `show-toggle-fullscreen`,
`debug-mode`, `timing-mode` and `pick-padding`.

Examples:

```text
?pdb=1tqn
?afdb=Q8W3K0&bg-color=white
?mvs-url=https://example.org/story.mvsj
?pdb=1cbs&illumination=1&collapse-left-panel
```

The URL is read by the app, not by the library. A page that creates its own viewer passes the same settings as
[options](/docs/v6/prism-viewer/embedding/#options) and loads with the API.

## The interface

- **Home** loads things: the registry's actions, drag and drop of files, and the database and URL forms.
- **Structure** shows the scene as a tree of structures, components and representations, with presets, quick styles and
  an inspector for the node picked.
- **Selection** holds what you can do with the current selection. Measurements and superposition add sections to it.
- **Sequence** is a band over the canvas with the sequence of the structure, in tracks or as text, linked to the 3D view
  in both directions: hover, selection and focus.
- **The canvas** has the camera, reset, export, fullscreen, canvas settings and animation controls, a hover label, and a
  switch for the selection mode.
- **History** lists what happened in the session, newest first, and undoes the latest step.
- **Settings** holds the viewer's and the canvas's settings. **Help** lists commands and mouse bindings. **Log** keeps
  the messages.
- **The command menu** lists every command with its form and opens any of them. The status column shows the result of
  what you ran, running tasks and warnings.
- **Snapshots** and the **State tree** are further tabs, loaded the first time they are used.
- **Assistant** is an optional tab and an Ask row in the command menu. It exists only when the app includes the
  [assistant](/docs/v6/prism-viewer/assistant/).

Hiding or replacing any of these is done with [overrides](/docs/v6/prism-viewer/embedding/#the-spec).
