---
title: Embedding the Prism viewer
description: "Create a Prism viewer in a page with createPrismViewer, use the single-file build, configure the spec and options, and run it headless in Node."
---

A page creates a viewer from a spec. The spec is plain data: lists of extensions, options, and what to leave out.

<!-- prettier-ignore-start -->
:::caution[Work in progress]
The packages are not published yet, so the import paths are the workspace's. The `<prism-viewer>` element is built, in the single-file build only; the viewport-only embedding that the design describes is not.
:::
<!-- prettier-ignore-end -->

## createPrismViewer

```ts
import { createPrismViewer, DefaultPrismViewerSpec } from '@molstar/prism-viewer';

const viewer = await createPrismViewer({
  target: document.getElementById('viewer')!,
  spec: DefaultPrismViewerSpec(),
});

await viewer.api.run('load.entry', { source: 'pdb', id: '1tqn' });
// ...
await viewer.dispose();
```

`target` is an element with a size. The viewer attaches an open shadow root to it and renders the shell and the canvas
inside, so the page's CSS does not reach in. To restyle the viewer, set the `--msp-*` custom properties on the target:
they cross the shadow boundary.

The promise resolves when the canvas is mounted and `viewer.state.ready` is true. A creation that fails undoes what it
made and rejects. An invalid spec rejects with a `PrismViewerSpecError` that names every problem at once.

`createPrismViewer` does not read the URL. See [URL parameters](/docs/v6/prism-viewer/#url-parameters) for what the app
does with it, and `runUrlLoads` below for the loading half.

### What a viewer gives you

| Member                                              | What it is                                                                                                                                                                                                    |
| --------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `viewer.api`                                        | Every operation, as a script calls them. See [the API](/docs/v6/prism-viewer/api/).                                                                                                                           |
| `viewer.run(label, (plugin, api) => ...)`           | Raw plugin work as one undoable step.                                                                                                                                                                         |
| `viewer.state.ready`, `.busy`, `.working`, `.empty` | `ReadonlyState`s: an rxjs `Observable` with a synchronous `value`. `busy` is the plugin's debounced flag; `working` is anything in progress now (an operation, a lazy extension's load, the plugin updating). |
| `viewer.events.hover`, `.click`, `.error`           | Observables. Hover and click carry the target under the pointer.                                                                                                                                              |
| `viewer.extensions`                                 | Which extensions are loaded, the status of each lazy one, `load(id)`, and `mount(ext)` to add one at runtime.                                                                                                 |
| `viewer.managers`                                   | `commands`, `history`, `navigation`, `viewState`, `notifications`, `transactions`.                                                                                                                            |
| `viewer.record`                                     | The session record that History shows.                                                                                                                                                                        |
| `viewer.plugin`                                     | The Mol\* plugin context. An escape hatch that the API does not promise to keep stable.                                                                                                                       |
| `viewer.dispose()`                                  | Tears the viewer down.                                                                                                                                                                                        |

### Load from the URL

The loaders feature has a function for the URL's loading parameters. The app calls it after creating the viewer:

```ts
import { runUrlLoads } from '@molstar/prism-viewer/features/loaders';

const results = await runUrlLoads(viewer.api, location.search);   // pdb=, afdb=, url=, mvs-url=, ...
```

Each load runs with the `url` origin. A failed load is in `results` and does not stop the others.

## The single-file build

The app builds a single ES module, `prism-viewer.js`, with every lazy part inlined and the CSS included. A page imports
it and needs no build step:

```html
<div id="viewer" style="position: absolute; inset: 0"></div>
<script type="module">
  import { createPrismViewer, FullPrismViewerSpec } from './prism-viewer.js';

  const viewer = await createPrismViewer({
    target: document.getElementById('viewer'),
    spec: FullPrismViewerSpec(),
  });
  await viewer.api.run('load.entry', { id: '1cbs' });
</script>
```

The module exports `createPrismViewer`, `FullPrismViewerSpec`, `DefaultPrismViewerSpec` and the element's
`createPrismViewerElementClass` and `definePrismViewerElement`. It has no assistant (see below), does not read the URL
and sets no globals: `molstar.lib`, Mol\*'s library for page scripts, is the app's
([Page globals](/docs/v6/prism-viewer/#page-globals)).

Importing it also defines the `<prism-viewer>` element, which wraps `createPrismViewer` with `FullPrismViewerSpec()`
until its `spec` property is set:

```html
<prism-viewer pdb="1cbs" theme="dark" style="height: 480px"></prism-viewer>
```

The attributes are `pdb` (an entry to load; a new value loads that one too), `theme` (`light`, `dark`, `blue`,
`contrast` or `auto`) and `hide-controls`; `theme` and `hide-controls` are read when the viewer is made.
`element.viewer` is the viewer once it is made, and `element.ready` a promise of it. Setting `element.spec` makes the
viewer again. Several on a page each own their viewer, and one leaving the page disposes its viewer.

`pnpm build:apps` writes it to `apps/prism-viewer/build/prism-viewer.js`. The `molstar` distribution carries it
(`build/prism-viewer/`, beside the app); it is not deployed yet.

## The spec

`DefaultPrismViewerSpec()` returns a new spec each call. It holds the viewer's base and its own features as extensions:

| Feature           | Brings                                                                                                               |
| ----------------- | -------------------------------------------------------------------------------------------------------------------- |
| `structure`       | The scene's API: the `structure`, `contents`, `component`, `representation` and `styles` areas.                      |
| `volume`          | The `volume` area.                                                                                                   |
| `particles`       | Lazy: Mol*'s particles (`ParticlesRegistry`) and the `particles` area, loaded by the first particle file or session. |
| `trajectory`      | The `trajectory` area and the trajectory bar.                                                                        |
| `annotations`     | The `annotations` area and Home's offer of a structure's annotations. Needs `structure`.                             |
| `home`            | Home.                                                                                                                |
| `structure-panel` | The Structure tab: the scene tree and its inspectors, Preset and Styles. Needs `structure`.                          |
| `selection-panel` | The Selection tab. Needs `structure`.                                                                                |
| `help`            | Help.                                                                                                                |
| `loaders`         | `load.entry`, `load.url`, `load.data`, `load.mvs`, `load.svs`, `load.snapshot`, `load.trajectory`, `load.volume`.    |
| `canvas-settings` | The canvas settings card and `canvas.lighting` (illumination too), `.camera`, `.clipping`.                           |
| `measurements`    | The `measure` area and the Measure and Label section.                                                                |
| `superposition`   | The `superpose` area and the Superpose section. Lazy.                                                                |
| `state-tree`      | The `state` area and the State tree tab. Lazy.                                                                       |
| `remote-states`   | Remote states in the `session` area. Lazy, and left out unless `options.ui.showRemoteStates` is on.                  |
| `snapshots`       | The `snapshot` area and the Snapshots tab.                                                                           |
| `animation`       | The `animation` area and the animation card.                                                                         |
| `export-image`    | The `export` area: images.                                                                                           |
| `export-session`  | `session.save`.                                                                                                      |

Leave a feature out and its code is not in your bundle: `examples/minimal-prism-viewer` lists only `structure` and is
under half the size of the full viewer. The panels keep their tab ids (`home`, `structure`, `selection`, `help`) when
they are there. Reading a built-in area or view model a spec left out (`viewer.api.volume`, `viewer.views.scene`) throws
an error naming the feature to add; `viewer.api.area('volume')` and `'scene' in viewer.views` ask without throwing.

`FullPrismViewerSpec()` is the default spec with Mol\*'s extensions added after the built-in ones. It lives in the app
(and in `prism-viewer.js`), because the extensions depend on the viewer: model quality
([the example](/docs/v6/prism-viewer/extensions/#an-example-model-quality)), model and geometry export, video export and
Zenodo (lazy), assembly symmetry, membrane orientation, Kinemage and the data-only extensions. The app's own spec is the
full spec with the URL's options, the assistant and `persist`.

A spec has these fields, all but `plugin` optional:

| Field        | What it is                                                                                                                                                                                                                                                                                                                    |
| ------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `plugin`     | The engine's base: registry entries, behaviors, config items and canvas props. Named parts in `@molstar/prism-viewer/sources/base`.                                                                                                                                                                                           |
| `options`    | Values for the options below, by group.                                                                                                                                                                                                                                                                                       |
| `extensions` | A list of extensions, in order. Leave one out to drop it. See [Extensions](/docs/v6/prism-viewer/extensions/).                                                                                                                                                                                                                |
| `languages`  | The query languages target text may be written in, PyMOL first if listed: `AllQueryLanguages` from `@molstar/prism-viewer/targets/languages/all`, or one each from `targets/languages/pymol`, `vmd`, `jmol`, `molscript`. None: keywords and address lists only, and no transpiler bundled. The default spec lists all four.  |
| `overrides`  | Replace or remove any part (a tab, a toolbar button, a command, a sequence behavior) by id.                                                                                                                                                                                                                                   |
| `sequence`   | `{ behaviors, style? }` for the one sequence viewer, `style` being `'tracks'` or `'text'`; `false` for none.                                                                                                                                                                                                                  |
| `agent`      | `{ core?, skills?, assistant? }`, or `false` for no tools at all.                                                                                                                                                                                                                                                             |
| `ui`         | `layout` (docks, tabs, the tab open at the start), `theme` and `keybindings`. Ignored headless.                                                                                                                                                                                                                               |
| `persist`    | `{ key, preferences?, restore?, rows? }`: a key names a viewer's stored data (two live viewers on a page cannot share one); `preferences` (`webStore('local', …)`) keeps the theme, `restore` (`webStore('session', …)`) this tab's light restore (`viewer.restore.run()` after a reload), `rows` where Keep writes a record. |
| `finalize`   | `(spec: PluginSpec) => PluginSpec`: the last word on the composed plugin spec.                                                                                                                                                                                                                                                |

Overrides name a part as `<owner>.<kind>.<name>` and a command by its own id. `null` removes, an object replaces fields:

```ts
const spec = {
  ...DefaultPrismViewerSpec(),
  overrides: {
    'core.toolbar.reset': null,                         // no reset button
    'camera.reset': null,                               // and no reset command
    'measurements.section.measure': { title: 'Distances' },
  },
};
```

An override of an id that nothing registers is an error at creation, naming it.

## Options

`spec.options` has four groups. Each is partial, and a missing value is the default.

| Group      | Values                                                                                                                                                                                                                         |
| ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `renderer` | `pixelScale`, `pickScale`, `transparency`, `resolutionMode`, `powerPreference`, `graphicsBackend`, `allowMajorPerformanceCaveat`, `disableAntialiasing`, `illumination`                                                        |
| `canvas`   | `background` (a CSS color; empty for the default), `focusBehavior`                                                                                                                                                             |
| `data`     | `pdbProvider`, `emdbProvider`, `volumeStreamingServer`, `volumeStreamingDisabled`, `stateServer`, `saccharideCompIdMapType`                                                                                                    |
| `ui`       | `showSequence`, `showHistory`, `showRemoteStates`, `collapsed`, `splash` (the empty canvas's splash; on), and `buttons`: `reset`, `screenshot`, `expand`, `fullscreen`, `settings`, `selectionMode`, `animation`, `trajectory` |

```ts
const viewer = await createPrismViewer({
  target,
  spec: {
    ...FullPrismViewerSpec(),
    options: {
      canvas: { background: '#ffffff' },
      ui: { collapsed: true, showSequence: false },
    },
    ui: { theme: 'dark' },
  },
});
```

`ui.theme` is `'light'`, `'dark'`, `'blue'`, `'contrast'` or `'auto'`.

## Several viewers on a page

Give each viewer its own `persist.key`, for example `persist: { key: 'left' }`. Keys go to the viewer that has focus,
not to the window.

A viewer without a `persist.key` keeps its records under the key `prism-default`, apart from the Prism Viewer app's
(`prism-viewer`). Viewers without a key on the same origin share one list of kept records: give a viewer its own key to
keep its records apart.

## The assistant

The assistant is not part of any default spec. The app adds it with `assistant()`, which returns a lazy loader:

```ts
import { assistant } from '@molstar/prism-assistant/spec';

const spec = {
  ...FullPrismViewerSpec(),
  agent: { assistant: assistant() },
};
```

The assistant becomes a chunk of its own, loaded when it is first used. See
[Assistant](/docs/v6/prism-viewer/assistant/).

## Headless

`createPrismViewerHeadless` makes the same viewer with no DOM: the same API, the same agent tools. It is for scripts,
tests and evaluations of the agent.

```ts
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';

const viewer = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec() });

const r = await viewer.api.load.files({ files: [new File([bytes], '1cbs.cif')] });
if (r.ok) console.log(r.did);
await viewer.dispose();
```

| Option    | Meaning                                                                                                                                                                       |
| --------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `spec`    | The spec, as above.                                                                                                                                                           |
| `plugin`  | `(spec: PluginSpec) => PluginContext`. Pass one that makes a `HeadlessPluginContext` (from `@molstar/plugin-headless`) to render images. Without it the plugin has no canvas. |
| `modules` | What the host offers: `'canvas'`, `'files'`, `'network'`, `'storage'`, `'page'`. The default is files and network, and a canvas when `plugin` is given.                       |

An operation that needs a module the host lacks does not run: it fails with the code `unavailable` and says why. The
canvas operations (`canvas.background`, `camera.frame`, `export.image`) are the common case. Headless viewers have no
shell, so `viewer.views.shell` is undefined.
