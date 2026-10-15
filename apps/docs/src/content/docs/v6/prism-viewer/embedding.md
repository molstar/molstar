---
title: Embedding the Prism viewer
description: "Create a Prism viewer in a page with createPrismViewer, use the single-file build, configure the spec and options, and run it headless in Node."
---

A page creates a viewer from a spec. The spec is plain data: lists of extensions, options, and what to leave out.

<!-- prettier-ignore-start -->
:::caution[Work in progress]
The packages are not published yet, so the import paths are the workspace's. The `<prism-viewer>` element and viewport-only embedding that the design describes are not built.
:::
<!-- prettier-ignore-end -->

## createPrismViewer

```ts
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewer } from '@molstar/prism-viewer/page';

const viewer = await createPrismViewer({
  target: document.getElementById('viewer')!,
  spec: DefaultPrismViewerSpec(),
});

await viewer.api.load.entry({ source: 'pdb', id: '1tqn' });
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

| Member                                    | What it is                                                                              |
| ----------------------------------------- | --------------------------------------------------------------------------------------- |
| `viewer.api`                              | Every operation, as a script calls them. See [the API](/docs/v6/prism-viewer/api/).     |
| `viewer.run(label, (plugin, api) => ...)` | Raw plugin work as one undoable step.                                                   |
| `viewer.state.ready`, `.busy`, `.empty`   | `ReadonlyState`s: an rxjs `Observable` with a synchronous `value`.                      |
| `viewer.events.hover`, `.click`, `.error` | Observables. Hover and click carry the target under the pointer.                        |
| `viewer.extensions`                       | Which extensions are loaded, the status of each lazy one, and `load(id)`.               |
| `viewer.managers`                         | `commands`, `history`, `navigation`, `viewState`, `notifications`, `transactions`.      |
| `viewer.record`                           | The session record that History shows.                                                  |
| `viewer.plugin`                           | The Mol\* plugin context. An escape hatch that the API does not promise to keep stable. |
| `viewer.dispose()`                        | Tears the viewer down.                                                                  |

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
  await viewer.api.load.entry({ id: '1cbs' });
</script>
```

The module exports `createPrismViewer`, `FullPrismViewerSpec` and `DefaultPrismViewerSpec`. It has no assistant (see
below) and does not read the URL.

`pnpm build:apps` writes it to `apps/prism-viewer/build/prism-viewer.js`. Nothing publishes it yet: it is not deployed
and not part of the `molstar` distribution.

## The spec

`DefaultPrismViewerSpec()` returns a new spec each call. It holds the viewer's base and its own features as extensions:

| Feature           | Brings                                                                                              |
| ----------------- | --------------------------------------------------------------------------------------------------- |
| `loaders`         | `load.entry`, `load.url`, `load.mvs`, `load.snapshot`, `load.trajectory`, `load.volume`.            |
| `canvas-settings` | The canvas settings card and `canvas.lighting`, `.illumination`, `.camera`, `.clipping`.            |
| `measurements`    | The `measure` area and the Measure section.                                                         |
| `superposition`   | The `superpose` area and the Superpose section. Lazy.                                               |
| `state-tree`      | The `state` area and the State tree tab. Lazy.                                                      |
| `remote-states`   | Remote states in the `session` area. Lazy, and left out unless `options.ui.showRemoteStates` is on. |
| `snapshots`       | The `snapshot` area and the Snapshots tab.                                                          |
| `animation`       | The `animation` area and the animation card.                                                        |
| `export-image`    | The `export` area: images.                                                                          |
| `export-session`  | `session.save`.                                                                                     |

`FullPrismViewerSpec()` is the default spec with Mol\*'s extensions added after the built-in ones. It lives in the app
(and in `prism-viewer.js`), because the extensions depend on the viewer. Today it adds
[model quality](/docs/v6/prism-viewer/extensions/#an-example-model-quality) (pLDDT and QMEAN). The app's own spec is the
full spec with the URL's options and the assistant.

A spec has these fields, all but `plugin` optional:

| Field        | What it is                                                                                                                          |
| ------------ | ----------------------------------------------------------------------------------------------------------------------------------- |
| `plugin`     | The engine's base: registry entries, behaviors, config items and canvas props. Named parts in `@molstar/prism-viewer/sources/base`. |
| `options`    | Values for the options below, by group.                                                                                             |
| `extensions` | A list of extensions, in order. Leave one out to drop it. See [Extensions](/docs/v6/prism-viewer/extensions/).                      |
| `overrides`  | Replace or remove any part (a tab, a toolbar button, a command, a sequence behavior) by id.                                         |
| `sequence`   | `{ behaviors, style? }` for the one sequence viewer, `style` being `'tracks'` or `'text'`; `false` for none.                        |
| `agent`      | `{ core?, skills?, assistant? }`, or `false` for no tools at all.                                                                   |
| `ui`         | `layout` (docks, tabs, the tab open at the start), `theme` and `keybindings`. Ignored headless.                                     |
| `persist`    | `{ key, rows? }`. A key names a viewer's stored data, and two live viewers on a page cannot share one.                              |
| `finalize`   | `(spec: PluginSpec) => PluginSpec`: the last word on the composed plugin spec.                                                      |

Overrides name a part as `<owner>.<kind>.<name>` and a command by its own id. `null` removes, an object replaces fields:

```ts
const spec = DefaultPrismViewerSpec();
spec.overrides = {
  'core.toolbar.reset': null,                         // no reset button
  'camera.reset': null,                               // and no reset command
  'measurements.section.measure': { title: 'Distances' },
};
```

An override of an id that nothing registers is an error at creation, naming it.

## Options

`spec.options` has four groups. Each is partial, and a missing value is the default.

| Group      | Values                                                                                                                                                                               |
| ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `renderer` | `pixelScale`, `pickScale`, `transparency`, `resolutionMode`, `powerPreference`, `preferWebgl1`, `allowMajorPerformanceCaveat`, `disableAntialiasing`, `illumination`                 |
| `canvas`   | `background` (a CSS color; empty for the default), `focusBehavior`                                                                                                                   |
| `data`     | `pdbProvider`, `emdbProvider`, `volumeStreamingServer`, `volumeStreamingDisabled`, `stateServer`, `saccharideCompIdMapType`                                                          |
| `ui`       | `showSequence`, `showHistory`, `showRemoteStates`, `collapsed`, and `buttons`: `reset`, `screenshot`, `expand`, `fullscreen`, `settings`, `selectionMode`, `animation`, `trajectory` |

```ts
const viewer = await createPrismViewer({
  target,
  spec: {
    ...FullPrismViewerSpec(),
    options: {
      canvas: { background: '#ffffff' },
      ui: { collapsed: true, buttons: { animation: false } },
    },
    ui: { theme: 'dark' },
  },
});
```

`ui.theme` is `'light'`, `'dark'`, `'blue'`, `'contrast'` or `'auto'`.

## Several viewers on a page

Give each viewer its own `persist.key`, for example `persist: { key: 'left' }`. Keys go to the viewer that has focus,
not to the window.

## The assistant

The assistant is not part of any default spec. The app adds it by passing a lazy import:

```ts
const spec = {
  ...FullPrismViewerSpec(),
  agent: { assistant: () => import('@molstar/prism-assistant/viewer') },
};
```

The import becomes a chunk of its own, loaded after the viewer is up. See [Assistant](/docs/v6/prism-viewer/assistant/).

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
