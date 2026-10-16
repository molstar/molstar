# The app, the single file, embedding

How the Prism viewer reaches people: the app on GitHub Pages, a single-file viewer a page imports, the `<prism-viewer>`
element, and the URL. Status: **accepted** (2026-10-15).

## 1. `apps/prism-viewer`

- **Private and thin** (README §2.4): `index.html`, `src/index.ts` (the page: the viewer, the URL's loads or the light
  restore), `src/spec.ts` (the spec: `FullPrismViewerSpec()`, options from the URL, the assistant's import, `persist`;
  viewer-spec.md §5), `src/extensions.ts` and `src/data-extensions.ts` (the extension list), a favicon.
- **The extension list** (`MolstarViewerExtensions`) is the app's: the viewer's own features come with
  `DefaultPrismViewerSpec()`; the app adds Mol\*'s extensions, the data-only ones as plain `{ id, title, plugin }`
  entries, heavy ones lazy (video export, Zenodo; the PAE view when it is built).
- **On the web** (2026-10-16): remote states are on (`ui.showRemoteStates`: Home's Remote states section, the card,
  Save and Open on the viewer's `state-server`); the library's default leaves them off for pages that embed it.
- **Build**: ESM with code splitting, one chunk per lazy feature and extension and one for the assistant, content-hashed
  file names (esbuild: `scripts/esbuild/app.mjs` builds IIFE today and gains a splitting mode; the `molstar`
  distribution's ESM build already has shared chunks). Prism's CSS is in its JS modules, so there is no stylesheet.
- **Deployment**: GitHub Pages, from the release workflow, on an origin of its own (a custom domain or a project
  site whose origin serves nothing else), since every page on an origin can use its stored keys and kept records
  (agent.md §4). A publish adds files rather than replacing the site, so the hashed chunks of the last few versions
  stay for tabs left open across a deploy (journal.md §7); `index.html` is not cached. Chunks load relative to the
  entry, so the app works from any base path.
- **Page globals**: `window.molstar = { viewer, api, lib }` once the viewer is made. `lib` (`src/lib.ts`) is v5's
  `molstar.lib` where the app has its parts (`structure`, `volume`, `shape`, `loci`, `math`, `plugin`, `extensions`;
  not the classic UI's `PluginUIContext`, `DefaultPluginUISpec` and view models, nor `StateTransforms`, which
  architecture.md §9.3 keeps for the classic Viewer only) and `lib.prism`: `createPrismViewer`,
  `DefaultPrismViewerSpec`, `FullPrismViewerSpec`, `node`, `MS` (the MolScript builder). It is the app's, in a chunk
  loaded beside the viewer (115 kB, none of it in the start bundle); `@molstar/prism-viewer` stays tree-shakable and
  the single file does not export it (it would add 32 kB gzipped to every page that imports it).
- **Size**: the first chunk is what opening an empty viewer needs; a budget per file (`molstar.budgets` in the app's
  manifest: `index.js` 3.2 MB, `prism-viewer.js` 8 MB, minified) fails a production build that outgrows it. Since
  `lib`'s import shares the start code, esbuild puts that code in shared chunks and `index.js` is a 6 kB list of
  imports: the start bundle is its static-import closure (65 files, 4.4 MB, 1.29 MB gzipped, 2026-10-16), which the
  `index.js` budget no longer bounds. Not yet: a budget over that closure, hashed entry names, the assistant's chunk
  loaded on first use rather than at start.

## 2. The single-file viewer

- **One ESM file** (`prism-viewer.js`) with every lazy part inlined, CSS included: a page imports it and runs it with no
  build step and no other request to us.
- **An artifact of the app** (review F27): `apps/prism-viewer` builds it from its own entry, which exports
  `createPrismViewer` and `FullPrismViewerSpec()`: `DefaultPrismViewerSpec()` with the app's extension list (Mol\*'s
  extensions) and without the assistant (a host that wants it uses the app or composes the packages). The name differs
  from `DefaultPrismViewerSpec()` (the built-ins only) on purpose.
- Published with the app (GitHub Pages) and copied into the `molstar` distribution (`build/prism-viewer/`, the file and
  its `images/`, as
  `build/viewer` is today; `scripts/distribution.mjs`), checked by the packed browser smoke (`smoke/browser/prism`: one
  script request, a load, a lazy extension). Publishing to Pages: not yet. Its images (backgrounds) sit beside it in `images/`, by URLs relative to the module, fetched only when a background uses them (inlined, they were 0.8 MB of every page); the agent's code (AI SDK, zod) comes only with an assistant's module, so the file has none. A host composing its own viewer from npm imports the extensions it wants instead.

```html
<script type="module">
  import { createPrismViewer, FullPrismViewerSpec } from 'https://…/prism-viewer.js';
  const viewer = await createPrismViewer({ target: document.getElementById('viewer'), spec: FullPrismViewerSpec() });
</script>
```

## 3. Embedding

- **`<prism-viewer>`**: the element wraps `createPrismViewer`, with a shadow root (host CSS stops at it; the `--msp-*`
  tokens cross it, README §2.20; the style registry adopts into it, library.md §3). Its `spec` property takes a `PrismViewerSpec`; attributes cover the common options
  (`pdb`, `theme`, `hide-controls`). Several on a page each own their viewer. Built: `@molstar/prism-viewer/element`
  (`definePrismViewerElement`, `createPrismViewerElementClass`): the element is `createPrismViewer`'s target (its
  shadow root holds the shell); `theme` and `hide-controls` are read when the viewer is made, `pdb` loads its entry
  (`load.entry`; a new value adds that entry); setting `spec` makes the viewer again; leaving the page disposes it,
  a move within the page keeps it; `ready` (a promise and an event) gives the viewer. The single file defines it once,
  with `FullPrismViewerSpec()` while `spec` is not set.
- **Embeddable parts** for hosts with their own layout: the viewport alone (canvas, hover label, toolbar), for
  mvs-stories and the docking viewer; more when asked.
- **React and Preact hosts** use the element, or `reactView` around it.

## 4. The URL

This is the one list (review F26). The loaders feature reads the load parameters, once the viewer is ready; the options
read the rest at creation.

| Parameter                                                                                  | Read by                         | Becomes                                                           |
| ------------------------------------------------------------------------------------------ | ------------------------------- | ----------------------------------------------------------------- |
| `pdb`, `pdb-ihm` (`pdb-dev`, deprecated, the same), `afdb`, `emdb`, `model-archive`        | loaders                         | `load.entry`                                                      |
| `url` with `url-format`, `url-is-binary`; `structure-url` with `structure-url-format`, `structure-url-is-binary` | loaders  | `load.url`                                                        |
| `all-models=1`, `assembly=<id>` (Prism's own; v5 has them as `loadAllModelsOrAssemblyFromUrl`) | loaders | `models: 'all'` and `assembly` of the `structure-url` and `url` loads; `assembly` of the `pdb`, `pdb-ihm`, `model-archive` entries; alone they load nothing |
| `mvs-url`, `mvs-data`, `mvs-format`                                                        | loaders                         | `load.mvs`                                                        |
| `snapshot-url`, `snapshot-url-type`, `snapshot-id`                                         | loaders                         | `load.snapshot`                                                   |
| `pixel-scale`, `pick-scale`, `transparency`, `illumination`, `power-preference`, `resolution-mode`, `prefer-webgl1`, `allow-major-performance-caveat` | options | `options.renderer`                              |
| `pdb-provider`, `emdb-provider`, `map-provider`                                            | options                         | `options.data` (`map-provider` sets the volume streaming server, as v5 does: `rcsb` is `https://maps.rcsb.org`, `pdbe` PDBe's; `volume-streaming-server` wins) |
| `hide-controls`, `collapse-left-panel`, `show-toggle-fullscreen`                          | options                         | `options.ui` (`collapsed`; `buttons.fullscreen` only when the parameter is given: Prism shows the button, v5 hid it unless `=1`) |
| `debug-mode`, `timing-mode`                                                                | the app                         | page-wide switches (`=1`: `setDebugMode`, `setTimingMode` before the viewer is created; `modesFromUrl`), not per viewer |
| `pick-padding`                                                                             | —                               | dropped (no config item; v5 drops it too)                         |

Unknown parameters are ignored and logged once.
