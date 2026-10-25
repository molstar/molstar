# The app, the single file, embedding

How the Prism viewer reaches people: `apps/prism-viewer`, `prism-viewer.js`, `<prism-viewer>`, the URL. The code is the
truth for how.

## 1. `apps/prism-viewer`

- **Private and thin** (README decision 4): the page (`src/index.ts`), its spec (`src/spec.ts`, viewer-spec.md §5), the
  extension list (`src/extensions.ts`, `data-extensions.ts`: Mol\*'s extensions, heavy ones lazy through stubs),
  `molstar.lib` (`src/lib.ts`).
- **On the web, so**: remote states on (the library's default is off), Home at start only when the URL loads nothing,
  the assistant a chunk loaded on first use.
- **Build**: ESM with splitting, content-hashed chunks, no stylesheet (CSS is in the JS), relative chunk URLs so any
  base path works (`scripts/esbuild/app.mjs`). Size budgets (`molstar.budgets`) fail a production build: the entry's on
  its start closure (`index.js` and every chunk it imports statically: what a page loads before the viewer runs; the
  entry itself is a stub of imports), the single file's on the file. A static import of a lazy part fails it.
- **Deployment** (`scripts/distribution.mjs`, `scripts/deploy.js`): GitHub Pages at `prism/`; a publish adds files, so
  earlier hashed chunks stay for open tabs (journal.md §7).
- **`window.molstar = { viewer, api, lib }`**: `lib` is v5's `molstar.lib` where the app has the parts (not the classic
  UI's, not `StateTransforms`: [architecture.md](../../../../.v6/designs/architecture.md) §9.3) plus `lib.prism`. Its
  own lazy chunk; neither `@molstar/prism-viewer` nor the single file exports it.

## 2. The single-file viewer

- **One ES module, `prism-viewer.js`**, every lazy part inlined: a page imports it with no build step and no other
  script request. Images beside it in `images/`, fetched only when used.
- **Built by the app** from `src/single.ts`: exports `createPrismViewer`, `FullPrismViewerSpec()` (the built-ins plus
  the app's extensions, no assistant; named apart from `DefaultPrismViewerSpec()` on purpose) and defines
  `<prism-viewer>`. Shipped with the app and the `molstar` distribution; checked by `smoke/browser/prism`.

## 3. Embedding

- **`<prism-viewer>`** (`@molstar/prism-viewer/element`) wraps `createPrismViewer` in a shadow root: host CSS stops,
  `--msp-*` tokens cross (README decision 20), styles adopt ([library.md](../../ui/_spec/library.md) §3).
- `spec` property (setting it remakes the viewer); attributes `theme`, `hide-controls` (read at creation), `pdb` (a new
  value adds that entry). Each element owns its viewer; leaving the page disposes it, a move keeps it; `ready` (promise
  and event), `error` event.
- Several elements on a page: each its own `persist.key` (README decision 29); `viewer.camera` and `syncCameras` keep
  their cameras together (api.md §12).
- React hosts use the element.

## 4. The URL

The one list. The loaders feature reads the load parameters once the viewer is ready and runs them as the `url` origin,
in this order (`features/loaders/url.ts`; names without case); the app reads the rest at creation
(`apps/prism-viewer/src/url.ts`, `spec.ts`). Unknown parameters are logged once.

- `snapshot-id`; `snapshot-url`, `snapshot-url-type` → `load.snapshot`
- `structure-url`, `structure-url-format` (default `mmcif`), `structure-url-is-binary` → `load.url`
- `mvs-url` or `mvs-data`, `mvs-format` → `load.mvs`
- `pdb`, `pdb-ihm` (`pdb-dev` deprecated), `emdb`, `afdb`, `model-archive` → `load.entry`
- `url`, `url-format`, `url-is-binary` → `load.url`
- `svs-url` → `load.svs` (last)
- `all-models=1`, `assembly=<id>` → options of the `structure-url`/`url` loads (`assembly` also of `pdb`, `pdb-ihm`,
  `model-archive`); alone they load nothing
- `pixel-scale`, `pick-scale`, `transparency`, `resolution-mode`, `power-preference`, `graphics-backend` (`auto`: WebGPU
  where it can, else WebGL; `webgpu`, `webgl`, `webgl2`, `webgl1`), `allow-major-performance-caveat`,
  `disable-antialiasing`, `illumination` → `options.renderer`
- `bg-color` (`viewport-background-color`), `viewport-focus-behavior`, `palette` → `options.canvas`
- `pdb-provider`, `emdb-provider`, `map-provider` (as v5), `volume-streaming-server` (wins over `map-provider`),
  `volume-streaming-disabled`, `state-server` → `options.data`
- `collapse-left-panel` (`hide-controls`) → `options.ui.collapsed`; `show-toggle-fullscreen` → the button, only when
  given (shown by default)
- `hide-rail` → the embedded layout (`ui.layout.compactRail`, [viewer-spec.md](viewer-spec.md) §2.5); `assistant` → the
  Assistant tab open at start; `restore` → this tab's light restore on reload when nothing loads (journal.md §6)
- `debug-mode=1`, `timing-mode=1` → page-wide, before the viewer is made
- `pick-padding` → dropped
