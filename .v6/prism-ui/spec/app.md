# The app, the single file, embedding

How the Prism viewer reaches people: the app on GitHub Pages, a single-file viewer a page imports, the `<prism-viewer>`
element, and the URL. Status: **accepted** (2026-10-15).

## 1. `apps/prism-viewer`

- **Private and thin** (README §2.4): `index.html`, `src/index.ts` (the spec: `DefaultPrismViewerSpec()`, the extension
  list, options from the URL, the assistant's import, `persist`; viewer-spec.md §5), a favicon.
- **The extension list** (`MolstarViewerExtensions`) is the app's: the viewer's own features come with
  `DefaultPrismViewerSpec()`; the app adds Mol\*'s extensions, data-only ones through `pluginExtension`, heavy ones lazy
  (video export, Zenodo, the PAE view).
- **Build**: ESM with code splitting, one chunk per lazy feature and extension and one for the assistant, content-hashed
  file names (esbuild: `scripts/esbuild/app.mjs` builds IIFE today and gains a splitting mode; the `molstar`
  distribution's ESM build already has shared chunks). Prism's CSS is in its JS modules, so there is no stylesheet.
- **Deployment**: GitHub Pages, from the release workflow; a version's files are immutable (hashed), `index.html` is
  not cached. Chunks load relative to the entry, so the app works from any base path.
- **Size**: the first chunk is what opening an empty viewer needs; a budget per chunk is set once the first build
  exists and checked in CI.

## 2. The single-file viewer

- **One ESM file** (`prism-viewer.js`) with every lazy part inlined, CSS included: a page imports it and runs it with no
  build step and no other request to us.
- Built from the same spec as the app: `DefaultPrismViewerSpec()` and the app's extension list, with the assistant left
  out (a host that wants it uses the app or the package).
- Published with the app (GitHub Pages) and in the `molstar` distribution.

```html
<script type="module">
  import { createPrismViewer, DefaultPrismViewerSpec } from 'https://…/prism-viewer.js';
  const viewer = await createPrismViewer({ target: document.getElementById('viewer'), spec: DefaultPrismViewerSpec() });
</script>
```

## 3. Embedding

- **`<prism-viewer>`**: the element wraps `createPrismViewer`, with a shadow root (host CSS stops at it; the `--msp-*`
  tokens cross it, README §2.20). Its `spec` property takes a `PrismViewerSpec`; attributes cover the common options
  (`pdb`, `theme`, `hide-controls`). Several on a page each own their viewer.
- **Embeddable parts** for hosts with their own layout: the viewport alone (canvas, hover label, toolbar), for
  mvs-stories and the docking viewer; more when asked.
- **React and Preact hosts** use the element, or `reactView` around it.

## 4. The URL

- The loaders feature reads v5's parameters at start and runs them as API calls with origin `url` (README §3.2):
  `pdb`, `pdb-ihm`, `afdb`, `emdb`, `model-archive`, `url` with `url-format` and `url-is-binary`, `structure-url`,
  `mvs-url`, `mvs-data`, `mvs-format`, `snapshot-url`, `snapshot-url-type`, `snapshot-id`.
- The options read theirs (viewer-spec.md §2.2): `pixel-scale`, `pick-scale`, `pick-padding`, `transparency`,
  `illumination`, `power-preference`, `resolution-mode`, `prefer-webgl1`, `allow-major-performance-caveat`,
  `pdb-provider`, `emdb-provider`, `map-provider`, `hide-controls`, `collapse-left-panel`, `show-toggle-fullscreen`,
  `debug-mode`, `timing-mode`.
- Unknown parameters are ignored and logged once.
