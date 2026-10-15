# Minimal Prism viewer

The smallest [Prism](../../packages/prism/AGENTS.md) viewer composition: it loads an mmCIF structure (a PDB ID, or the
URL of an `.cif` or `.bcif` file) and shows it as **cartoon** or **ball-and-stick**. A concrete proof that a viewer can
be composed from nothing, and a fixed point to benchmark bundle size and startup against.

```sh
pnpm -C examples/minimal-prism-viewer dev      # http://localhost:5173, from the workspace sources
pnpm -C examples/minimal-prism-viewer build    # vite build → build/
pnpm -C examples/minimal-prism-viewer build:single  # one JS file → build/single/
pnpm -C examples/minimal-prism-viewer preview  # serves build/
pnpm -C examples/minimal-prism-viewer typecheck
```

It is a Vite app (the other examples build with esbuild through `scripts/esbuild/app.mjs`;
`"molstar": { "builder": "vite" }` in the manifest keeps it out of `pnpm build:apps` and `pnpm dev`). `vite.config.ts`
resolves workspace packages to their sources with the `molstar-src` export condition, so nothing needs building first.

## The composition

[`src/spec.ts`](src/spec.ts) is a `PrismViewerSpec` written out from nothing, not `DefaultPrismViewerSpec()` with parts
removed ([viewer-spec.md](../../.v6/prism-ui/spec/viewer-spec.md)):

| Spec field   | Here                                                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------------------------------------ |
| `plugin`     | registry: `Mmcif`, `DefaultHierarchyPresetEntry`, `PolymerCartoonPresetEntry`, `BallAndStickPresetEntry`                       |
| `plugin`     | behaviors: none (no highlight, selection, focus, camera helpers, custom properties, MolViewSpec)                               |
| `extensions` | none: no loaders, canvas settings, measurements, snapshots, animation, export, state tree, superposition                       |
| `sequence`   | `false`: Seq\* is not in the bundle                                                                                            |
| `agent`      | `false`: no tools, no assistant, no agent code (AI SDK, zod)                                                                   |
| `persist`    | none: nothing is stored in the browser                                                                                         |
| `overrides`  | every core tab, rail button and card removed by id; `options.ui.buttons` turns off the pick-mode switch and the trajectory bar |

The page ([`src/index.ts`](src/index.ts), [`index.html`](index.html)) has a header of its own outside the viewer: an ID
or URL field and the two looks. Loading is raw plugin work through `viewer.run` (download, parse as mmCIF, apply the
default hierarchy preset with the chosen representation preset), because the `loaders` feature (`load.entry`) is one of
the things left out. Switching runs the core's `structure.preset` operation (`polymer-cartoon` or `ball-and-stick`).
`performance.mark` records `minimal-prism-viewer:created` and `minimal-prism-viewer:structure-shown` for startup
measurements.

### What it includes

The viewer core (the plugin, the canvas, the Prism shell's frame, the API with all of its core areas, the history), the
mmCIF reader, the two presets and the cartoon and ball-and-stick representations with their default themes, and the RCSB
URL it loads from.

### What it leaves out

Hover highlight, picking and selection, focus, the Home, Structure, Selection, Log, History, Settings and Help tabs, the
toolbar (reset, fullscreen, export, theme, palette), the sequence viewer, every feature and extension, the assistant,
other formats (PDB, SDF, volumes, trajectories), other representations and themes, and persistence.

## Bundle size

Vite 8.3, production, minified, no source maps, 2026-10-16, nothing loaded at runtime. Two builds:

**Single file** (`pnpm build:single`, `vite build --mode single` → `build/single/`): every lazy part inlined into one JS
file, the counterpart of the Prism viewer app's single file (`apps/prism-viewer/build/prism-viewer.js`).

| Build                                          | JS size    | gzip (level 9) | brotli     |
| ---------------------------------------------- | ---------- | -------------- | ---------- |
| minimal example, single file (Vite)            | 3,253.8 kB | 899.9 kB       | 717.5 kB   |
| full Prism viewer, `prism-viewer.js` (esbuild) | 6,243.1 kB | 1,888.5 kB     | 1,498.1 kB |

The minimal composition is 52% of the full viewer's single file (48% gzipped). The bundlers differ (the full app is
built by esbuild without identifier minification), so the ratio is approximate; a like-for-like number needs the full
spec built with Vite too.

**Default build** (`pnpm build` → `build/`, split into chunks): `assets/index-*.js` 3,236.5 kB (894.6 kB gzip, 712.9 kB
brotli) and a 2.0 kB chunk (1.0 kB gzip), `index.html` 1.1 kB.

**Engine alone** (no Prism: a bare `PluginContext` with the same four registry entries, a canvas, and the same load; one
file): 2,187.8 kB, 608.7 kB gzipped, 491.0 kB brotli. So about two thirds of the minimal Prism viewer is the engine, and
the Prism viewer and its shell add about 1.07 MB (about 0.29 MB gzipped).

The page renders 1CBS about 1 s after load on a laptop (`structure-shown` mark: 0.95 s from the single-file preview, 1.1
s from the dev server; network included, a rough first number, not a benchmark).

## Findings: where the viewer API forces more in than needed

Leaving a part out of the spec does not remove it from the bundle when the viewer's core imports it anyway. Measured on
an unminified build (6.0 MB of module code in all):

- **The core API areas are always built** (`CoreAreas` in `packages/prism/viewer/src/api/index.ts`): `volume`,
  `particles`, `trajectory`, `annotations`, `styles`, `select`, `contents` and the rest come with the viewer. Together
  `packages/prism/viewer/src/api` is 214 kB, of which `representation` 28 kB and `volume` 23 kB, and they pull in the
  volume and particle models, transforms and representations that the composition never shows (volume and particle code
  together are about 265 kB unminified).
- **The shell's views come with the shell**: `packages/prism/viewer/src/views` (173 kB, among them the volume transfer
  and levels forms) and `packages/prism/ui/src/views` and `components` (360 kB) are bundled although the overrides
  remove every tab and button; overrides take parts out at run time, not out of the bundle.
- **No way to drop the rail**: with every tab and button removed, the shell still draws its empty 44 px right rail.
  `compactRail` only hides it behind a button in the canvas's corner.
- **MolScript transpilers**: the viewer's `targets/compile.ts` imports `@molstar/query-language/compile`, which
  registers the Jmol, PyMOL and VMD transpilers (118 kB unminified, more with their shared code), although only
  MolScript is needed. The bare engine does not have them.
- **Engine code only the viewer reaches**: about 350 kB more unminified than the bare engine has, from what the viewer's
  API imports: the interaction computations (hydrogen bonds, charges, refinement) and their cylinder representations,
  the volume streaming behavior, isosurface and marching cubes, the focus representation and camera behaviors, the
  valence model.
- **The engine floor** (2,188 kB minified, 609 kB gzipped alone, see below): `@molstar/graphics` is 1.8 MB unminified
  (shaders 300 kB, the render passes 290 kB, geometry 260 kB, every structure visual 230 kB), `@molstar/model`'s
  structure 490 kB, `@molstar/io`'s CIF reader 270 kB. Those are the cost of any Mol\* canvas; narrowing them is the
  engine packages' work, not the viewer spec's.
- **`PrismViewerBehaviors` and `DefaultRegistry` are not the floor**: this example imports neither (`sources/base.ts`
  and `default-spec.ts` are the only modules that do, and nothing in the composition path imports them).
- **`@molstar/ai-fibers`** (44 kB unminified) is in the bundle through the viewer's imports, although `agent: false`.

So the smallest composition today is the viewer core plus the engine, about 3.2 MB (half of the full viewer's single
file). The viewer spec trims the features, extensions, the sequence and the assistant (each its own chunk or off the
bundle), which is where the full viewer's remaining weight is; the core's areas and views would need to become optional
(areas by spec, views behind the shell's parts) to go lower.
