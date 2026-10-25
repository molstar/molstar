# Plan: Backend selection, runtime switching and builds with a subset of backends

## Status

- **Current**: Default backend decided and done (`auto`, see Default backend): the Viewer's `graphics-backend` and the
  plugin default (`PluginConfig.General.GraphicsBackends` unset, `@molstar/graphics/gl/backends/default` the `all`
  set). Steps 1
  (descriptors, async creation, the plugin config item), 2 (import hygiene: the guard test, per-backend program
  registries of [glsl-in-webgl.md](glsl-in-webgl.md)), 3 (GPU compute behind the backends, chosen by capability:
  [backend-narrowing.md](backend-narrowing.md) step 3) and 4 (backend sets, the bundle check and Viewer variants with
  one backend, below) done. Steps 5–8 (runtime switching) done, see Runtime switching status.
- **Runtime switching status**: `Canvas3DContext.setBackend(backends)` creates the context of the first backend that can
  be created on a new canvas element (each backend on its own fresh canvas), which takes the place of the old one
  (attributes and size copied, the input observer moves its listeners with `InputObserver.setElement`, its subjects
  stay), replaces `graphics`, `canvas` and `passes`, emits `contextChanged` while the previous context is alive, then
  destroys it; it refuses while presenting XR. `Canvas3D` keeps its identity and rebuilds on it: the scene, helpers,
  camera, controls, interaction and XR manager stay (`Scene.setContext`, the helpers' and `XRManager.setContext`), the
  renderer, Hi-Z, shader manager, pick, ray and multi-sample helpers are created again with their props. The render
  objects stay (their values, markers included): textures loaded from data are created again on the new context from the
  data of their last load (`Texture.getSource`, `transferTexture`, `transferTextures`); objects with rendered or
  computed textures (GPU surfaces, GPU colors) are left out until their representation updated them: representations
  continue with the new context (`Representation.setGraphics`), Canvas3D updates those with left-out objects, and
  visuals whose render object holds textures or storage buffers of another context (`hasForeignResources`, by
  `resources.ownsTexture` and `compute.ownsBuffer`)
  create their geometry again (the same render object); caches skip textures of another context (texture-mesh double
  buffers with their draw arguments, color smoothing). Image passes created by `getImagePass` follow
  (`ImagePass.setContext`). Nothing to do for apps: the plugin only offers `setGraphicsBackends` and a Graphics Backend
  setting in the viewport settings where it was created with more than one backend (`graphicsBackends`). Context-keyed
  caches: the named textures and framebuffers live on the context, the volume storage buffer is keyed by its compute
  context, the GPU isosurface volume texture is now keyed by the context. Validation:
  `pnpm test:viewer-smoke --backends auto --switch` switches WebGL2, WebGPU, WebGL, WebGL2 at runtime: each renders with
  the coverage of the first render, hovering the grid gives the same labels, and switched to WebGL2 a second time the
  context holds as many resources as after the first switch. Found on the way: WebGL2 uniform blocks missed the
  per-frame updates of window rectangles (the flipped cells were observed instead of their sources), so passes drew with
  the viewport of their creation until another member of the block changed; SMAA drew a blank image when the canvas was
  not resized after the context was created (fixed: the flipped and depth-range cells forward observers to their
  sources, `pixels` rectangles are volatile and checked on every write; unit test). And the SSAO pass created at its
  final size kept the input sizes in its level copies and blur instead of those of its targets (striped occlusion after
  a switch, until the next resize); its constructor ends with `setSize` now. Checked by rendering after a switch with
  SSAO, SMAA, outlines, shadows, DOF and bloom, then again after resizing all passes and back: the same image on WebGL2,
  WebGL and WebGPU. Pending pick reads across a switch (hovering while the camera moves, e.g. the rock animation): a
  read started by the previous pick helper resolved through it, drawing the scene (already on the new context) with the
  previous renderer, every frame; disposed pick and ray helpers now resolve their pending reads to nothing and the
  interaction helper requests them again (`resetPicks`). Not done: context restore on top of the rebuild (the restore in
  place stays), the visual-harness switching case (the smoke test covers it), automatic fallback on device loss; the
  sb-ncbr tunnels keep the context in their transform params (computed with the previous context after a switch).
- Step 4 done: the only import that put WebGL into bundles without it was `Canvas3DContext.fromCanvas`, now
  `createWebGLCanvas3DContext` in `@molstar/graphics/canvas3d/webgl` (the Canvas3D layer's WebGL convenience; the `gl`
  layer doesn't import `canvas3d`), with `Canvas3DContext.fromGraphics` public. Backend sets
  `@molstar/graphics/gl/backends/{webgl,webgpu,all}` export `getDefaultGraphicsBackends`, `.../default` re-exports the
  `all` set and the plugin takes its default from it, so a lean build aliases one import (esbuild `alias`). No
  `sideEffects` needed: esbuild drops the other backend without it. `pnpm check:backend-bundles`
  (`scripts/check-backend-bundles.mjs`, not in CI) builds Canvas3D with each backend and the plugin with the default and
  with it aliased to each backend, and fails on code of an excluded backend (GLSL and WGSL sources included), naming the
  import: Canvas3D with WebGPU 1504 KiB (was 1865 KiB with 353 KiB of WebGL), with WebGL 1585 KiB, the plugin 2003 KiB
  with WebGL and 1923 KiB aliased to WebGPU.
- Viewer variants: the Viewer takes its backends from the backend sets (`graphics-backend`: `webgpu` the WebGPU and
  WebGL sets, `auto` the `all` set (the plugin default), `webgl` the WebGL set), and `molstar.variants` in its manifest (`webgl`,
  `webgpu`) builds it with one backend (`node scripts/esbuild/app.mjs apps/viewer --variants` into `build/<variant>/`,
  the other set modules aliased to the kept one, `scripts/esbuild/backend-sets.mjs`): 4.55 MB (WebGL) and 4.47 MB
  (WebGPU) against 4.83 MB, no code of the other backend (`pnpm check:backend-bundles`);
  `pnpm test:viewer-smoke --variant webgpu` checks one. The import guard resolves imports by path (relative ones through
  `gl/` too, side effect and dynamic imports).
- Started as the second WebGPU milestone of [backend-split.md](backend-split.md) (Order and validation). Extends
  [backend-split.md](backend-split.md) step 6 (context creation and loss), which this plan details. Related:
  [wgsl-strategy.md](wgsl-strategy.md), [neutral-conventions.md](neutral-conventions.md),
  [core-gpu-compute.md](core-gpu-compute.md).

## Goals

1. **Runtime switching**: change the render backend (e.g. WebGL2 → WebGPU, or back after a WebGPU device loss) without
   reloading, keeping the plugin state (structures, representations, camera, selection, props).
2. **Builds with a subset of backends**: consumers such as the Viewer can be built with WebGL only, WebGPU only, or
   both, and the excluded backend's code (context, resources, shader sources, compute paths) is not in the bundle.

## Current state

- A canvas element gives out one context type for its lifetime (`webgl`, `webgl2` and `webgpu` exclude each other), so
  every switch needs a new canvas element.
- `Canvas3DContext.fromCanvas` is synchronous and WebGL-only: `canvas3d.ts` statically imports `gl/webgl/context`,
  `Canvas3DContext.webgl` is a `WebGLContext`; WebGPU creation is async (`requestAdapter`/`requestDevice`).
- Renderables import their `ShaderCode` (GLSL sources) from `gl/shader-code.ts`, so all GLSL is reachable from neutral
  code; WGSL sources are only reachable from `gl/webgpu/` (registry by program name).
- WebGL-typed modules outside `gl/webgl/`: 20 files in `canvas3d/passes`, 6 in `canvas3d/helper`, GPU compute
  (`gl/compute`, gaussian density, texture mesh), plugin core/ui/headless and extensions (geo-export, alpha-orbitals,
  anvil, assembly-symmetry, sb-ncbr tunnels). Most are tracked by backend-split step 5 (readback and sync), the guard
  test lists them per file.
- The plugin creates `canvas3dContext` and `canvas3d` once (`_initViewer`), both are readonly, and behaviors subscribe
  to `canvas3d` events for the lifetime of the plugin.
- Context loss and restore exist for WebGL: resources reset and re-upload from their retained data, render objects with
  GPU-derived values (texture meshes from GPU isosurfaces and Gaussian surfaces) recompute through `values.meta.reset`,
  helpers reset (`pick`, `ray`, `hi-z`).
- Caches holding backend objects outside the scene: `namedTextures`/`namedFramebuffers`/`namedComputeRenderables` on the
  context, volume textures in `Volume._propertyData` (GPU isosurface, segments), asset-backed textures, representation
  contexts passing the context to visuals.

## Design

### Backend descriptors, injected instead of imported

```ts
interface GraphicsBackend {
    readonly id: string; // 'webgl', 'webgl2', 'webgpu'
    /** Cheap availability check without creating a context (e.g. `navigator.gpu`, WebGL2 constructor) */
    isAvailable(): boolean;
    /** Creates a context on a fresh canvas; resolves to undefined if the device/context can't be created */
    create(canvas: HTMLCanvasElement, attribs: GraphicsAttribs): Promise<GraphicsContext | undefined>;
}
```

- `gl/webgl/backend.ts` exports `WebGL2Backend` and `WebGLBackend`, `gl/webgpu/backend.ts` exports `WebGPUBackend`;
  nothing in neutral code imports them.
- `Canvas3DContext.create(container, { backends, preference })` (async) tries the backends in order and owns the canvas
  element (it creates it, so it can replace it). `Canvas3DContext.graphics: GraphicsContext` replaces `webgl`; the
  WebGL-only paths narrow with `isWebGLContext` as today.
- Consumers choose the backend list: the plugin from a config item (`PluginConfig.Graphics.Backends`, default from the
  build entry), the visual test harness per profile, headless (Node) uses the WebGL backend of the `gl` package.

### Neutral code reaches no backend code

- **Program sources per backend**: renderables reference a neutral program description (name, shader extensions, output
  types, ignored defines) instead of `ShaderCode` with sources; `gl/webgl/programs.ts` maps names to GLSL (all
  `gl/shader/*.glsl.ts` imports move behind it), `gl/webgpu/shader/programs.ts` to WGSL (exists). The same for compute
  programs. Detailed in [glsl-in-webgl.md](glsl-in-webgl.md).
- **GPU compute as backend services**: GPU isosurface, Gaussian density, histogram pyramid, color smoothing, grid3d are
  provided by the backend (`graphics.compute?.extractIsosurface(...)`) instead of being imported by representations;
  representations fall back to the CPU paths when a service is missing. This is also the direction of
  [core-gpu-compute.md](core-gpu-compute.md) and wgsl-strategy decision 6 (compute planned).
- **Passes** take a `GraphicsContext` (backend-split step 5) and use no backend imports.
- **Guard test**: extend `backend-neutral.test.ts` so modules outside `gl/webgl/` and `gl/webgpu/` don't import either
  folder (per-file allowlist from today's list, shrinking to the backend entry files), and backend folders don't import
  each other.

### Builds

- Entry modules per backend set, e.g. `@molstar/graphics/gl/backends/webgl`, `.../webgpu`, `.../all`, each exporting the
  backend list; an app imports one (the Viewer passes it into the plugin spec), lean builds use a different entry file
  or an esbuild alias for the same import.
- `sideEffects: false` (or an explicit list) in the package manifests so unreferenced modules drop out.
- Bundle check: build test bundles for each entry with an esbuild metafile and assert the module sets (WebGPU-only: no
  `gl/webgl/`, no GLSL sources; WebGL-only: no `gl/webgpu/`, no WGSL), report sizes.

### Runtime switching

`Canvas3D` keeps its identity (the object the plugin and behaviors hold) and rebuilds its backend-bound internals:

1. Create the new context on a new canvas element (async); on failure keep the current backend and report it.
2. Pause rendering and the commit queue.
3. Release the old internals: renderables (render items) of all render objects, passes, renderer, pick/ray/hi-z helpers,
   context caches; dispose the old context and canvas.
4. Reset backend-bound values in render objects through the existing `meta.reset` mechanism (GPU-derived geometry
   recomputes with the new context); invalidate context-keyed caches (volume textures keyed by context instead of stored
   per volume, named resources).
5. Rebuild the internals on the new context: passes, renderer, scene renderables from the same render objects, helpers;
   the input observer listens on the container (or is rebound) so its subjects survive.
6. Re-apply camera, props and pixel scale, request a draw, emit `canvas3d.contextChanged`; consumers that hold the
   context beyond a call (representation contexts, plugin transforms passing `graphics` to visuals) read it from
   `canvas3d` at use time.

- Internally this means splitting `Canvas3D` into a stable shell (subjects, props, camera, scene of render objects,
  interaction) and a backend runtime (context, renderables, passes, helpers). The same split simplifies context restore,
  which becomes "rebuild the runtime on the restored context".
- WebGL1 ↔ WebGL2 switching exercises the whole mechanism before WebGPU reaches parity (it also needs a new canvas).
- Optional automatic fallback: on WebGPU device loss, or when a scene needs a feature the active backend lacks, switch
  to the next backend of the list and notify.

## Default backend

Recommendation for the open question of the Viewer's default backend, now that WebGPU renders everything the visual
tests cover.

### Where WebGPU stands

- **Rendering**: all 118 visual cases pass on WebGPU with none skipped (the WebGL profiles skip only the compute
  checks), most with 0 differing pixels against WebGL; the Viewer smoke test passes on both backends, including a device
  loss restored in place ([webgpu-viewer.md](webgpu-viewer.md)); headless rendering takes Dawn's Node.js bindings.
- **GPU compute**: Gaussian density, marching cubes (drawn indirectly, no count round trip), color smoothing, volume
  segments and grid3d (alpha orbitals) run as compute kernels ([webgpu-compute.md](webgpu-compute.md),
  [indirect-draws.md](indirect-draws.md)); checks assert that the GPU paths are taken and match the CPU.
- **Speed** ([marching-cubes-benchmark.md](marching-cubes-benchmark.md), median ms):

  | Case                                                                      | WebGL | WebGPU |
  | ------------------------------------------------------------------------- | ----- | ------ |
  | Gaussian surface representation of 7qpd (11.5 k atoms, 7 chains), created | 35    | 39     |
  | the same, updated                                                         | 31    | 29     |
  | volume isosurface representation (128³), created                          | 17    | 9.2    |
  | Gaussian surface (100 k atom globule at 0.5 Å)                            | 76    | 36     |

  First runs (shader and pipeline compilation) are shorter on WebGPU (about 100–150 ms against 250–320 ms on WebGL).

- **Rendering** (`render/*` cases of `pnpm benchmark`: 30 forced frames at 1280×720 issued back to back, then waited
  for; frame time = total / 30, `submit` = the CPU time to issue a frame; RTX 4070 Ti, Chrome, ANGLE D3D11; median ms):

  | Case                                                                                      | WebGL2  | WebGL2 without UBOs | WebGL1  | WebGPU            |
  | ----------------------------------------------------------------------------------------- | ------- | ------------------- | ------- | ----------------- |
  | 1000 boxes, each its own representation (draw-call bound), frame                          | 5.5     | 6.1                 | 6.1     | 7.3 (was 10.8)    |
  | the same, submit (CPU)                                                                    | 4.4     | 6.1                 | 6.0     | 7.1 (was 10.5)    |
  | 100 boxes, frame                                                                          | 0.5–0.7 | 1.0                 | 0.6–0.9 | 0.7 (was 1.0–1.1) |
  | 300 k atoms, spacefill impostors (fragment bound)                                         | 4.3     | 4.3                 | 4.3     | 4.6               |
  | 300 k atoms, spacefill meshes (vertex bound)                                              | 12.9    | 12.9                | 13.0    | 10.9              |
  | 7qpd cartoon, ball-and-stick, spacefill, Gaussian surface; cartoon with SSAO and outlines | 0.3–0.9 | 0.2–0.8             | 0.2–0.7 | 0.3–0.9           |

  Small scenes (up to about 25 draw calls) are within noise on all backends (sub-millisecond, run-to-run spread from GPU
  clocks larger than the differences). GPU-bound frames are on par, WebGPU 16% faster for vertex-heavy meshes and 7%
  slower for sphere impostors. The gap is the CPU cost per draw call: about 7 µs on WebGPU against about 4.4 µs on
  WebGL2 (10.5 µs before [webgpu-draw-overhead.md](webgpu-draw-overhead.md)), so scenes of many separately drawn objects
  render at about three quarters of WebGL2's frame rate on WebGPU.

### Gaps that users would notice

| Gap                                                                                                            | Effect on WebGPU                                                                                                                  | Fix                                                                                         |
| -------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- |
| Geometry export of GPU surfaces (`geo-export` reads texture meshes back with WebGL)                            | throws; a regression since GPU surfaces became the default on WebGPU (they were CPU meshes before)                                | async texture readback (`readTextureAsync`) in the exporter                                 |
| Vertex-granularity colors on GPU surfaces (the position iterator of texture meshes reads them back with WebGL) | themes that sample at vertex positions (e.g. coloring a surface by an external volume) get no locations                           | the CPU surface when the theme needs positions, or async readback                           |
| Device limits                                                                                                  | done: the adapter's maxima are requested (were the spec defaults: 2d textures of 8192, storage bindings of 128 MiB)               | [webgpu-device-limits.md](webgpu-device-limits.md)                                          |
| Frame times                                                                                                    | on par for typical and GPU-bound scenes; about 1.6× the CPU time per draw call (scenes of hundreds of separately drawn objects)   | [webgpu-draw-overhead.md](webgpu-draw-overhead.md): done steps 1–4, candidates for the rest |
| Runtime switching                                                                                              | done: `Canvas3DContext.setBackend`, `PluginContext.setGraphicsBackends`, the viewport setting; a lost device is replaced in place | Runtime switching status above                                                              |

Synchronous readbacks are already handled: picking and screenshots are asynchronous on WebGPU (`identify` finds nothing
until a read has finished).

### Browsers

WebGPU is enabled by default in Chrome and Edge (Windows, macOS, ChromeOS, Android) and in recent Safari (macOS, iOS);
Firefox ships it on Windows first, Linux support varies by browser and driver. This changes quickly: check the current
state when deciding, and run the Viewer smoke test on each browser available (`--channel`). Where `navigator.gpu` is
missing or no adapter is returned, the descriptors fall back to WebGL2, then WebGL, at creation.

### Options

- **A. WebGL by default, WebGPU opt-in** (today: `graphics-backend=webgpu`): no risk, but WebGPU gets little use and its
  advantages (compute surfaces, indirect draws) stay unused.
- **B. WebGPU where available** (`[WebGPUBackend, WebGL2Backend, WebGLBackend]`), `graphics-backend=webgl` to opt out:
  the faster compute paths for everyone with WebGPU; any WebGPU-only bug reaches all of them.
- **C. Staged**: B in the Viewer first (an app, easy to override and to revert), the plugin and library default
  (`PluginConfig.General.GraphicsBackends`) stays WebGL for a release, then follows; headless stays WebGL by default
  (the `gl` package, Dawn optional).

Recommendation: C, once the blockers below are done.

### Blockers and steps

1. **Export of GPU surfaces** (done, [backend-narrowing.md](backend-narrowing.md) step 2): geo-export reads with
   `readTextureData`; `gpu-paths/export-gpu-surface` compares with the CPU surface on both backends.
2. **Vertex-granularity colors** (done, same step): the vertices are read back before the colors
   (`TextureMesh.withVertexData`); `gpu-paths/vertex-colors-gpu-surface`.
3. **Render benchmark** (measured): frame times while the camera rotates 2° per frame (`render/*-rotating`, 1280×720,
   RTX 4070 Ti, Chrome 154, median of 9 runs of 30 frames, two rounds), WebGL2 against WebGPU:

   | Case                                                                        | WebGL2       | WebGPU                |
   | --------------------------------------------------------------------------- | ------------ | --------------------- |
   | 7qpd as cartoon, ball-and-stick, molecular (alpha 0.5) and Gaussian surface | 0.94–0.99 ms | 1.03–1.04 ms (+5–10%) |
   | 60 instanced copies of 7qpd as cartoon                                      | 3.69 ms      | 3.86 ms (+5%)         |
   | the same as spacefill                                                       | 0.98–0.99 ms | 1.15 ms (+17%)        |

   The switch needs WebGPU within 10% or faster: met for the Viewer scene and cartoon, not for sphere impostors (the GPU
   time of the sphere shader, which [webgpu-draw-overhead.md](webgpu-draw-overhead.md) traces to Chrome's WGSL
   translation; the CPU part is 0.06 ms higher). GPU surface creation and updates are about 2× faster on WebGPU (Speed
   above), which the frame times do not show. Decided: switch with the known sphere gap (the Viewer default is `auto`,
   step 6). The switch makes the Viewer's `graphics-backend=auto` (`AutoWebGPUBackend`, then WebGL2 and WebGL) the
   default: WebGPU only where the adapter can blend into float32 targets (`float32-blendable`; Firefox cannot yet, so it
   stays on WebGL2, which has wboit, dpoit and illumination there); `graphics-backend=webgpu` asks for WebGPU
   regardless.

4. **Adapter limits** (done, [webgpu-device-limits.md](webgpu-device-limits.md)): devices request the adapter's maxima;
   a `webgpu-core` visual profile keeps the minimum spec covered.
5. **Browser matrix** (done on Windows): `pnpm test:viewer-smoke` (`--channel`, `--browser`) passes on Chrome 154,
   Chrome Beta 156 and Edge 155 with both backends (same hover hits, context loss restored), and on Playwright's Firefox
   155 and WebKit with WebGL, where a requested WebGPU backend falls back to WebGL (no adapter there) and works the
   same; the visual suite (`webgl2`, `webgpu`, `webgpu-core`) passes on Chrome Beta and Edge. Not covered: Firefox's and
   Safari's own WebGPU (Playwright's builds have no adapter), macOS/Metal, Linux, integrated GPUs.
6. **Switch** (done for the Viewer): the Viewer's `graphics-backend` defaults to `auto` (WebGPU where the adapter has
   `float32-blendable`, else WebGL2, then WebGL), `graphics-backend=webgl` to opt out, recorded in the v6 changelog; the
   plugin default (`PluginConfig.General.GraphicsBackends` unset) followed: the `all` set. `PreferWebGl1` is removed: WebGL1 by listing
   only its backend. Headless stays WebGL.

## Steps

1. **Backend descriptors and async creation** (partly done: `GraphicsBackendDescriptor` and `createGraphicsContext` in
   `gl/backend.ts`, the name `GraphicsBackend` being the backend id type of the capabilities;
   `WebGL2Backend`/`WebGLBackend` in `gl/webgl/backend.ts`, `WebGPUBackend` in `gl/webgpu/backend.ts`; async
   `Canvas3DContext.create(canvas, assetManager, backends, attribs, props)` and `Canvas3DContext.graphics`, both
   `create` and `fromCanvas` share `fromGraphics`, which accepts only WebGL contexts until the passes take a
   `GraphicsContext`; the visual test profiles create their contexts through the descriptors. Open: `fromCanvas` moving
   to a WebGL module, the plugin config item): `GraphicsBackend`, WebGL/WebGL2/WebGPU descriptors, async
   `Canvas3DContext.create` owning the canvas, `Canvas3DContext.graphics`; `fromCanvas` stays as a WebGL convenience in
   a WebGL module (breaking change recorded). The visual test harness creates its profiles through the descriptors.
2. **Import hygiene**: guard test for backend imports; program descriptions without sources and per-backend source
   registries; passes on `GraphicsContext` (with backend-split step 5).
3. **GPU compute services**: compute entry points move behind the backends, CPU fallbacks where a backend lacks them.
4. **Build entries and bundle check**: backend set entries, `sideEffects`, metafile assertions for WebGL-only,
   WebGPU-only and both; Viewer variants (default and lean builds).
5. **Canvas3D shell and runtime**: split, context restore reimplemented on top of the runtime rebuild; visual tests
   unchanged.
6. **Context-keyed caches**: audit and fix (volume property textures, named resources, asset textures, representation
   contexts, plugin and extension caches).
7. **Runtime switching**: `canvas3dContext.setBackend(id)` with canvas replacement and `contextChanged`; plugin config
   and a viewport setting; WebGL1 ↔ WebGL2 first, WebGPU once its passes exist; optional automatic fallback.
8. **Validation**: a switching test in the visual harness (render, switch, render, compare against a fresh instance,
   switch back) and resource counts back to the baseline after switching (no leaks of the old context's objects).

## Order

- Steps 1, 2 and 4 are cheap and keep regressions out early (bundle check, import guard); step 3 follows the compute
  work.
- Steps 5-7 are useful with WebGL1 ↔ WebGL2 already; switching to WebGPU is only useful for users once WebGPU renders
  the passes (render targets, postprocessing, OIT, picking).

## Open questions

- Default backend set of the Viewer bundle: all with WebGL preferred until WebGPU reaches parity, plus lean
  single-backend builds? Decided: `auto` (WebGPU where capable) for the Viewer and the
  plugin default, single-backend builds by aliasing the backend sets.
- Automatic fallback on device loss or missing features, or only explicit switching?
- Should the sync readback APIs (screenshots, `getImageData`) be async everywhere before switching lands (backend-split
  step 5), since WebGPU readback is async only?
- Headless rendering: WebGL in Node only, or a WebGPU backend for Node (e.g. Dawn bindings) later? Answered: both; the
  caller passes Dawn's `webgpu` package as an external module, as `gl` for WebGL.
