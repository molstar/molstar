# Plan: Backend selection, runtime switching and builds with a subset of backends

## Status

- Planning only; step 1 is part of the second WebGPU milestone of [backend-split.md](backend-split.md) (Order and
  validation). Extends [backend-split.md](backend-split.md) step 6 (context creation and loss), which this plan details.
  Related: [wgsl-strategy.md](wgsl-strategy.md), [neutral-conventions.md](neutral-conventions.md),
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
  single-backend builds?
- Automatic fallback on device loss or missing features, or only explicit switching?
- Should the sync readback APIs (screenshots, `getImageData`) be async everywhere before switching lands (backend-split
  step 5), since WebGPU readback is async only?
- Headless rendering: WebGL in Node only, or a WebGPU backend for Node (e.g. Dawn bindings) later? Answered: both; the
  caller passes Dawn's `webgpu` package as an external module, as `gl` for WebGL.
