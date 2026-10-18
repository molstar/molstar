# Plan: Resource readiness (async loads and program compilation)

## Status

- **Current**: steps 1–4 done; open: 5 (optional). `ImagePass` rendering once moved to
  [pipeline-readiness.md](pipeline-readiness.md) (prepare from params without rendering; a dry run of the passes in
  collect mode reached it, 63 of 142 WebGPU visual cases rendering again to none, but was not kept).
- `ImagePass` updates the background before rendering (`updateBackground`): `DrawPass` chose a transparent clear from
  the background of the previous render (`background.isEnabled` reads the loaded skybox or image, which the
  postprocessing pass updated only later in the render), so an image after one with a skybox lost its opaque background
  (35 WebGL images across profiles in full visual runs, which matched their isolated runs). The warm-up is to be revisited
  for better coverage: it prepares hover marking only, while first uses of other features still compile on first draw
  (the first selection's programs, illumination, transparency modes, postprocessing; see
  [program-switch-frames.md](program-switch-frames.md)).
- Steps 3 and 4 (the rest): collect mode and the marking warm-up. `resources.preparePrograms(run)` runs `run` with draws
  only requesting their pipelines (WebGPU: no draws, no clears, no timer passes, the drawing buffer's formats without
  taking the canvas texture) and resolves when the pending work settled; `run` is synchronous (the mode holds for the
  whole context). Canvas3D's idle warm-up prepares the draw pass as if objects were marked (`DrawPass.prepareMarking`:
  the scene offscreen, the marking depth and mask passes over all objects, the edge, overlay, composite and copy quads),
  so on WebGPU the first hover compiled 7 pipelines and showed its highlight about 200 ms after the move (160 ms
  compiling), and compiles none now; the Viewer smoke test checks that hovering compiles no pipelines (fails without the
  warm-up). The first selection still compiles the pipelines of the programs its define changes create, which their
  families' known configurations start at creation (4b); the structure stays drawn. WebGL tracks the links
  `linkPrograms` starts as pending work until they finished (`KHR_parallel_shader_compile` `COMPLETION_STATUS`, polled),
  so `pending.whenReady()` covers program compilation on both backends; frames behave as before (WebGL skips frames
  until the required programs are linked). `ImagePass.render` keeps rendering again while a render started work: it
  renders asynchronously (waits for pending work and the GPU between illumination iterations), which a synchronous
  collect mode can't run.
- Step 2 done: the skybox (cube textures on both backends) and the background image (`resources.loadImage`, the same
  decode on both backends; the asset, the upload and the mipmaps) are pending work with their `onload` as part of it, so
  that work waiting for them sees the background marked loaded. `ImagePass.render` renders again while a render started
  work, up to three more times (a finished load enables draws whose pipelines then compile);
  `ImagePass.updateBackground` updates and waits for the pending work. The visual harness no longer calls it: without
  the tracking, the background cases render blank on both backends.
- Step 1 done: `resources.pending` (`gl/pending-work.ts`) and `resources.loadImage` on both backends (WebGPU decodes
  with `createImageBitmap` without color space conversion), `loadImageTexture` removed, SMAA lookups loaded through it,
  `ImagePass.render` waits for `whenReady`, `Canvas3D` redraws on `settled`. The isolated SMAA visual case now matches
  full runs (SMAA applied on WebGL, WebGPU identical).
- Prompted by the SMAA lookup textures (area, search) that load asynchronously from images: renders before they load run
  SMAA with empty lookups (no effect), on both backends; isolated visual test runs of the SMAA case hit this. Related:
  [webgpu-viewer.md](webgpu-viewer.md) (milestone 3), [backend-selection.md](backend-selection.md),
  [gpu-timer.md](gpu-timer.md).

## Goal

One neutral way to know whether the graphics resources a frame needs are ready, for both backends:

- **Asynchronous loads**: textures from images (SMAA lookups, background image and skybox, user images), later other
  asset-backed resources.
- **Program compilation**: WebGL parallel shader compilation (today `ShaderManager` with
  `linkPrograms`/`finalizePrograms`); WebGPU shader modules and render pipelines (today compiled lazily at first use, no
  readiness).

Consumers wait or skip as they do for programs today: `Canvas3D.draw` skips frames until ready (non-blocking) and
redraws when something becomes ready; `ImagePass.render` (screenshots, visual tests) waits until ready, so the result
doesn't depend on timing.

## Current state

- **Image textures**: `loadImageTexture(src, cell, texture)` (`gl/texture.ts`) creates an `Image`, uploads in `onload`
  and updates the cell; it returns nothing, nothing tracks it, and no redraw is requested. Used by SMAA
  (`getWeightsRenderable`). Until loaded, WebGL samples an empty texture (zeros); WebGPU binds a zero placeholder (since
  milestone 3).
- **Background**: `BackgroundPass` loads the skybox (WebGL cube texture, `gl/webgl/texture.ts` loads six `Image`s) and
  the image variant with its own `loaded` flags, `isReady()` and an `onload` callback; `ImagePass.updateBackground()`
  returns a promise the visual harness awaits.
- **Programs, WebGL**: `resources.linkPrograms(variants)` starts linking pending programs of the variants (parallel with
  `KHR_parallel_shader_compile`), `finalizePrograms(variants, isSynchronous)` reports whether they are linked (or
  blocks). `ShaderManager` decides the required variants from props and scene (color, tracing, marking, emissive, depth;
  pick on demand); `Canvas3D.draw` skips the frame while `finalizeRequired` is false; `ImagePass.render` calls
  `ShaderManager.ensureRequired` (synchronous).
- **Programs, WebGPU**: `linkPrograms` does nothing, `finalizePrograms` returns true. Shader modules are created at a
  render item's first draw (`getWebGPUProgram`), pipelines at first use per render state, target formats and draw mode
  (`getPipeline`, synchronous `createRenderPipeline`). First use of a feature (e.g. the first hover: pick and marking
  pipelines) costs a frame hitch; compile errors surface asynchronously as uncaptured errors.
- **GPU completion**: `waitForGpuCommandsComplete` (WebGPU: `queue.onSubmittedWorkDone`) is for readback; uploads need
  no waiting (queue order).

## Design

### Pending-work tracker (`GraphicsResources`)

Neutral, on the graphics context:

```ts
interface PendingWork {
    /** Tracks `promise` until it settles; `label` for debugging, rejections are reported once */
    track<T>(promise: Promise<T>, label: string): Promise<T>;
    /** Whether nothing tracked is pending */
    isReady(): boolean;
    /** Resolves when nothing tracked is pending (including work tracked while waiting) */
    whenReady(): Promise<void>;
    /** Emits when tracked work settles, e.g. for a redraw */
    readonly settled: Subject<void>;
}
```

- `resources.pending: PendingWork` on both backends (shared implementation in `gl/`); counters in `stats` for debugging.
- Context loss and `destroy` clear the tracked set (the promises' results are ignored).

### Image textures

- `loadImageTexture` becomes `resources.loadImage(texture, source)`: decode, upload, update the cell, tracked as
  `image <label>`.
- Decoding per backend:
  - WebGL: `Image` + `onload` as today (upload with `UNPACK_COLORSPACE_CONVERSION_WEBGL = NONE`);
  - WebGPU: `createImageBitmap(blob or image, { colorSpaceConversion: 'none', premultiplyAlpha: 'none' })`, then
    `copyExternalImageToTexture`, so data textures (SMAA lookups) are exact by construction instead of relying on the
    decoded `<img>` having no color conversion.
- Background skybox and image go through the same tracker; `BackgroundPass.isReady` stays (it also covers "nothing to
  load"), and `ImagePass.updateBackground` becomes a special case of `whenReady`.

### Programs

`ShaderManager` stays the policy (which variants the scene and props need); the tracker is the mechanism.

- WebGL: `linkPrograms` tracks a promise per started link that resolves when `finalize` succeeds (polled from
  `finalizePrograms` as today, or with a timer); behavior unchanged.
- WebGPU:
  - `linkPrograms(variants)` creates the shader modules of pending render items for those variants (the program key is
    known: variant, defines, linear textures) and tracks `module.getCompilationInfo()`; compile errors are reported from
    the info (with the WGSL line), not only as uncaptured errors.
  - Pipelines: `createRenderPipelineAsync` for pipeline keys seen before (a per-program set of recently used keys) when
    the module is ready; a draw whose pipeline is still compiling is skipped in that frame (the frame is redrawn when
    the pipeline settles), falling back to synchronous creation when `isSynchronous` (screenshots, picking).
  - `finalizePrograms(variants, isSynchronous)` returns whether the tracked modules of the variants are ready;
    synchronous finalization creates what's missing synchronously.

### WebGPU pipelines (step 4, revised)

Done first (4a, the frame hitch): missing render pipelines compile in the background (`createRenderPipelineAsync`, on
Dawn's worker threads) and their draws are skipped until ready; the compilation is pending work, so `Canvas3D` redraws
when it settles. The first hover on the Viewer went from a 340–440 ms frame to at most 63 ms. Pipeline keys use
content-based render state ids (states built per draw, e.g. the solid interior stencil passes, used to create a pipeline
per draw). `ImagePass.render` renders a second time when the first one started work (`pending.trackedCount` changed),
also when that work settled before the render ended (e.g. between illumination iterations).

Done next (4b, key reuse): `linkPrograms` creates the WebGPU programs (shader modules) of the registered render items
for the variants (`ShaderManager` links the required ones at each commit), ahead of their first draws; every pipeline
configuration (render state, target formats, draw mode) drawn with is remembered per program family (program name and
render variant), and a new program of a family starts compiling all of them when it is created. Programs recreated after
device loss do the same. `ImagePass` waits for this pending work before rendering, so over the WebGPU visual profile its
second renders went from 42 to 26 of 90 cases; the rest use configurations their family has not drawn with yet. Key reuse
was removed later, when preparing covered what it did and it only added compiles
([pipeline-readiness.md](pipeline-readiness.md), 6g).

Open: `ShaderManager.ensureRequired` means "everything is ready" on WebGL (programs linked synchronously); on WebGPU it
covers the configurations seen before (key reuse), not new ones, and the second render remains for those. Pipelines are
keyed by render state, target formats and draw mode, which are only known when a pass draws, so full readiness needs the
passes' draw sequence:

- **Collect mode** (option B, planned): `resources.preparePrograms(run)` runs `run` (the draws of a pass sequence) with
  draws only requesting their pipelines (no commands recorded), then resolves when the requested modules and pipelines
  are compiled. `ShaderManager.ensureRequired` becomes async and, on WebGPU, prepares the `ImagePass` render for the
  required variants; `ImagePass.render` then renders once. WebGL resolves right after linking, unchanged.
- **Interactive warm-up** (picking done): after a scene change, `Canvas3D` links the pick and marking variants and
  renders the pick pass without reading it, in idle time (`ShaderManager.warmUp`, `PickHelper.warmUp`; attrib `warmUp`,
  off in the visual test harness, where it slowed the runs by a third). The first hover in the Viewer starts 8 pipeline
  compiles instead of 13. The rest are the marking pass (the `marking` variant, `edge`, `overlay`, `copy`), which draws
  only marked objects, so preparing it needs collect mode; likewise emissive and depth when enabled.
- **Shader compilation info**: `getCompilationInfo` of the modules created in `linkPrograms` tracked and reported with
  WGSL locations.

### Consumers

- `Canvas3D.draw`: skip the frame while `shaderManager.finalizeRequired()` is false (as today); image loads don't hold
  back frames (zero placeholders until loaded); request a draw on `pending.settled`.
- `ImagePass.render`: `await resources.pending.whenReady()` before rendering (after `ShaderManager.ensureRequired`),
  replacing the separate background wait.
- Visual harness: nothing case-specific; isolated runs of the SMAA case match full runs.
- Picking: unchanged (synchronous finalization of the pick variant), on WebGPU with synchronous pipeline creation for
  missing pick pipelines.

## Steps

1. **Tracker and image textures**:
   - `resources.pending` on both backends;
   - `loadImageTexture` tracked, with the redraw on settle;
   - `ImagePass.render` waits for `whenReady`;
   - WebGPU decodes with `createImageBitmap` (`colorSpaceConversion: 'none'`).

   Check: the SMAA visual case is identical in isolated and full runs on WebGL and WebGPU, and the Viewer's first frames
   apply SMAA once the lookups are loaded.

2. **Background through the tracker**: skybox and image loads tracked; `updateBackground` built on `whenReady`; visual
   background cases unchanged.
3. **WebGL programs through the tracker**: link promises; `ShaderManager` unchanged in policy; WebGL visual profiles
   byte-identical, Viewer behavior unchanged (frames skipped while compiling).
4. **WebGPU program readiness** (4a and 4b done, see WebGPU pipelines; the rest planned there):
   - modules created in `linkPrograms`, compilation info tracked and reported;
   - async pipelines for known keys, skipped draws while compiling, synchronous fallback;
   - `finalizePrograms` real.

   Check: no frame hitch when enabling marking or illumination-free features (timer groups around the first frames),
   compile errors with WGSL locations, visual suite unchanged.

5. **Optional**: asset-backed resources beyond images (e.g. volume data streamed into textures) tracked the same way.

## Validation

- WebGL and WebGPU visual profiles unchanged; the SMAA case identical in isolated and full runs.
- Viewer smoke test unchanged; manual check that toggling features doesn't stall frames on WebGPU (timer groups).
- Unit tests for the tracker (nested tracking while waiting, rejections, `settled` events) and for the WebGPU pipeline
  key reuse.

## Decisions

- While images load, `Canvas3D` draws with the zero placeholders and redraws when a load settles (an image appearing a
  moment later rather than a blank canvas); frames are still skipped while programs compile, as today.
  `ImagePass.render` always waits for `whenReady`.
- Image loading moves into `resources` (`resources.loadImage(texture, source)`); `loadImageTexture` is removed, the call
  sites already have the context.
