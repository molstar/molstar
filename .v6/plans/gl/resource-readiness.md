# Plan: Resource readiness (async loads and program compilation)

## Status

- Step 1 done: `resources.pending` (`gl/pending-work.ts`) and `resources.loadImage` on both backends (WebGPU decodes
  with `createImageBitmap` without color space conversion), `loadImageTexture` removed, SMAA lookups loaded through it,
  `ImagePass.render` waits for `whenReady`, `Canvas3D` redraws on `settled`. The isolated SMAA visual case now matches
  full runs (SMAA applied on WebGL, WebGPU identical). Next: step 2 (background loads).
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

Open: `ShaderManager.ensureRequired` means "everything is ready" on WebGL (programs linked synchronously); on WebGPU it
covers neither shader modules nor pipelines, and the second render is a workaround. Pipelines are keyed by render state,
target formats and draw mode, which are only known when a pass draws, so readiness needs the passes' draw sequence:

- **Collect mode**: `resources.preparePrograms(run)` runs `run` (the draws of a pass sequence) with draws only
  requesting their pipelines (no commands recorded), then resolves when the requested modules and pipelines are
  compiled. `ShaderManager.ensureRequired` becomes async and, on WebGPU, prepares the `ImagePass` render for the
  required variants; `ImagePass.render` then renders once. WebGL resolves right after linking, unchanged.
- **Interactive warm-up**: after the scene settles, `Canvas3D` prepares the passes of variants not drawn yet but likely
  next (pick, marking; emissive, depth when enabled) in collect mode, in idle time, so the first hover or feature toggle
  draws completely in its first frame.
- **Key reuse**: pipeline keys seen per program name and variant are compiled for new programs of the same name and
  variant (new defines, new representations) as soon as they are created (`linkPrograms`), so new representations don't
  skip draws either.
- **Shader modules**: created in `linkPrograms` for the required variants, `getCompilationInfo` tracked and reported
  with WGSL locations.

Decisions (open): collect mode (exact, needs a dry run of the pass sequence) and key reuse (cheap, covers repeats only)
together or one of them; `ensureRequired` async (internal to `ImagePass`, `Canvas3D` keeps the non-blocking
`finalizeRequired`).

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
4. **WebGPU program readiness** (4a done, see WebGPU pipelines; the rest planned there):
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
