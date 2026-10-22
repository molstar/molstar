# Plan: Minimize narrowing to a graphics backend

## Status

- **Current**: done, steps 1–5; decisions 1A, 2A and 3 taken.
- Step 1 done: the unconditional `gl.finish()` calls are gone (`active-voxels.ts` is neutral; the rest of the WebGL
  marching cubes follows the histogram pyramid's sum readback, step 3); `GaussianDensityGPU` and `fieldFromTexture2d`
  removed (recorded as breaking); loop bounds are `UniformSpec.loopBound`, turned into defines by the WebGL render item
  on WebGL1 (the direct-volume renderable no longer narrows). Dropping `gl.finish()` gained no speed: end to end (the
  representation benchmarks) the times are unchanged within the run-to-run variation, per-stage benchmark times of small
  grids are about 1.5× slower (their own sync wait after each stage replaces the finish).
- Step 2 done: `GraphicsContext.readTextureData` reads whole 2d color textures on both backends (data order, four values
  per texel; ubyte, float, half float); geo-export reads its grids and texture meshes with it (indirectly drawn surfaces
  with their read-back count), the texture mesh position iterator uses vertices read back before
  (`TextureMesh.withVertexData` in the visuals, when the color theme is per vertex), the sync helpers of
  `gl/compute/util` and `TextureMesh.meta.webgl` are gone. Check cases `gpu-paths/read-texture-data`,
  `export-gpu-surface` and `vertex-colors-gpu-surface` (both fail with the read back removed).
- Step 3 done: the fragment-shader marching cubes (active voxels, histogram pyramid, isosurface buffers) and the WebGL
  half of grid3d live in `gl/webgl/compute/` and register with each WebGL context (`registerAtlasIsosurfaceExtraction`,
  `registerGrid3dFragmentImplementation`, explicit calls in `createContext`); the volume isosurface, segment and
  Gaussian surface visuals choose compute kernels or the registered atlas extraction by capability, without backend
  checks; grid3d's WGSL travels with its `KernelCode` (`sources.wgsl`) instead of the WebGPU kernel registry. Outside
  the backends only the passes (picking, Hi-Z, step 4), `canvas3d.ts` and headless creation (step 5) narrow.
- Step 4 done. `graphics.createTextureReader({ format })` returns a `TextureReader` (`gl/texture-reader.ts`) for
  repeated reads of the same rectangle of textures, `rgba8` (bytes) or `r32f` (floats), rows bottom to top. WebGL reads
  with `readPixels` and, on WebGL2, asynchronously with pixel pack buffers and a polled fence
  (`gl/webgl/texture-reader.ts`, its own framebuffer; `r32f` reads rgba/float and keeps red where red/float isn't
  readable, which Hi-Z used to give up on). WebGPU reads asynchronously with staging buffers
  (`gl/webgpu/texture-reader.ts`, rows flipped). Readers report pending/resolved/failed only; when to give up a pending
  read is the caller's (picking: `maxAsyncReadLag` checks and a frame, now on WebGPU too; Hi-Z: `maxFrameLag` frames).
  The design's `createReadbackBuffer(byteLength)` became a reader of textures (the pick targets and the Hi-Z levels are
  textures on all backends).
  - Picking: `PickBuffers` is one neutral class over a reader, `checkAsyncPickingSupport` takes the buffers, the
    `PickPass.bind*` methods are gone. Viewer smoke: 8/25 hover hits on WebGL and WebGPU, as before.
  - Hi-Z: neutral, and on WebGPU (`hi-z.wgsl.ts`, `textureLoad` of the previous level or the depth texture, the level
    size in `uTexSize`); it needs float render targets and asynchronous reads (WebGL2, WebGPU). It starts no read of
    levels with skipped draws (work pending, e.g. a WebGPU pipeline compiling), which read back as depth 0 and occluded
    everything. Check case `gpu-paths/hi-z-occlusion` (a near quadrant of a depth target, spheres behind and beside it;
    fails with the WGSL row flip or the pending-work guard removed).
  - The guard test's file allowlist is empty: outside `gl/webgl/` no file uses WebGL specifics.
- Step 5 done. `Canvas3D.webgl` and `Canvas3DContext.webgl` are removed (decision 2). Context loss is neutral: the WebGL
  context handles the loss and restore events of its canvas itself and implements `restore` (resolved by the restore
  event), so `Canvas3DContext` recovers both backends with `graphics.restore(resetPasses)`; the debug shortcut
  (shift+ctrl+alt and a click) calls `graphics.simulateContextLoss()` (WebGL: `WEBGL_lose_context`, restored after a
  second; WebGPU: `device.destroy()`), on WebGPU too now; `destroy` takes `GraphicsDestroyOptions` on all backends. The
  Viewer smoke test loses and restores the WebGL context as well (fails with the restore handler removed). Guard test:
  no `isWebGLContext`/`isWebGPUContext` calls or casts outside the backends (none left, no allowlist), and the import
  check covers `gl/webgpu/` too; backend imports outside the backends are only context creation (`fromCanvas`, headless)
  and backend descriptors. Examples and tests narrow where they are WebGL-specific.
- Decisions taken (1A, 2A, 3). Follows [backend-selection.md](backend-selection.md) (import hygiene, GPU compute as
  backend services) and [backend-split.md](backend-split.md) step 5 (readback and sync); the guard test is
  `gl/_test/backend-neutral.test.ts`.

## Principle

Neutral code decides by **capability or service**, not by backend type: `graphics.compute` (kernels available),
`capabilities.*`, neutral services that pick their implementation inside. Narrowing to a backend type (`isWebGLContext`,
`isWebGPUContext`, `WebGLContext`/`WebGPUContext` types, imports of `gl/webgl/` or `gl/webgpu/`) is fine only for:

1. code inside `gl/webgl/` and `gl/webgpu/`;
2. creating and choosing backends (descriptors, headless context creation, apps choosing their backend list);
3. tests that compare backends (benchmark paths, visual profiles).

Capability checks are not narrowing: they keep working when a backend gains a feature (e.g. compute on WebGL through the
map-kernel tier of [core-gpu-compute.md](core-gpu-compute.md)).

## Inventory

Outside `gl/webgl/` and `gl/webgpu/`, excluding tests:

| Site                                                                                                                | Narrows for                                                                    | Category                                 |
| ------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------ | ---------------------------------------- |
| `apps/viewer/src/plugin-spec.ts`, `plugin/core/src/context.ts`                                                      | backend descriptors (default list)                                             | A: construction                          |
| `plugin/headless/src/screenshot.ts`                                                                                 | headless context creation (`gl` package, Dawn)                                 | A: construction; one compat accessor (F) |
| `canvas3d/canvas3d.ts`                                                                                              | `fromCanvas` (sync WebGL creation)                                             | A: construction                          |
| `canvas3d/canvas3d.ts`                                                                                              | `webgl` accessor, debug context loss (`WEBGL_lose_context`), `destroy` options | F: compat and debug                      |
| `repr/structure/visual/gaussian-surface-mesh.ts` (6), `repr/volume/isosurface.ts` (3), `repr/volume/segment.ts` (3) | WebGL histogram pyramid or compute kernels                                     | B: implementation chosen in visuals      |
| `gl/compute/histogram-pyramid/reduction.ts`, `marching-cubes/active-voxels.ts`, `marching-cubes/isosurface.ts`      | an unconditional `gl.finish()` only                                            | C1: needless                             |
| `gl/compute/histogram-pyramid/sum.ts`                                                                               | synchronous one-pixel `readPixels`                                             | C2: WebGL implementation                 |
| `gl/compute/util.ts` (`readTexture`, `readSingleChannelTexture`)                                                    | synchronous float readback (`asGLTexture`)                                     | C3: missing neutral readback             |
| `gl/compute/grid3d.ts`                                                                                              | GLSL fragment path or WGSL kernel                                              | C4: two implementations in one module    |
| `gl/renderable/direct-volume.ts`                                                                                    | WebGL1 loop bound as a define                                                  | D: WebGL quirk in neutral code           |
| `canvas3d/passes/pick.ts`                                                                                           | pixel pack buffers and fences (WebGL2), sync reads (WebGL1)                    | E: missing neutral async readback        |
| `canvas3d/passes/hi-z.ts`                                                                                           | fences and pixel pack buffers; unsupported elsewhere                           | E: same; Hi-Z is off on WebGPU           |
| `geo/geometry/texture-mesh/texture-mesh.ts`                                                                         | `meta.webgl` for the position iterator's `readPixels`                          | G: missing neutral float readback        |
| `geo/gaussian-density/gpu.ts`                                                                                       | `GaussianDensityGPU` and `fieldFromTexture2d` (sync readback)                  | H: no callers                            |
| `packages/extensions/geo-export/src/mesh-exporter.ts`                                                                        | texture readback of GPU surfaces                                               | C3 (also the default-backend blocker)    |

## Findings

- **B**: three visuals carry twelve backend checks to choose between two implementations of the same algorithm (and
  repeat the choice in `createGeometry` and `mustRecreate`). The choice belongs behind one neutral entry point.
- **C1**: three WebGL marching-cubes stages are typed `WebGLContext` only to call `gl.finish()`, unconditionally; the
  histogram sum's `readPixels` already synchronizes (removing them turned out neutral for speed, see Status).
- **C3, G**: the neutral `readTextureAsync` reads rgba8 only; float and single-channel reads exist on WebGPU
  (`readTexture`) and as WebGL-only sync helpers, so GPU surfaces can't be exported or iterated neutrally (the two
  blockers of the default-backend recommendation).
- **D**: the WebGL1 loop-bound workaround in a neutral renderable contradicts "WebGL adapts, the neutral API follows
  WebGPU"; grid3d's `loopBounds` is the same quirk at the spec level.
- **E**: picking and Hi-Z implement async readback themselves with WebGL objects; a neutral async readback buffer
  (WebGL2: pixel pack buffer and fence; WebGL1: sync; WebGPU: staging buffer) would make both neutral and turn Hi-Z on
  for WebGPU.
- **H**: dead code with WebGL types, exported (removing it is a breaking change to record).

## Design

1. **Neutral GPU surface entry points**: `computeIsosurfaceTextureMesh(graphics, input, …)` (volumes, segments) and the
   Gaussian surface counterpart in `gl/compute` (or `geo/`), with `isGpuIsosurfaceSupported(graphics)`; they choose
   compute kernels when `graphics.compute` has them, else the fragment-shader pipeline (WebGL). Visuals call one
   function and keep their fallback reasons.
2. **The fragment-shader pipeline as WebGL code**: histogram pyramid, active voxels, isosurface buffers and the GLSL
   half of grid3d move to `gl/webgl/compute/`, reached only through the entry points; `gl/compute` keeps neutral entry
   points and kernels (decision 1).
3. **Neutral readback**: `readTextureAsync` for all color formats of 2d textures (rgba8, r8, float, half float), as
   WebGPU's `readTexture` does; `gl/compute/util` sync helpers go; geo-export and the texture mesh position iterator
   read asynchronously (`meta.graphics` instead of `meta.webgl`).
4. **Neutral async readback buffer** for repeated reads (picking, Hi-Z): `graphics.createReadbackBuffer(byteLength)`
   with `read(target, rect)` and a fence-like `isReady`/promise; WebGL2 pixel pack buffers and fences, WebGL1
   synchronous, WebGPU staging buffers.
5. **WebGL quirks in the backend**: the WebGL render item turns uniforms marked as loop bounds into defines on WebGL1
   (`UniformSpec` option, or the program description), so neutral schemas don't branch.
6. **Compat accessors**: `Canvas3D.webgl` and `Canvas3DContext.webgl` (and headless) removed in favor of `graphics`
   (decision 2, a breaking change); debug context loss through a neutral `graphics.simulateLoss?()` (WebGL:
   `WEBGL_lose_context`, WebGPU: `device.destroy()`, as the smoke test does).
7. **Guard test**: besides imports, `isWebGLContext`/`isWebGPUContext` outside the backends need an allowlist entry; the
   allowlist shrinks to category A.

## Decisions

1. **Fragment-shader compute**: A. move it to `gl/webgl/compute/` (honest: it is the WebGL implementation, and the
   neutral entry points hide it); B. keep it in `gl/compute`, typed `GraphicsContext` where possible. Chosen: A.
2. **`webgl` accessors**: A. remove (v6 breaking change, `graphics` exists); B. deprecate for a release. Chosen: A,
   recorded in the breaking changes.
3. **`GaussianDensityGPU`**: removed (with `fieldFromTexture2d`), recorded as breaking.

## Steps

1. Quick wins: drop the unconditional `gl.finish()` calls (C1; benchmark the WebGL isosurface before and after), remove
   `GaussianDensityGPU` (H), loop bounds in the WebGL backend (D, 5).
2. Neutral readback of all texture formats (3): geo-export and the position iterator on both backends, the sync helpers
   removed; closes the two blockers of the default backend.
3. Neutral GPU surface entry points (1) and the fragment-shader pipeline behind them (2, decision 1); grid3d split.
4. Neutral async readback buffer (4): picking and Hi-Z neutral, Hi-Z on WebGPU (a visual case).
5. Compat accessors and debug context loss (6); the guard test for narrowing (7), allowlist down to construction.
