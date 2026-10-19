# WebGPU branch: follow-ups from the review

For Alex. A review of `gl-refactor-webgpu` (`1194e63af`, squashed; master merged on top), done 2026-10-18 by six
reviewers by area, the strongest findings checked again against the code by two more. Line numbers are as of
`1194e63af`. Status: **confirmed** (defect and scenario checked in the code), **plausible** (defect exists, scenario
unproven), **observed** (seen in the running viewer).

Fixed since the review, in the same branch:

- Graphics: correctness 1 to 16, memory 1, 2 and 3 (each marked below).

- Prism: `graphics-backend` option (`auto` by default), as the Viewer; canvas presence counts a released context at once
  (WebGPU reports the loss a task later) and never clears `released` on a refused second release; mounts are admitted
  into the context budget before the plugin makes its context (`admitCanvas`, async creation); long press follows a
  replaced canvas.
- `PluginContext`: a context made while the plugin was disposed is freed and nothing is built on it; overlapping mounts
  (React StrictMode) share one container and context.
- Docs: the `prefer-webgl1` removal in `breaking-v6-changes.md`, the `PickBuffers` contradiction there, the v5
  `instance.md` snippet restored, `graphics-backend` in the Prism viewer docs.

## Decision

- Decided: `PluginFeatureDetection.preferWebGl1` removed, `PluginConfig.General.PreferWebGl1` off by default (changelog
  and `breaking-v6-changes.md`).
## Correctness

1. Fixed: the Hi-Z reduction reads the GLSL block (`p + (-1,-1)`, `p + (0,-1)`); `renderer/hi-z/viewport` now matches
   WebGL.
2. Fixed: cylinder impostors return the far hit for `uDepthBack`, solid-interior back caps included, as GLSL.
3. Fixed: the direct-volume WBOIT weight is `alpha * alpha * w`, as GLSL.
4. Fixed: the compute isosurface wraps grid positions past the data (`uDataDim`, periodic volumes one point larger per
   axis), as the vertex estimate does; checked against the CPU marching cubes of the wrapped tensor.
5. Fixed: storage buffers of another compute context count as foreign (`hasForeignResources`, `compute.ownsBuffer`),
   so GPU surfaces in the vertex arena are computed again after a restore or switch and not transferred; a texture mesh
   reused across the two paths drops the other path's values (`isoChunks` and the arena, or its textures); the Gaussian
   density volume reuses only its own texture; reads failing with a released device are not reported. Visual case
   `renderer/context/release-restore-gpu-surface`; a backend switch is not covered by a visual case (the harness can't
   switch), its paths by unit tests.
6. Fixed: a pick whose render skipped draws (`GraphicsStats.skippedDraws`, pipelines compiling) is reported
   `'pending'` and drawn again once tracked work settles. A second cause on every backend: `PickHelper` built its pixel
   spiral before `pickPadding` was set, so a new helper (restore, backend switch, headless) found nothing until a resize.
7. Fixed: the read lag counts frames, not checks; a given-up read is started again and a click stays `'pending'`.
8. Fixed: `setBackend` and fallbacks check a `disposed` flag after the await and destroy the context they made.
9. Fixed: WebGPU loads fp16 textures (float data converted to half floats).
10. Fixed: WebGPU `recreate` asks for a device again while none is given (5 retries with backoff, about 3 s, as long as
    `Canvas3D` waits for a restore).
11. Fixed: `setBackend` while released keeps the backends for `restore`, which creates the context from them
    (`contextRestored`), or restores the current one if none can be created; no live context for a released canvas.
    Not covered by a test (the harness can't switch backends).
12. Fixed: prepare mode plans the background of the props (`BackgroundPass.prepare`, its variant's program) without
    updating the live one; images are neither loaded nor freed.
13. Fixed: `ImagePass` renders and reads again on the new context when it changed during a render or read
    (`setContext`), at most 3 times. Not covered by a test.
14. Fixed: a representation's updates run one after another (`Representation.createUpdateQueue` in every
    `createOrUpdate`), so `Canvas3D`'s update after a context change no longer interleaves with a state update.
15. Fixed: a restore keeps the canvas and needs no parent; a fallback (a new canvas in its place) still does.
16. Fixed: packed groups have the top byte 0x3f, so the f32 they are stored as is normal (never NaN, infinite or
    subnormal); the readback gives alpha 255 as group textures.

## Memory and performance

1. Fixed: `destroy()` unconfigures the canvas and destroys the device the context owns (`ownsDevice`, set by
   `createGraphicsContext` and the headless screenshot; not by tests passing their own device).
2. Fixed: WebGPU 2D and cube textures keep a `destroyed` flag; loads after `destroy` do nothing.
3. Fixed: the canvas is configured with `getPreferredCanvasFormat()` (`bgra8unorm` on Windows and macOS), the drawing
   buffer's target descriptor per context; drawing-buffer reads swap red and blue back.
4. **Per-unit Gaussian surfaces gated by the whole structure** (confirmed). `gaussian-surface-mesh.ts:112-131` checks the
   structure's box cubed against `maxStorageBufferSize / 4`, but each unit computes on its own box: assemblies above
   ~300 Å at resolution ≤ 1 fall back to the CPU on WebGPU (core limits; higher with the adapter's). The small-surface
   cutoff (cells and atoms) has the same flaw the other way: many small chains go to the GPU, each paying its fixed
   dispatch and readback cost. Gate on the unit's box and atoms. Deferred: the visual factory and `mustRecreate` see the
   structure only, a per-unit decision moves into `UnitsVisual`.
5. **Pass targets allocated up front** (measured; pre-existing on master). Deferred: part of rethinking the pass
   framework. Prism at 1872×1536 (DPR 2): ~836 MiB of textures on both backends (master WebGL 797 MiB). dpoit targets (132 MiB) exist with transparency `wboit`, tracing
   targets (167 MiB) with illumination off. Allocating passes on use would save ~300 MiB at this size. The branch adds
   ~40 MiB in the marking pass (93 → 137 MiB, target not identified) and, on WebGPU, ~58 MiB of depth-stencil
   (`depth24plus-stencil8`/`depth32float-stencil8` sizes, presumably).
6. WebGL, minor (from the review), done: the program cache memoizes the declaration hash per schema and defines, and
   WebGL1 items share one adapted schema per schema (`getWebGL1Schema`); emulated vertex ids past 2^24 (WebGL1) warn, a
   known limit (`minimum-capabilities.md`). Accepted as is: `setSolidInteriorPass` invalidating the globals per cap draw.
   Checked, fine: `histogram-pyramid/sum.ts` reads through `ctx.readPixels`, which ends the render pass.

## Checked and fine

No callers of removed APIs (`canvas3d.webgl`, `fromCanvas`) left in the repository; `identify` returning `'pending'`
handled where used; package dependencies and the build aliasing of backend sets; std140 layouts vs the WGSL structs, bind
group and pipeline cache keys; reversed-depth and top-left conventions; WebGL pass states equal to master. Rejected:
gating OIT on `EXT_float_blend` (correct per spec; the old check allowed invalid blending).
