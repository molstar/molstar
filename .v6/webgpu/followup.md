# WebGPU branch: follow-ups from the review

For Alex. A review of `gl-refactor-webgpu` (`1194e63af`, squashed; master merged on top), done 2026-10-18 by six
reviewers by area, the strongest findings checked again against the code by two more. Line numbers are as of
`1194e63af`. Status: **confirmed** (defect and scenario checked in the code), **plausible** (defect exists, scenario
unproven), **observed** (seen in the running viewer).

Nothing under `packages/graphics` was changed after the review (the GL freeze while the backend lands). Fixed outside it,
in the same branch:

- Prism: `graphics-backend` option (`auto` by default), as the Viewer; canvas presence counts a released context at once
  (WebGPU reports the loss a task later) and never clears `released` on a refused second release; mounts are admitted
  into the context budget before the plugin makes its context (`admitCanvas`, async creation); long press follows a
  replaced canvas.
- `PluginContext`: a context made while the plugin was disposed is freed and nothing is built on it; overlapping mounts
  (React StrictMode) share one container and context.
- Docs: the `prefer-webgl1` removal in `breaking-v6-changes.md`, the `PickBuffers` contradiction there, the v5
  `instance.md` snippet restored, `graphics-backend` in the Prism viewer docs.

## Decision

- **`auto` ignores `PluginFeatureDetection.preferWebGl1`** (`plugin/core/src/util/graphics-backends.ts:26`). That
  detection is true on all iOS and on Safari 15.1–15.3, where the Viewer used WebGL1 before; now they get WebGPU or
  WebGL2. The comment says it is on purpose; the changelog does not say so, and its "`webgl`: WebGL2, else WebGL1" is
  wrong on iOS (WebGL1 only). Keep it, or pass `preferWebGl1` to the `all` set in `auto` (the conservative choice).

## Correctness

1. **Hi-Z reduction reads a misaligned block in y** (confirmed). `gl/webgpu/shader/hi-z.wgsl.ts:27-29` reads rows `p`,
   `p+1`; GLSL reads `p`, `p-1`. With power-of-two levels `p.y = 2j+1`, so WGSL reads {2j+1, 2j+2}; the shift compounds
   to about a texel per level and `HiZPass.isOccluded` can cull visible objects on WebGPU. Fix: offsets `(-1,-1)` and
   `(0,-1)`. One line; check with the Hi-Z visual cases.
2. **Cylinder impostors ignore `uDepthBack`** (confirmed). `cylinders.wgsl.ts:205-345` always returns the front hit;
   GLSL `cylinders.frag.ts:144-150` (and the solid-interior back caps) return the far one, as `spheres.wgsl.ts:243` does.
   The thickness pass of illumination (`thicknessMode: 'auto'`, the default) gets ~0 thickness for bonds: light and
   shadow leaks with ball-and-stick on WebGPU. Port of the GLSL branch; the illumination visual cases should cover it.
3. **Direct-volume WBOIT weight** (confirmed, live). `chunks/output.wgsl.ts` `wboit_output_wgsl` writes
   `select(1, alpha, uTransparentBackground) * alpha * w`; GLSL `wboit-write.glsl.ts:24-28` writes `alpha * alpha * w`
   for `dGeometryType_directVolume` (premultiplied color). Volumes do go through the WBOIT targets in the default mode
   (`draw.ts` `_renderWboit`), so with an opaque background WebGPU volumes are weighted by an extra 1/alpha. One line.
4. **Periodic volumes on the compute isosurface** (confirmed). `repr/volume/isosurface.ts:331-336` uploads the raw `n³`
   data with `gridDim = n+1` and linear strides; `marching-cubes.wgsl.ts:21-24` only clamps. At `x = n` it reads the next
   row, at `z = n` past the buffer: wrong triangles on the wrapped faces of full-cell CCP4 maps (`wrap: 'auto'`). WebGL
   builds the atlas from the wrapped `dataOffset`. Fix: modulo in `value()` or upload wrapped data; needs a periodic test
   case.
5. **GPU surfaces across a context change** (confirmed). Larger: a rule for rebuilding compute meshes.
   - Chunked texture meshes (`texture-mesh.ts:152-246`) hold null textures, so `hasForeignTextures` is false, and
     TextureMesh has no `isCompatible`: after WebGPU→WebGL (when WebGL can still do GPU surfaces, `mustRecreate` false)
     the object is re-added as is and draws nothing; after a WebGPU restore it binds the dead device's buffers.
   - `createChunked(..., textureMesh)` reusing a mesh last filled on WebGL keeps its WebGL textures: `hasForeignTextures`
     stays true, so `canvas3d.ts:1871-1879` never re-adds it (invisible for good). `TextureMesh.create(..., old)` never
     resets `isoChunks`/`isoTable` either.
   - `gaussian-density-volume.ts:50, :144` passes `gridTexture.ref.value` without the `ownsTexture` check the other paths
     have: throws on WebGL→WebGPU (WebGPU kernels on a WebGL texture), empty and still foreign after a WebGL restore.
6. **First pick after a restore or load finds nothing** (confirmed, observed). Pick pipelines compile asynchronously and
   their draws are skipped, but `PickHelper.plan` clears `dirty` (`pick-helper.ts:116`), so `asyncIdentify().wait()`
   reads cleared buffers and resolves `undefined` instead of waiting (Prism: first pick right after `restore()`; works
   200 ms later). Fix: keep `dirty` and report `'pending'` while `graphics.resources.pending` is not ready, or await
   `prepare({ interaction: true })` first. Matters for code that picks once (agent tools, tests).
7. **A slow pick read clears the selection** (confirmed). `PickBuffers.check` (`passes/pick.ts:380-392`) fails after 5
   calls (≥16 ms), counted per call: a pending click and the hover pick poll the same buffers each tick, so ~3 frames.
   Failed falls through to `undefined` in `asyncIdentify().tryGet` (`pick-helper.ts:227-229`), and WebGPU clicks are
   always async: under GPU load (illumination accumulating) a click emits empty loci and `clickDeselectAllOnEmpty`
   clears the selection. Count frames, or report `'pending'` on Failed for clicks.
8. **Dispose during `setBackend` or a fallback** (confirmed). `canvas3d.ts:540-595`: no disposed check after
   `await createGraphicsContext`; the canvas is swapped, a new context and Passes are made and leak, and
   `previous.destroy()` runs again on the destroyed context. A `disposed` flag checked after the await.
9. **fp16 direct volumes throw on WebGPU** (confirmed, non-default). `dataType: 'halfFloat'`: WebGPU reports
   `halfFloatTexture: true` (`webgpu/context.ts:978`) but `webgpu/texture.ts:262` `load` is `notSupported` for fp16.
   Either load fp16 or fall back to float/byte on WebGPU.
10. **Restore asks for a device once** (plausible). `webgpu/context.ts:1295-1301` with `canvas3d.ts:505-513`: no device
    right after a real loss (GPU process restarting) → `failed('not restored')` at once, WebGPU never tried again on the
    page; WebGL waits out `ContextRestoreTimeoutMs`.
11. **`setBackend` while released** (plausible). `canvas3d.ts:576` has no `released` guard: a live context, `released`
    reset, `contextChanged` but no `contextRestored`; Prism's canvas presence keeps a stale `lostSeen`. Only through
    `setGraphicsBackends` on a released viewer.
12. **Prepare mode changes the live background** (from the review, not re-checked). `precompile`/`prepare({ props })`
    reach `BackgroundPass.update` with the hypothetical props (`postprocessing.ts:500`, `illumination.ts:593`):
    `clearSkybox`/`clearImage` destroy the current texture and the other snapshot's image loads.
13. **A screenshot across a context switch** (from the review). `passes/image.ts:154, :212-227` captures the old renderer
    and helper before the awaited illumination loop; `setContext` swaps the passes underneath.
14. **`setContext` runs `repr.createOrUpdate()` outside the state system** (from the review, `canvas3d.ts:2072-2081`):
    may interleave with a running state update (trajectory playback).
15. **Restore needs `canvas.parentNode`** (from the review, `canvas3d.ts:541`): a restore keeps the canvas, but a canvas
    without a parent can't be restored (nor fall back).
16. **Group ids as NaN bit patterns** (plausible, low risk). `marching-cubes-chunks.wgsl.ts:147` stores
    `bitcast<f32>(pack4x8unorm(vec4(rgb, 1.0)))` (top byte 0xFF: NaN/-Inf for about half the groups); only moved and
    bitcast back today, but WGSL lets implementations assume no NaN. A u32 lane removes it.

## Memory and performance

1. **`destroy()` keeps the device** (confirmed). `webgpu/context.ts:1303-1317` frees the ring, timer and drawing buffer
   only; never `device.destroy()` nor `unconfigure()`. Every switch away from WebGPU, each disposed viewer and the
   superseded context in `recreateLost` (`canvas3d.ts:507`) keep a whole device until GC. Fix: destroy the device the
   context owns (not one passed in, e.g. headless tests). The biggest memory win.
2. **A destroyed texture reallocates on load** (confirmed). `webgpu/texture.ts:211-227, :357-365`: no destroyed flag, so
   `loadImage` finishing after `destroy` (background changed while decoding) creates a texture never freed, plus mips
   (`background.ts:538`); byte counts stay inflated. WebGL keeps a `destroyed` flag.
3. **Canvas format** (confirmed, observed: Chrome's "different format than is preferred" warning). The canvas is
   configured `rgba8unorm` (`webgpu/context.ts:77-86, :295-300`); Mac and Windows prefer `bgra8unorm`, one extra copy per
   frame. Use `getPreferredCanvasFormat()` per context; make `DrawingBufferTarget` per context (pipeline keys, prepared
   pipelines); swap R/B in `readDrawingBuffer` and in `readRenderTarget` of draw targets.
4. **Per-unit Gaussian surfaces gated by the whole structure** (confirmed). `gaussian-surface-mesh.ts:112-131` checks the
   structure's box cubed against `maxStorageBufferSize / 4`, but each unit computes on its own box: assemblies above
   ~300 Å at resolution ≤ 1 fall back to the CPU on WebGPU. Gate on the unit's box.
5. **Pass targets allocated up front** (measured; pre-existing on master). Prism at 1872×1536 (DPR 2): ~836 MiB of
   textures on both backends (master WebGL 797 MiB). dpoit targets (132 MiB) exist with transparency `wboit`, tracing
   targets (167 MiB) with illumination off. Allocating passes on use would save ~300 MiB at this size. The branch adds
   ~40 MiB in the marking pass (93 → 137 MiB, target not identified) and, on WebGPU, ~58 MiB of depth-stencil
   (`depth24plus-stencil8`/`depth32float-stencil8` sizes, presumably).
6. WebGL, minor (from the review): points are 6-vertex quads now, so where the vertex id is emulated (WebGL1) more than
   ~2.8M points exceed float precision of `aVertex`; `setSolidInteriorPass` (`renderer.ts:839`) invalidates all globals
   per cap draw; the program cache hashes the full declarations per lookup (`webgl/resources.ts:181-199`), and
   `getWebGL1Schema` copies the schema per direct-volume item; `histogram-pyramid/sum.ts:84` binds the drawing buffer
   without `endRenderPass`.

## Checked and fine

No callers of removed APIs (`canvas3d.webgl`, `fromCanvas`) left in the repository; `identify` returning `'pending'`
handled where used; package dependencies and the build aliasing of backend sets; std140 layouts vs the WGSL structs, bind
group and pipeline cache keys; reversed-depth and top-left conventions; WebGL pass states equal to master. Rejected:
gating OIT on `EXT_float_blend` (correct per spec; the old check allowed invalid blending).
