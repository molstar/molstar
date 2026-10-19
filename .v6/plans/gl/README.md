# Graphics plans

Plans of the graphics layer (`packages/graphics/src/gl`, `canvas3d`): the backend split, the WebGPU backend and the work
around it. Each plan's Status section has the details; this index lists where each stands. Rules for all of them:
[minimum-capabilities.md](minimum-capabilities.md) (targets and minimum specs) and
[neutral-conventions.md](neutral-conventions.md) (the neutral API follows WebGPU, WebGL adapts).

Open follow-ups from the review of the WebGPU branch (correctness, memory, one decision):
[../../webgpu/followup.md](../../webgpu/followup.md).

## Done

- [backend-split.md](backend-split.md): backend-neutral interfaces, so that WebGL and WebGPU backends sit side by side.
- [backend-narrowing.md](backend-narrowing.md): no WebGL specifics outside the backends; neutral readback, picking, Hi-Z
  and context loss.
- [uniform-buffers.md](uniform-buffers.md), [texture-bindings.md](texture-bindings.md),
  [attribute-bindings.md](attribute-bindings.md), [compute-shader-declarations.md](compute-shader-declarations.md):
  uniform, texture and attribute declarations and layouts from the schemas.
- [glsl-in-webgl.md](glsl-in-webgl.md): GLSL inside the WebGL backend, programs as descriptors (step 7 waits for build
  entries).
- [wgsl-strategy.md](wgsl-strategy.md): hand-written WGSL with a shared preprocessor; every program ported (step 8,
  reversed-Z, open).
- [webgpu-viewer.md](webgpu-viewer.md): the Viewer on WebGPU with every feature WebGL renders.
- [webgpu-compute-shaders.md](webgpu-compute-shaders.md), [webgpu-compute.md](webgpu-compute.md),
  [indirect-draws.md](indirect-draws.md): compute kernels and the GPU surface paths on WebGPU, at parity with WebGL or
  faster for Gaussian surface representations.
- [webgpu-device-limits.md](webgpu-device-limits.md): devices with the adapter's limits, a core-limits test profile.
- [webgpu-draw-overhead.md](webgpu-draw-overhead.md): per-draw CPU work of WebGPU render items reduced.
- [gpu-timer.md](gpu-timer.md): GPU timing per pass on both backends.
- [neutral-conventions.md](neutral-conventions.md): the neutral API follows WebGPU, WebGL adapts; reversed depth and the
  top-left rect origin in [depth-conventions.md](depth-conventions.md).
- [marching-cubes-benchmark.md](marching-cubes-benchmark.md): `pnpm benchmark`, including the render cases.
- [program-switch-frames.md](program-switch-frames.md): adaptive tracing compiles nothing (A); WebGPU frames held while
  a props change's pipelines compile, at most 1 s (B, through pipeline-readiness); new representations and define
  changes appear when their pipelines are ready (kept).
- [ssao-performance.md](ssao-performance.md): SSAO at half resolution by default with depth-aware upsampling, reading
  the full-resolution depth (occlusion independent of the resolution), 4x4 ordered kernel rotation; 32 samples kept.
- [gpu-surface-buffers.md](gpu-surface-buffers.md): WebGPU isosurfaces in 4 K vertex chunks of a vertex arena,
  allocated on the GPU without readback (atomic emission), scratch resources kept; multi-draw indirect (step 4)
  postponed.
- [resource-readiness.md](resource-readiness.md): pending work for images, backgrounds, WebGL links and WebGPU
  pipelines, collect mode and the marking warm-up; `ImagePass` rendering once continues in
  [pipeline-readiness.md](pipeline-readiness.md); optional step 5 (asset-backed resources beyond images) not done.
- [depth-conventions.md](depth-conventions.md): reversed-Z on both backends (WebGL with `EXT_clip_control` at WebGPU's
  precision, reversed without the gain otherwise) and the top-left neutral origin; the remaining SSAO banding (far from
  the camera) accepted ([ssao-blur-depth.md](ssao-blur-depth.md) dropped).
- [geometry-preparation.md](geometry-preparation.md): `GeometryUtils.prepare` (GPU surface vertices read back for the
  position iterator, for color and size themes) and one update path in the visuals, replacing
  `TextureMesh.withVertexData`.
- [backend-switch-followups.md](backend-switch-followups.md): geometry that depends on capabilities, the tunnels'
  context in params, fallback to the next backend when a context is lost for good (without loops), and restore as a
  rebuild on a new context (`GraphicsContext.recreate`, in-place resets removed).

## In progress or paused

- [backend-selection.md](backend-selection.md): the Viewer defaults to `auto` (WebGPU where capable), the plugin default
  follows after a release; backend sets, the bundle check, Viewer variants with one backend and runtime switching
  keeping the scene done.
- [webgpu-draw-submission.md](webgpu-draw-submission.md): paused after measuring the upper bound of render bundles (CPU
  time of `render/boxes-1000` from 5.5 to 0.1 ms per frame).
- [pipeline-readiness.md](pipeline-readiness.md): `prepare` from params without rendering (pipeline requests, renderer
  collect, passes as plans, `Canvas3D.prepare`, frames held after props changes, prepared compute kernels) done;
  snapshot replay without compile hitches (step 6): precompile lists (object signatures, recorded and prepared) and
  their derivation from a state's params, prepared when a state or story is loaded, done and measured (WebGPU stories play
  without compiling); fewer variants per signature done (kinase story: 24 signatures, 1.6 s on WebGPU); WebGL's first draws into several
  targets (an ANGLE compile once per browser) accepted.
- [fewer-defines.md](fewer-defines.md): paused, good enough for now: the color, size and instance-granularity types
  as runtime selectors (kinase story 24 to 18 precompiled signatures); the clip object count and shading switches
  open.
- [sphere-impostors.md](sphere-impostors.md): analysed, paused: sphere impostors are fragment-bound and get no early-Z
  (depth writes); conservative depth without gain on NVIDIA; WebGPU's WGSL shading costs 3.8 ms more at 2560x1440
  (generated code).

## Dropped

- [ssao-blur-depth.md](ssao-blur-depth.md): linear depth in the SSAO target is good enough; reading the scene depth in
  the blur (slower) only for banding beyond about 400 Å from the camera.
- [sphere-depth-prepass.md](sphere-depth-prepass.md): a depth pre-pass for sphere impostors; dropped after measuring,
  the depth-only pass costs about as much as the impostors without lighting; its read-only depth passes kept.
- [webgpu-sync-readback.md](webgpu-sync-readback.md): synchronous picking reads on WebGPU through a canvas copy; worked
  in Chrome and Firefox, dropped because clicks already wait for the asynchronous read (API callers await
  `Canvas3D.asyncIdentify(…).wait()`).

## Not started

- [webgpu-xr.md](webgpu-xr.md): XR on WebGPU through `XRGPUBinding` behind the neutral `GraphicsXR`.
- [storage-buffers.md](storage-buffers.md): data arrays in storage buffers on WebGPU (layer B of draw submission).
- [primitive-expansion.md](primitive-expansion.md): cylinders, lines and text expanded in the renderer.
- [core-gpu-compute.md](core-gpu-compute.md): a neutral map-kernel tier and GPU algorithms in core.
- [attribute-color-size.md](attribute-color-size.md): per-vertex color and size attributes, to be considered.
