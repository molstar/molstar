# Graphics plans

Plans of the graphics layer (`packages/graphics/src/gl`, `canvas3d`): the backend split, the WebGPU backend and the work
around it. Each plan's Status section has the details; this index lists where each stands. Rules for all of them:
[minimum-capabilities.md](minimum-capabilities.md) (targets and minimum specs) and
[neutral-conventions.md](neutral-conventions.md) (the neutral API follows WebGPU, WebGL adapts).

Follow-ups from the review of the WebGPU branch (correctness, memory, one decision), all fixed or decided:
[../../webgpu/followup.md](../../webgpu/followup.md).

## Done

- [backend-split.md](backend-split.md): backend-neutral interfaces, so that WebGL and WebGPU backends sit side by side.
- [backend-narrowing.md](backend-narrowing.md): no WebGL specifics outside the backends; neutral readback, picking, Hi-Z
  and context loss.
- [uniform-buffers.md](uniform-buffers.md), [texture-bindings.md](texture-bindings.md),
  [attribute-bindings.md](attribute-bindings.md), [compute-shader-declarations.md](compute-shader-declarations.md):
  uniform, texture and attribute declarations and layouts from the schemas.
- [glsl-in-webgl.md](glsl-in-webgl.md): GLSL inside the WebGL backend, programs as descriptors, the bundle check
  (`pnpm check:backend-bundles`).
- [wgsl-strategy.md](wgsl-strategy.md): hand-written WGSL with a shared preprocessor; every program ported, reversed-Z
  (step 8) with [depth-conventions.md](depth-conventions.md).
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
- [pass-targets.md](pass-targets.md): pass targets allocated on use (review item M5): lazy storage, released
  when the props and scene don't draw into them, fewer and smaller targets (a default frame 77 MiB opaque, 132 MiB
  transparent instead of 932 MiB at 1872×1536); sharing intermediates between passes measured, not needed for now.

- [shared-textures.md](shared-textures.md): render items share one texture per image object for specs that opt in
  (`tFont`): one GPU copy of a font atlas per context, uploaded once per change.
- [screen-space-lod.md](screen-space-lod.md): sphere levels of detail by their size on screen (`lodPixelSize`, CSS
  pixels, coarser strides 1, 10, 100, 1000), the mesoscale preset's and the particle representations' default.

## In progress or paused

- [backend-selection.md](backend-selection.md): the Viewer and the plugin default to `auto` (WebGPU where capable);
  backend sets, the bundle check, Viewer variants with one backend and runtime switching
  keeping the scene done.
- [webgpu-draw-submission.md](webgpu-draw-submission.md): render bundles as draw lists for every pass but illumination's
  inputs, data changes replayed, layers A and caps done; changing chunks drawn directly, `drawLists` on by default
  (no slower than without in the measured cases); step 8 open
  (upper bound: `render/boxes-1000` CPU 5.5 to 0.1 ms per frame).
- [draw-table.md](draw-table.md): objects as flat rows (scene facts, renderer pass masks and per-view loops) replacing
  per-object planning, done (steps 1-5); its GPU side superseded by GPU culling per merged object (its unused row
  buffer removed).
- [shared-draw-resources.md](shared-draw-resources.md): geometry sub-allocated in shared buffers, theme layers in
  storage pages updated by range, shared object groups done (WebGPU replay of 10k draws 37 to about 1 ms per frame
  measured); the CPU geometry store deferred (measured: little to gain).
- [merged-render-objects.md](merged-render-objects.md): CPU-side merging (`renderMerged`) the main path for many units
  of a representation: `auto` from 30 units, update costs, geometry and layers without copies, the costs per member, GPU
  culling with occlusion (`gpuCulling: 'auto'`, from 20k instances) done; left: the visuals' own update work (the GPU
  culling's instance copy on weak GPUs measured, about 2% to gain, not taken up).
- [shader-f16.md](shader-f16.md): paused after measuring: the impostor lighting in half precision 6-11% faster on an
  integrated AMD GPU, 4-12% slower on NVIDIA (the conversions); per-GPU choice not worth it now.
- [pipeline-readiness.md](pipeline-readiness.md): `prepare` from params without rendering (pipeline requests, renderer
  collect, passes as plans, `Canvas3D.prepare`, frames held after props changes, prepared compute kernels) done;
  snapshot replay without compile hitches (step 6): precompile lists (object signatures, recorded and prepared) and
  their derivation from a state's params, prepared when a state or story is loaded, done and measured (WebGPU stories play
  without compiling); fewer variants per signature done (kinase story: 24 signatures, 1.6 s on WebGPU); WebGL's first draws into several
  targets (an ANGLE compile once per browser) accepted.
- [fewer-defines.md](fewer-defines.md): paused, good enough for now: the color, size and instance-granularity types
  as runtime selectors (kinase story 24 to 18 precompiled signatures), the clip object count on WebGL (WebGPU keeps it
  a pipeline constant); the shading switches measured, not converted (a frame cost for programs rarely saved); the
  layers' presence open.
- [storage-buffers.md](storage-buffers.md): done, in part by other work: the theme layers read from shared storage pages
  on WebGPU ([shared-draw-resources.md](shared-draw-resources.md)), an item's layers beyond the binding limit across up
  to 4 pages; the layers `DataSpec`s with accessors in GLSL and WGSL, WebGL laying them out.
- [sphere-impostors.md](sphere-impostors.md): indexed quads done (18-28% with many small spheres, one shared index
  buffer per context); before: analysed, paused: sphere impostors are fragment-bound and get no early-Z
  (depth writes); conservative depth without gain on NVIDIA; WebGPU's WGSL shading costs 3.8 ms more at 2560x1440
  (generated code).

- [direct-volume-performance.md](direct-volume-performance.md): the ray marcher only inside the volume's box (step 1,
  about 2x for maps, 4x for tilted slabs, truncated rays fixed); next: skipping transparent samples, jitter,
  empty-space skipping, reduced resolution.

- [gizmo-layer.md](gizmo-layer.md): proposed: overlay layers replacing `HandleHelper` and `CameraHelper` (gizmos above
  the scene, scene depth as occluded alpha, screen-constant size; the corner axes on a screen-anchored layer).

- [shader-analysis.md](shader-analysis.md): proposed: every reachable variant compiled on every backend, GLSL and WGSL
  reviewed side by side for parity and numerical hazards, the costliest passes measured and analysed.

## Dropped

- [batched-draws.md](batched-draws.md): many objects of separate representations per draw (merged geometry, object and
  material blocks by index, solo and merged programs); built, then removed: Mol* has at most about 1k such objects,
  many units of a representation are served better by merging (`units-1000-spacefill` merged 0.43 ms, batched
  2.43 ms per frame).
- [ssao-blur-depth.md](ssao-blur-depth.md): linear depth in the SSAO target is good enough; reading the scene depth in
  the blur (slower) only for banding beyond about 400 Å from the camera.
- [sphere-depth-prepass.md](sphere-depth-prepass.md): a depth pre-pass for sphere impostors; dropped after measuring,
  the depth-only pass costs about as much as the impostors without lighting; its read-only depth passes kept.
- [webgpu-sync-readback.md](webgpu-sync-readback.md): synchronous picking reads on WebGPU through a canvas copy; worked
  in Chrome and Firefox, dropped because clicks already wait for the asynchronous read (API callers await
  `Canvas3D.asyncIdentify(…).wait()`).

## Not started

- [webgpu-xr.md](webgpu-xr.md): XR on WebGPU through `XRGPUBinding` behind the neutral `GraphicsXR`.
- [primitive-expansion.md](primitive-expansion.md): cylinders, lines and text expanded in the renderer.
- [core-gpu-compute.md](core-gpu-compute.md): a neutral map-kernel tier and GPU algorithms in core.
- [attribute-color-size.md](attribute-color-size.md): per-vertex color and size attributes, to be considered.
- [dual-source-blending.md](dual-source-blending.md): dual-source blending where available, for features (colored
  transmission in blended transparency, per-pixel accumulation weights); not for performance, and not for the passes
  with multiple render targets (WBOIT, DPOIT).
