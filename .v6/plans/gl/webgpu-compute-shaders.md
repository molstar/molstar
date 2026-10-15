# Plan: WebGPU compute shaders (infrastructure)

## Status

- Planned. The infrastructure for running WGSL compute shaders on the WebGPU backend: kernels, storage buffers and
  storage textures, compute passes, timing, readback, limits and tests. The algorithm ports built on it (marching cubes,
  Gaussian density, color smoothing, grid3d) are in [webgpu-compute.md](webgpu-compute.md). Related:
  [core-gpu-compute.md](core-gpu-compute.md) (a later, core-level kernel interface),
  [storage-buffers.md](storage-buffers.md) (data arrays of render programs), [gpu-timer.md](gpu-timer.md),
  [resource-readiness.md](resource-readiness.md), [wgsl-strategy.md](wgsl-strategy.md) (section Compute).

## Current state

- All GPU compute is fragment-shader compute: `createComputeRenderItem` draws a quad into a render target, the WGSL
  ports so far (`copy`, `ssao`, the passes) are render pipelines. There are no compute pipelines, storage buffers or
  storage textures on the WebGPU backend.
- What exists to build on:
  - WGSL declarations generated from schemas (uniform blocks in std140 = WGSL uniform layout, texture and sampler
    bindings, `getBindingLayout`), preprocessing with defines, override constants from number defines
    (`getOverrideConstants`), program registry by `ShaderCode` name (`webgpu/shader/programs.ts`) and
    `isProgramSupported`.
  - The context's command encoder with `encode` (ends the open render pass for copies), `submit` per task, timestamp
    writes per pass encoder (`timerSource.getTimestampWrites`), `readGPUTexture` (staging buffer + `mapAsync`),
    `resources.pending` for asynchronous work.
  - The backend-neutral rule: code outside `gl/webgl/` and `gl/webgpu/` imports neither backend
    (`backend-neutral.test.ts`); WGSL lives in `gl/webgpu/`.

## Design

### Neutral API (follows WebGPU, absent on WebGL)

- `GraphicsContext.compute?: ComputeContext`, defined only where compute shaders exist (WebGPU); `capabilities.compute`
  mirrors it. Algorithms branch on it and keep their WebGL fragment-shader paths unchanged (byte-identical).
- Kernels are named like programs: `KernelCode(name, { workgroupSize })` is neutral (no source); the WGSL source is
  registered by name in `gl/webgpu/compute/kernels.ts`, like `programs.ts`. `compute.isKernelSupported(code)` like
  `isProgramSupported`.
- `compute.createKernel(code, schema, values) => ComputeKernel` with `dispatch(x, y?, z?)`,
  `dispatchIndirect(buffer, offset)`, `getByteCount()`, `destroy()`; values are `ValueCell`s as for render items
  (uniform and binding updates by version, as render items do).
- Resources: `compute.buffer(byteLength, usage) => StorageBuffer` (`write(data, byteOffset?)`, `copyTo(...)`,
  `readAsync(byteOffset?, byteLength?) => Promise<ArrayBuffer>`, `destroy`), counted in `stats.resourceCounts` and
  `getByteCounts`. Textures gain a `storage` usage flag at creation (only for formats with storage support:
  `rgba8unorm`, `rgba16float`, `rgba32float`, `r32float`, `r32uint`, ...).

### Schema and declarations

- New specs next to the existing ones: `StorageBufferSpec(element, access)` with `access` `'read' | 'read-write'` and
  WGSL element types (`u32`, `i32`, `f32`, `vec4<f32>`, `atomic<u32>`, or a struct name declared in the kernel),
  declared as runtime-sized arrays; `StorageTextureSpec(format, dimension, access)` (`'write' | 'read' | 'read-write'`).
  Uniforms (`UniformSpec`) and sampled textures (`TextureSpec`) as for programs.
- Declarations (`var<storage, ...>`, `texture_storage_*`, uniform blocks, samplers) and the bind group layouts are
  generated from the schema, as for render programs; kernels contain no hand-written binding declarations. One bind
  group per kernel for now (no tiers), rebuilt when a bound resource changes.
- Workgroup size from the `KernelCode` (`@workgroup_size` from an override constant, so one module serves several
  sizes); number defines as override constants as for programs.

### Encoding, ordering and timing

- Compute passes follow render passes: `ctx.beginComputePass({ label })` ends an open render or compute pass and begins
  a compute pass encoder (lazily, on the first dispatch); dispatches go into the current compute pass;
  `beginRenderPass`, `encode` and `submit` end it. Labels are required, as for render passes.
- Timing: compute passes get `timestampWrites` exactly like render passes (pass leaves, summed by `group`), no
  per-dispatch timestamps (consistent with [gpu-timer.md](gpu-timer.md)).
- Ordering: everything is recorded on the one command encoder, so dispatches, copies and render passes run in recording
  order; `queue.writeBuffer`/`writeTexture` run before the submit of the encoder that uses them (as for textures today).
  No explicit barriers exist in WebGPU (usage scopes per pass handle hazards); a buffer can't be both writable storage
  and read in the same dispatch (validation), which the binding generator checks early.

### Readback and indirect use

- `StorageBuffer.readAsync` submits, copies into a pooled staging buffer and maps it (as `readGPUTexture`); tracked as
  pending work so screenshots wait. No synchronous readback (WebGPU has none, see
  [webgpu-sync-readback.md](webgpu-sync-readback.md)).
- Counts produced on the GPU (e.g. marching cubes vertex count) can stay on the GPU: `dispatchIndirect` for follow-up
  kernels and, as a later step, indirect draws for render items (`drawCount` from a buffer), avoiding the round trip.

### Limits and validation

- WebGPU core defaults are the minimum spec ([minimum-capabilities.md](minimum-capabilities.md)): 256 invocations per
  workgroup (x, y ≤ 256, z ≤ 64), 16 KiB workgroup storage, 65535 workgroups per dimension (larger grids dispatch in 2D
  and linearize), 8 storage buffers and 4 storage textures per stage, 128 MiB per storage binding. Unit tests check
  every registered kernel against them, like `minimum-limits.test.ts` for programs.
- No float atomics in WGSL: kernels use `atomic<u32>`/`atomic<i32>` (fixed-point sums, `atomicMin` on non-negative float
  bits).
- Errors: shader compilation info reported with WGSL locations; validation errors through the existing error handling
  (debug mode).

### Tests

- Unit tests: declaration generation for each spec kind, bind group layouts, limits per kernel, preprocessing of each
  kernel with its define variants.
- GPU conformance in the browser: visual harness cases of a new `compute/*` group that run small kernels (fill, prefix
  sum, atomic min/add, storage texture write) and compare the read-back results with CPU references (pass/fail, no
  image); skipped on WebGL.

## Steps

1. Neutral types (`ComputeContext`, `KernelCode`, `StorageBuffer`, `StorageBufferSpec`, `StorageTextureSpec`,
   `capabilities.compute`) and the kernel registry; WebGL: `compute` undefined.
2. WebGPU storage buffers (create, write, copy, `readAsync` with pooled staging buffers, stats) and storage usage on
   textures.
3. Kernels: declarations and bind group layouts from schemas, pipeline cache (module per preprocessed source, pipeline
   per override constants), `dispatch`/`dispatchIndirect`, value updates.
4. Compute passes on the context (`beginComputePass`, interplay with render passes, `encode` and `submit`) and timestamp
   writes; timer unit tests.
5. Limits tests and the `compute/*` conformance cases (fill, scan, atomics, storage texture).
6. Optional: indirect draws for render items (`drawCount` from a storage buffer).

## Decisions (open)

1. Placement of the API: `GraphicsContext.compute` (neutral, optional, WGSL by kernel name in `gl/webgpu/`), or a
   WebGPU-only module that the algorithms import behind `isWebGPUContext` (breaks the backend-neutral rule for the
   algorithm files).
2. Pass grouping: explicit `beginComputePass({ label })` grouping dispatches (as render passes, one timestamp pair per
   pass), or one compute pass per dispatch (simpler, more passes and timestamps).
3. GPU-produced counts: read back first and indirect draws later, or indirect draws from the start.
