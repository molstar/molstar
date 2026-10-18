# Plan: GPU compute interface in core

## Status

- **Current**: not started. Its prerequisites are done: the WebGPU compute interface and the texture mesh color
  smoothing ([webgpu-compute.md](webgpu-compute.md)); grid3d runs as a kernel on WebGPU (decision 6A there), so this
  plan is about the neutral map-kernel tier and core-level algorithms, not about missing GPU paths.
- Planned after the WebGPU texture mesh color smoothing ([webgpu-compute.md](webgpu-compute.md)). The WebGPU-only
  compute interface (`GraphicsContext.compute`, [webgpu-compute-shaders.md](webgpu-compute-shaders.md)) exists; this
  plan adds the neutral map-kernel tier both backends implement (see Tiers), which is also where the deferred grid3d
  decision of [webgpu-compute.md](webgpu-compute.md) (decision 6) lands.
- Background: the workspace split ([workspace-prototype.md](../workspace-prototype.md), layer direction
  `plugin -> graphics -> model -> io -> core`) moved GPU helpers out of `core/math` into graphics because they depend on
  WebGL, e.g. the GPU Gaussian density (`computeGaussianDensityTexture*`, `GaussianDensityTextureData` now in
  `@molstar/graphics/geo/gaussian-density`; the CPU version stayed in `@molstar/core/math/geometry/gaussian-density`,
  see [breaking-v6-changes.md](../breaking-v6-changes.md)).
- Goal: a small, ideally type-only GPU compute interface in core, so GPU-accelerated algorithms can live in `core/math`
  (and model, io) again. Hard rule: core never imports code from graphics (no runtime imports and no type imports;
  core's package does not depend on graphics). Graphics implements the interface; callers pass the implementation in.

## Inventory

| Algorithm                                                                 | Today                                                         | Result                                                                        | Fit for core                                                                                                      |
| ------------------------------------------------------------------------- | ------------------------------------------------------------- | ----------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------- |
| Gaussian density                                                          | CPU in core, GPU in `graphics/geo/gaussian-density{,/gpu}.ts` | density/min-distance/group grids; as texture for rendering, or read back      | Grid values (read back) yes; texture output for GPU isosurface stays in graphics                                  |
| Grid evaluation (`createGrid3dComputeRenderable`, `gl/compute/grid3d.ts`) | graphics                                                      | per-cell values, read back (`Float32Array`)                                   | Yes, the generic primitive (used by the alpha-orbitals extension, which would then not need graphics for compute) |
| Marching cubes (GPU: histogram pyramid, active voxels, isosurface)        | graphics                                                      | texture-mesh buffers for rendering                                            | No, render-integrated                                                                                             |
| Marching cubes (CPU, `geo/util/marching-cubes`)                           | graphics (produces `Mesh`)                                    | mesh                                                                          | Separate question (output type is graphics geometry)                                                              |
| Color smoothing (mesh/texture-mesh)                                       | graphics                                                      | color textures                                                                | No, render-integrated                                                                                             |
| Future candidates                                                         | -                                                             | e.g. molecular surface fields, distance/contact computations, spatial lookups | Yes if results are read back or consumed through opaque handles                                                   |

Rule of thumb: algorithms with CPU-visible results (or GPU results consumed only through opaque handles) belong in core;
anything producing render resources stays in graphics.

## Tiers

Neutral where both backends can deliver, and no further: WebGL fragment-shader compute can do per-element maps, scatter
with blending (sum, min, max) and multi-pass reductions, but not arbitrary writes, atomics, workgroup memory and
barriers, read-write buffers in one pass, indirect dispatch, or integer math on WebGL1.

1. **`ComputeContext`** (exists, WebGPU only, `capabilities.compute`): kernels, storage buffers and textures, compute
   passes. New GPU algorithms are written against it. WebGL does not implement it: emulating storage buffers with
   textures would leave every algorithm to find out which operations don't work.
2. **Map kernels** (`GpuCompute` below, both backends): a grid of elements (1d or 3d), uniforms, input arrays behind
   opaque handles (a storage buffer on WebGPU, a data texture on WebGL), one output array of up to 4 floats per element
   (more with multiple render targets on WebGL), asynchronous readback. WebGPU dispatches a compute kernel writing a
   storage buffer, WebGL renders a fragment shader into a float target (the grid3d mechanism, generalized). WebGL's
   limits are the contract: loop bounds as defines, output at the invocation's own index only, floats only. Kernel
   bodies in GLSL and WGSL until there is a translation step.
3. **Neutral algorithms with an implementation per backend**, for everything needing scatter, atomics or workgroup
   memory: prefix sum, marching cubes, Gaussian density, color smoothing. One entry point, the implementation chosen by
   capability (as the isosurface chooses the histogram pyramid or the compute kernels today).

What moves to tier 2: grid3d and its callers (the alpha orbitals extension), and possibly the normalize step of color
smoothing; marching cubes, the density and the accumulation stay tier 3.

## Design

- `@molstar/core/gpu` (type-only module, `export type`/`interface` only, no runtime code):
  - `GpuCompute`: a subset of the capabilities (`GraphicsCapabilities` keys relevant for compute, e.g. float render
    targets, 3D textures), limits, and kernel/data creation.
  - Data: create from flat typed arrays (consistent with `DataSpec`, [storage-buffers.md](storage-buffers.md)), opaque
    handles for GPU-resident inputs/outputs, async readback (`Promise<TypedArray>`, as WebGPU requires, see
    backend-split step 5).
  - Kernels: a neutral kernel spec close to grid3d: an index domain (1D/3D grid), typed inputs (uniforms, data arrays,
    grid textures), and a per-element body. The body is shader source supplied by the algorithm, per language
    (`{ glsl, wgsl }`) until the WGSL strategy decides on translation; WebGL constraints (constant loop bounds via
    defines, as grid3d's `loopBounds`) are part of the spec.
  - Cancellation/progress via the existing `RuntimeContext`.
- Graphics implements `GpuCompute` on its `GraphicsContext`: WebGL by rendering into float or packed RGBA8 targets (the
  grid3d mechanism, generalized), WebGPU with compute shaders and storage buffers. Kernel declarations are generated
  from the kernel spec like the schema-generated declarations today.
- Algorithms in `core/math` take an optional `gpu?: GpuCompute` argument (explicit injection, no global registry), check
  capabilities, and fall back to the CPU path; results match the CPU path within a documented tolerance (float32 on
  GPU).
- Graphics keeps thin adapters where rendering needs GPU-resident results (e.g. density texture into GPU marching
  cubes), built on the same handles.

## Steps

1. Prerequisites: `GraphicsContext` with capabilities, neutral resources, async readback (done in the backend split),
   the WebGPU compute interface (done).
2. Guard test: no imports of `@molstar/graphics` (or `plugin`) from core, io and model sources, including `import type`;
   the package dependency graph already forbids it at build time.
3. Define the type-only `GpuCompute` interface and kernel spec in core; conformance tests in graphics (headless-gl for
   WebGL; WebGPU when available).
4. Implement on WebGL by generalizing grid3d, and on WebGPU with a kernel per map writing a storage buffer (on
   `ComputeContext`); port the alpha-orbitals compute to it (the extension no longer imports `gl/compute/grid3d`, and
   runs on the GPU with WebGPU too).
5. Move the Gaussian density grid computation back to `core/math` with the GPU path via `GpuCompute`; graphics keeps the
   texture-output variant for GPU surfaces.
6. Document the moves in [breaking-v6-changes.md](../breaking-v6-changes.md) (or the v7 equivalent).

## Open questions

- Strictly type-only, or a tiny runtime part in core (e.g. kernel spec helpers/validation)?
- Kernel body language: GLSL + WGSL per kernel vs. one language with translation (depends on the WGSL strategy plan).
- Opaque GPU handles across packages: lifetime/ownership (who destroys), context loss.
