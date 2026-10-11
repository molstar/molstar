# Plan: GPU compute interface in core

## Status

- Not an immediate goal; keep in mind for the backend split ([backend-split.md](backend-split.md)) and WebGPU work.
- Background: the workspace split ([workspace-prototype.md](../workspace-prototype.md), layer direction `plugin -> graphics -> model -> io -> core`) moved GPU helpers out of `core/math` into graphics because they depend on WebGL, e.g. the GPU Gaussian density (`computeGaussianDensityTexture*`, `GaussianDensityTextureData` now in `@molstar/graphics/geo/gaussian-density`; the CPU version stayed in `@molstar/core/math/geometry/gaussian-density`, see [breaking-v6-changes.md](../breaking-v6-changes.md)).
- Goal: a small, ideally type-only GPU compute interface in core, so GPU-accelerated algorithms can live in `core/math` (and model, io) again. Hard rule: core never imports code from graphics (no runtime imports and no type imports; core's package does not depend on graphics). Graphics implements the interface; callers pass the implementation in.

## Inventory

| Algorithm | Today | Result | Fit for core |
|---|---|---|---|
| Gaussian density | CPU in core, GPU in `graphics/geo/gaussian-density{,/gpu}.ts` | density/min-distance/group grids; as texture for rendering, or read back | Grid values (read back) yes; texture output for GPU isosurface stays in graphics |
| Grid evaluation (`createGrid3dComputeRenderable`, `gl/compute/grid3d.ts`) | graphics | per-cell values, read back (`Float32Array`) | Yes, the generic primitive (used by the alpha-orbitals extension, which would then not need graphics for compute) |
| Marching cubes (GPU: histogram pyramid, active voxels, isosurface) | graphics | texture-mesh buffers for rendering | No, render-integrated |
| Marching cubes (CPU, `geo/util/marching-cubes`) | graphics (produces `Mesh`) | mesh | Separate question (output type is graphics geometry) |
| Color smoothing (mesh/texture-mesh) | graphics | color textures | No, render-integrated |
| Future candidates | - | e.g. molecular surface fields, distance/contact computations, spatial lookups | Yes if results are read back or consumed through opaque handles |

Rule of thumb: algorithms with CPU-visible results (or GPU results consumed only through opaque handles) belong in core; anything producing render resources stays in graphics.

## Design

- `@molstar/core/gpu` (type-only module, `export type`/`interface` only, no runtime code):
  - `GpuCompute`: a subset of the capabilities (`GraphicsCapabilities` keys relevant for compute, e.g. float render targets, 3D textures), limits, and kernel/data creation.
  - Data: create from flat typed arrays (consistent with `DataSpec`, [storage-buffers.md](storage-buffers.md)), opaque handles for GPU-resident inputs/outputs, async readback (`Promise<TypedArray>`, as WebGPU requires, see backend-split step 5).
  - Kernels: a neutral kernel spec close to grid3d: an index domain (1D/3D grid), typed inputs (uniforms, data arrays, grid textures), and a per-element body. The body is shader source supplied by the algorithm, per language (`{ glsl, wgsl }`) until the WGSL strategy decides on translation; WebGL constraints (constant loop bounds via defines, as grid3d's `loopBounds`) are part of the spec.
  - Cancellation/progress via the existing `RuntimeContext`.
- Graphics implements `GpuCompute` on its `GraphicsContext`: WebGL by rendering into float or packed RGBA8 targets (the grid3d mechanism, generalized), WebGPU with compute shaders and storage buffers. Kernel declarations are generated from the kernel spec like the schema-generated declarations today.
- Algorithms in `core/math` take an optional `gpu?: GpuCompute` argument (explicit injection, no global registry), check capabilities, and fall back to the CPU path; results match the CPU path within a documented tolerance (float32 on GPU).
- Graphics keeps thin adapters where rendering needs GPU-resident results (e.g. density texture into GPU marching cubes), built on the same handles.

## Steps (when picked up)

1. Prerequisites from the backend split: `GraphicsContext` with capabilities (step 1), neutral resources (step 3), async readback (step 5).
2. Guard test: no imports of `@molstar/graphics` (or `plugin`) from core, io and model sources, including `import type`; the package dependency graph already forbids it at build time.
3. Define the type-only `GpuCompute` interface and kernel spec in core; conformance tests in graphics (headless-gl for WebGL; WebGPU when available).
4. Implement on WebGL by generalizing grid3d; port the alpha-orbitals compute to it (extension no longer imports `gl/compute/grid3d`).
5. Move the Gaussian density grid computation back to `core/math` with the GPU path via `GpuCompute`; graphics keeps the texture-output variant for GPU surfaces.
6. Document the moves in [breaking-v6-changes.md](../breaking-v6-changes.md) (or the v7 equivalent).

## Open questions

- Strictly type-only, or a tiny runtime part in core (e.g. kernel spec helpers/validation)?
- Kernel body language: GLSL + WGSL per kernel vs. one language with translation (depends on the WGSL strategy plan).
- Opaque GPU handles across packages: lifetime/ownership (who destroys), context loss.
