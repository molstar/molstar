# Plan: GPU compute on WebGPU (Gaussian surface, marching cubes, grid3d)

## Status

- **Current**: done, steps 1–6: marching cubes, Gaussian density, texture mesh color smoothing and grid3d on WebGPU,
  chosen by capability (no backend checks), with indirect draws ([indirect-draws.md](indirect-draws.md)); visual cases
  `*-gpu*`, `grid3d-compute*` and `fallback/*-gpu` pass. The 3d Gaussian density, once the slow part, runs as compute
  kernels and is faster than WebGL's from medium sizes on; current numbers in
  [marching-cubes-benchmark.md](marching-cubes-benchmark.md) (Current numbers). Nothing open.
- Gaussian density splatting done: the `gaussian-density` WGSL program renders per atom an instanced quad into the
  slices of a 3d texture (render attachments of 3d texture slices, `depthSlice`, on the WebGPU backend); WebGL draws the
  same quads instead of point sprites (byte-identical on all profiles, and no longer limited by the maximum point size).
  The GPU gaussian density takes a `GraphicsContext`; a first computation whose pipelines were still compiling renders
  again once they are ready (`renderUntilComplete`). `gaussian-density-volume` runs on WebGPU (visual case
  `geometry/direct-volume/gaussian-volume`, within one level of WebGL). The compute infrastructure is done (steps 1–5 of
  [webgpu-compute-shaders.md](webgpu-compute-shaders.md), with the prefix sum). Marching cubes done: count, prefix sum
  and emit kernels (`gl/compute/marching-cubes/compute.ts`), the total read back to size the output textures; the volume
  isosurface (values in a storage buffer cached on the volume) and the gaussian surface (unit and structure visuals, on
  the 3d density with the group packed in rgb) use them on WebGPU; visual cases `isosurface-gpu`,
  `isosurface-gpu-negative` and `gaussian-surface-gpu` match WebGL. Benchmark in
  [marching-cubes-benchmark.md](marching-cubes-benchmark.md): ×1.4–3.4 faster than WebGL on volumes from 128³; the 3d
  density is the slow part of WebGPU gaussian surfaces. Follows [webgpu-viewer.md](webgpu-viewer.md) milestone 4
  (everything else of the Viewer runs on WebGPU). Related: [wgsl-strategy.md](wgsl-strategy.md) (section Compute),
  [core-gpu-compute.md](core-gpu-compute.md), [storage-buffers.md](storage-buffers.md), [gpu-timer.md](gpu-timer.md).

## Limits of the WebGPU paths

WebGPU's 2d texture limit is often half of WebGL's (8192 against 16384), so WebGL criteria sized by it (the 2d atlas of
slices) wrongly exclude or coarsen WebGPU paths, which use 3d textures and storage buffers, and they do so silently (a
CPU fallback or a coarser surface still renders). Audit of the size limits that WebGPU paths pass through:

| Site                                                                                 | WebGPU                                                                                                                                                                    | Status                                                |
| ------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------- |
| Gaussian surface suitability (`suitableForGpu`)                                      | 3d grid: each side within `max3dTextureSize`, at most `maxStorageBufferSize / 4` cells (the marching cubes counts, `getComputeDensityMaxCells`; 2²⁵ with the core limits) | fixed (7qpd took the CPU path)                        |
| Gaussian density resolution (`ensureReasonableResolution` with `getTextureMaxCells`) | all cells of the 3d grid (`ensureComputeResolution`)                                                                                                                      | fixed (units coarsened)                               |
| Volume isosurface suitability                                                        | the volume values in one storage buffer, `compute.maxStorageBufferSize`                                                                                                   | was the core default; now the device limit            |
| Color smoothing grid (2d atlas, both backends)                                       | coarsened until the atlas fits `maxTextureSize` (`getFittingResolution`)                                                                                                  | added (larger grids failed to allocate)               |
| Marching cubes output textures (`getOutputTextureSize`)                              | 2d textures, `maxTextureSize`                                                                                                                                             | correct                                               |
| Direct volume (3d texture, half of `max3dTextureSize`)                               | as WebGL2                                                                                                                                                                 | correct                                               |
| Volume segmentation                                                                  | the smoothed segment mask (f32) in a storage buffer, `compute.maxStorageBufferSize`                                                                                       | ported: compute marching cubes, drawn indirectly      |
| grid3d (alpha orbitals)                                                              | a kernel per spec with WGSL code (`wgsl`), values in a storage buffer                                                                                                     | ported: `canComputeGrid3dOnGPU(graphics, computable)` |
| 2d Gaussian density                                                                  | WebGL without 3d textures                                                                                                                                                 | not needed (WebGPU has 3d textures)                   |

The check cases `gpu-paths/representations` (including volume segments), `gpu-paths/alpha-orbitals` (grid3d against the
CPU collocation, angular momenta 0–4) and `gpu-paths/gaussian-surface-units` assert that the GPU paths are taken and
match the CPU per unit (support as reported without strict mode).

## Inventory (WebGL, fragment-shader compute)

| Stage                                                        | Today                                                                                                              | Output                                                 | Consumers                                                     |
| ------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------ | ------------------------------------------------------------- |
| Gaussian density (`geo/gaussian-density/gpu.ts`, 667 lines)  | splats per atom into 2d atlas or 3d texture slices, additive blending for density, min blending for distance/group | rgba texture (density, packed group) + transform, dims | GPU marching cubes; `gaussian-density-volume` (direct volume) |
| Active voxels (`gl/compute/marching-cubes/active-voxels.ts`) | per cell: marching cubes case and vertex count                                                                     | texture                                                | histogram pyramid                                             |
| Histogram pyramid (`gl/compute/histogram-pyramid/*`)         | log2 reduction levels, sum readback                                                                                | pyramid textures, total count                          | isosurface                                                    |
| Isosurface (`gl/compute/marching-cubes/isosurface.ts`)       | per output vertex: pyramid traversal, edge interpolation, normal, group (MRT)                                      | `tPosition`, `tNormal`, `tGroup` float textures        | `TextureMesh`                                                 |
| Texture mesh rendering                                       | `mesh` program with `dGeometryType_textureMesh` (vertices read from textures by vertex id)                         |                                                        | not ported to WGSL                                            |
| Color smoothing (`gl/shader/compute/color-smoothing`)        | accumulate/normalize per-vertex colors into a grid                                                                 | color texture                                          | `TextureMesh` with smooth colors                              |
| grid3d (`gl/compute/grid3d.ts`)                              | generic per-cell kernel, caller-supplied GLSL body, packed or float output, read back                              | `Float32Array`                                         | alpha orbitals extension, visual tests                        |

## Design

- **Backend-specific implementations behind the existing algorithm entry points**: `computeGaussianDensityTexture`,
  `extractIsosurface`, the texture mesh color smoothing and grid3d keep their signatures (taking `GraphicsContext`) and
  dispatch to the WebGL implementation (unchanged, byte-identical) or a WebGPU one in `gl/webgpu/compute/`. The
  `isWebGLContext` gates become capability checks (`graphics.capabilities.compute` or per-algorithm support). The
  generic kernel interface of [core-gpu-compute.md](core-gpu-compute.md) stays a separate, later step.
- **WebGPU compute infrastructure**: kernels, storage buffers and textures, compute passes, timing and readback are
  planned separately in [webgpu-compute-shaders.md](webgpu-compute-shaders.md); this plan builds on it.
- **Data stays on the GPU between stages** as storage buffers; only the stage outputs consumed by rendering are textures
  (density for direct volumes, the texture mesh `tPosition`/`tNormal`/`tGroup`), written as storage textures
  (`rgba32float`, `rgba8unorm`) or copied from buffers.
- **Marching cubes**: per-cell case/count pass, a prefix sum (work-efficient scan over workgroups) for compaction, an
  emit pass writing vertices, normals and groups; the total vertex count read back asynchronously (the texture mesh
  needs the draw count; the WebGL path reads the pyramid sum synchronously).
- **Gaussian density**: see decision 4.
- **Texture mesh rendering**: the WGSL `mesh` program reads position, normal and group from the textures by vertex index
  (`dGeometryType_textureMesh`), no vertex buffers.
- **Readiness**: dispatched work is ordered before rendering by the queue; the asynchronous count readback is tracked as
  pending work ([resource-readiness.md](resource-readiness.md)), so screenshots wait for it.

## Steps

1. Texture mesh rendering in WGSL (done: `geometry/texture-mesh/primitives`); the steps below need the compute
   infrastructure of [webgpu-compute-shaders.md](webgpu-compute-shaders.md), or the port of decision 1.
2. Marching cubes on volumes (`geometry/texture-mesh/isosurface-gpu`): active voxels/count, scan, emit; `isosurface.ts`
   takes the GPU path on WebGPU.
3. Gaussian density (`geometry/texture-mesh/gaussian-surface-gpu`, `gaussian-density-volume`).
4. Texture mesh color smoothing (`gaussian-surface-gpu-smooth-colors`).
5. grid3d (`geometry/mesh/grid3d-compute*`) and the visual cases on the neutral context.
6. Viewer check (timer groups, large structures) and the `fallback/*-gpu` cases.

Check per step: the visual case matches the WebGL reference within tolerance (vertex order may differ from the histogram
pyramid's, so pixel equality is expected but not vertex equality); WebGL profiles byte-identical.

## Decisions

Decided:

1. Approach: mixed (1C), Gaussian density as render-pipeline splatting, marching cubes with compute shaders.
2. API: neutral and optional, `GraphicsContext.compute` (2A).
3. Passes: explicit `beginComputePass({ label })` grouping (3A).
4. Gaussian density: render-pipeline splatting (4A), follows from 1C.
5. Counts: read back first (5A); indirect draws are wanted as soon as they are available (5B as a follow-up, not
   optional).
6. grid3d: deferred (6B), then done as 6A: specs with WGSL code (`wgsl`) run as kernels on WebGPU.

The options and trade-offs considered, with the recommendation each, for reference. Decision 1 shapes the others.

Context: on WebGPU the Gaussian surface, isosurface and Gaussian density volume representations take the CPU paths
today, so this is about speed on large structures and volumes, not about missing features. Texture mesh rendering in
WGSL (step 1) is done.

### 1. Approach: port the fragment-shader compute, or use compute shaders

WebGL computes with render passes into float textures (splatting with blending, the histogram pyramid as reduction
passes, the isosurface emitted per output vertex with multiple render targets).

- **A. Port as is**: the same passes as WGSL render programs, like the postprocessing passes. No new infrastructure
  (render targets, float blending and render passes exist); results close to WebGL (same algorithm, vertex order
  included). Downsides: many passes (one per pyramid level), the pyramid sum is read back synchronously on WebGL and
  would become an asynchronous readback; WebGPU's strengths (scans, atomics, storage buffers) unused; the code stays as
  complex as the WebGL version.
- **B. Compute shaders**: the infrastructure of [webgpu-compute-shaders.md](webgpu-compute-shaders.md) (kernels, storage
  buffers and textures, compute passes, timing, readback), then marching cubes as count, scan and emit kernels. Simpler
  and faster algorithms, and the base for a core GPU compute interface ([core-gpu-compute.md](core-gpu-compute.md)).
  Downsides: the infrastructure is the largest piece of new backend code since the render passes; vertex order differs
  from WebGL (pixel equality expected, not vertex equality).
- **C. Mixed**: Gaussian density as render-pipeline splatting (option 4A, no compute needed), marching cubes with
  compute shaders (B).

Recommendation: C. Splatting is a natural fit for the rasterizer and ports directly; marching cubes is where compute
shaders pay off (one scan instead of a pyramid of passes), and it brings the infrastructure in with its first real user.

### 2. Placement of the compute API (with compute shaders)

- **A. Neutral and optional**: `GraphicsContext.compute?: ComputeContext` (WebGPU only), kernels by name like programs,
  WGSL registered in `gl/webgpu/`; algorithms branch on `graphics.compute` and keep the WebGL paths unchanged.
- **B. WebGPU-only module**: algorithms import it behind `isWebGPUContext`, which breaks the backend-neutral rule for
  the algorithm files (`backend-neutral.test.ts`).

Recommendation: A. It follows the existing split (programs by name, WGSL in `gl/webgpu/`) and keeps the guard test
meaningful; it is also the shape a core interface would take.

### 3. Grouping dispatches into passes

- **A. Explicit**: `beginComputePass({ label })` groups dispatches, like render passes (one timestamp pair per pass,
  labels required).
- **B. Implicit**: one compute pass per dispatch.

Recommendation: A. Same model and timer behavior as render passes; B would add a timestamp pair and pass overhead per
dispatch (marching cubes has several per surface).

### 4. Gaussian density on WebGPU

- **A. Render-pipeline splatting** as on WebGL: per atom a quad on each slice in its radius, additive blending for the
  density, max blending for the distance, then the group pass; `rgba16float` targets (always blendable on WebGPU; WebGL
  uses half floats when available too), 3D textures rendered slice by slice (`depthSlice`). Closest to WebGL's results;
  no compute needed.
- **B. Compute scatter**: per atom, atomics into a fixed-point grid (no float atomics in WGSL), `atomicMin` on distance
  bits, then the group. Deterministic and blending-free, but the fixed-point scale limits precision and atoms in dense
  regions contend for the same voxels.
- **C. Compute gather**: per voxel, sum the atoms of the neighboring cells of a spatial grid (sorted atoms, built on the
  CPU or the GPU). Deterministic, no atomics, the nearest atom gives the group directly; needs the spatial grid and is
  the most work.

Recommendation: A now (with decision 1C); B or C later if profiling shows splatting is the bottleneck.

### 5. GPU-produced counts (marching cubes vertex count)

- **A. Read back first**: the total count read back asynchronously (tracked as pending work, so screenshots wait); the
  texture mesh is drawn once it arrives (one frame later than on WebGL).
- **B. Indirect draws from the start**: render items take their draw count from a GPU buffer, no round trip; needs
  indirect draws in the neutral render item API, and CPU code that uses the vertex count (bounding data, picking sizes,
  color smoothing) has to cope with not knowing it.

Recommendation: A, then B as a separate optimization (step 6 of the infrastructure plan).

### 6. grid3d (alpha orbitals extension)

- **A. Optional WGSL body**: callers supply a WGSL kernel body next to the GLSL one; the alpha orbitals extension gets
  one.
- **B. WebGL only for now**: WebGPU keeps the CPU path (as today).

Recommendation: B until the compute infrastructure exists, then A as the last step (it is the generic kernel that a core
interface would expose).

### Order

1. Gaussian density splatting in WGSL (render programs): `gaussian-density-volume` and the density input of the GPU
   Gaussian surface on WebGPU (done).
2. Compute infrastructure, steps 1–5 of [webgpu-compute-shaders.md](webgpu-compute-shaders.md) (done).
3. Marching cubes (count, scan, emit) on volumes, then on Gaussian density: isosurface and Gaussian surface on the GPU
   (done). Measured against the CPU and WebGL paths with [marching-cubes-benchmark.md](marching-cubes-benchmark.md).
4. Texture mesh color smoothing.
5. Indirect draws for the marching cubes texture mesh (draw count from the GPU buffer, no readback).
6. grid3d with WGSL bodies (done).
