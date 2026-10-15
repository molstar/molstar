# Plan: GPU compute on WebGPU (Gaussian surface, marching cubes, grid3d)

## Status

- Planned. Follows [webgpu-viewer.md](webgpu-viewer.md) milestone 4 (everything else of the Viewer runs on WebGPU).
  Today WebGPU takes the CPU paths: the GPU paths are gated by `isWebGLContext` (`useGpu` in `gaussian-surface-mesh.ts`,
  `isosurface.ts`, `gaussian-density-volume.ts`; `canComputeGrid3dOnGPU`), so this is a performance port, not a
  correctness gap. Related: [wgsl-strategy.md](wgsl-strategy.md) (section Compute),
  [core-gpu-compute.md](core-gpu-compute.md), [storage-buffers.md](storage-buffers.md), [gpu-timer.md](gpu-timer.md).

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
- **Gaussian density**: see decision 1.
- **Texture mesh rendering**: the WGSL `mesh` program reads position, normal and group from the textures by vertex index
  (`dGeometryType_textureMesh`), no vertex buffers.
- **Readiness**: dispatched work is ordered before rendering by the queue; the asynchronous count readback is tracked as
  pending work ([resource-readiness.md](resource-readiness.md)), so screenshots wait for it.

## Steps

1. Texture mesh rendering in WGSL (testable now with CPU-built texture meshes or a WebGL-free case); needs the compute
   infrastructure of [webgpu-compute-shaders.md](webgpu-compute-shaders.md) for the steps below.
2. Marching cubes on volumes (`geometry/texture-mesh/isosurface-gpu`): active voxels/count, scan, emit; `isosurface.ts`
   takes the GPU path on WebGPU.
3. Gaussian density (`geometry/texture-mesh/gaussian-surface-gpu`, `gaussian-density-volume`).
4. Texture mesh color smoothing (`gaussian-surface-gpu-smooth-colors`).
5. grid3d (`geometry/mesh/grid3d-compute*`) and the visual cases on the neutral context.
6. Viewer check (timer groups, large structures) and the `fallback/*-gpu` cases.

Check per step: the visual case matches the WebGL reference within tolerance (vertex order may differ from the histogram
pyramid's, so pixel equality is expected but not vertex equality); WebGL profiles byte-identical.

## Decisions (open)

1. Gaussian density on WebGPU: compute scatter per atom with integer atomics (fixed-point density sum, `atomicMin` on
   distance bits, then group), deterministic and blending-free; or render-pipeline splatting with blending as on WebGL
   (needs blendable float targets).
2. grid3d: callers supply a WGSL kernel body next to the GLSL one (alpha orbitals extension updated), or grid3d stays
   WebGL-only for now.
