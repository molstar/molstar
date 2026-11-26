# Smoothed colors sampled on the GPU (merged surfaces)

Status: done (WebGPU; WebGL keeps the CPU path). Follows [merged-smoothed-layers-gpu.md](merged-smoothed-layers-gpu.md), whose path it gives the color
theme's smoothing.

## Why

A merged surface smooths its colors on the CPU (`processValues` without a context, `units-visual.ts`: a grid per object
does not merge, so data per vertex), where per-group surfaces splat on the GPU and sample their grid in the fragment
shader; the GPU path built for the layers (step 4 of the layers' plan) left the colors out. Measured (60 copies of 7qpd
as a root structure, CPU meshes at resolution 1, colored by element-index, WebGPU; `merging/assembly-gaussian-smooth-*`,
the representations' build):

| Colors       | merged   | per group |
| ------------ | -------- | --------- |
| smoothed     | 2711 ms  | 215 ms    |
| not smoothed | 208 ms   |           |

So about 2.5 s of a merged surface's build is its colors' smoothing, paid again by every build: a theme change, a
surface's props, a trajectory frame. Per group the same smoothing costs about 7 ms. Frames the same (2.1 to 2.7 ms).

## Design

- **The layers' path for `tColor`**: the splat on the GPU (`computeMeshColorSmoothing`, the mesh's vertices on the GPU,
  `GpuMeshVertices`), the grid sampled per vertex by `layer-sample-vertices` into the storage the merged object draws
  `tColor` from (`DataArray.gpu`, `createGpuLayerData`), its array the shared zeros. WebGPU with compute only
  (`isComputeLayerSamplingSupported`); WebGL and contexts without compute keep the CPU path.
- **Colors by group too**: a theme by group (`group`, the copies alike) smooths with one instance into `vertex` data,
  not `vertexInstance`: the sampling with one identity transform and `vertexCount` elements, not times the instances.
  Instance colors (`groupInstance`) as the layers (`vertexInstance`).
- **Where the context is**: the merged visual smooths in `processValues` while it makes a segment and as an update
  writes it, without a context (`createSegment`, `updateSegment`). It gives the context and asks per vertex
  (`perVertexOnGpu`, as `applyMesh*Smoothing` take it): `processValues` takes the option, `applyMeshColorSmoothing` a GPU
  branch as `applyMeshOverpaintSmoothing`'s. The per-group visuals unchanged (their grids).
- **What the merged object reads as before**: the color type `vertex` or `vertexInstance` (the merge key the same as on
  the CPU path), the data a segment's record cell, the export reading it back (`readGpuDataArrays`).
- **Released with the segment** (`disposeSegment`: `tColor`'s writer, as the layers'), on a context change made again
  (the segments' rebuild).

## Steps

1. **The kernel's vertex case**: `createGpuLayerData` (or its sampling input) for data per vertex alone (`vertex`: one
   identity transform, `vertexCount` elements); tested against the CPU's `getTrilinearlyInterpolated` on a grid.
2. **`applyMeshColorSmoothing` per vertex on the GPU** (`perVertexOnGpu`): the splat and the sampling as the layers';
   tested as the overpaint's (the same per-vertex colors within a byte as the CPU path, by group and by instance).
3. **The merged visual**: `processValues` with the context and `perVertexOnGpu` where it makes and updates a segment;
   `tColor`'s writer released with the segment; tested: a merged surface's smoothed colors drawn as before (an image
   against the CPU path's), read back by the export, a context restored.
4. **Measured** (done): the merged build 2711 to 236 ms on WebGPU, as without smoothing (260) and per group (233);
   frames the same (1.9 to 2.7 ms). WebGL2 on the CPU path, 4171 ms. The textures 99 to 132 MB: each segment's writer
   keeps its grid (for the data's writes again, as the layers'), beside the per-vertex data; a grid sized to the
   reached instances (`regionGrid`, [merged-smoothing-efficiency.md](merged-smoothing-efficiency.md)), which colors
   by element or chain reach all of.
