# Smoothed layers sampled on the GPU (merged surfaces)

Status: steps 1 to 3 and the follow-ups done; step 4 (the splat on the GPU): 4a to 4c done, 4d next. Follows [merged-units-core.md](merged-units-core.md) step
4b.

## Why

A smoothed layer (overpaint, transparency, emissive, substance) on a merged surface was sampled per vertex and instance
on the CPU, where per-group visuals sample their grid in the fragment shader. Measured (60 copies of 7qpd as a root
structure: 7 groups of 60 instances, CPU meshes at resolution 1, WebGPU; `merging/assembly-gaussian-smooth-layer-*`,
the layer set each frame): on all units 2.75 s a frame merged against 2.04 s per group, the same 35 MB of layer data; on
one unit 18.8 against 12.9 ms.

- **Shared by both paths**: the splat into the grid, on the CPU (about 2 s).
- **The merged path's own**: the trilinear sampling at every vertex of every reached instance (8.9 M samples), and the
  35 MB written to the CPU arrays and uploaded.

## Design

- **WebGPU with compute only** (`isComputeLayerSamplingSupported`). WebGL and contexts without compute keep the CPU
  path (`getTrilinearlyInterpolated`).
- **The grid as before**: splatted on the CPU (`calcMeshColorSmoothing` with the context), uploaded as a texture (as
  the per-group grids), with the instances it reaches. The merged visual asks for it per vertex
  (`Visual.SmoothingContext.perVertexOnGpu`, `applyMesh*Smoothing`).
- **A kernel samples the vertices** (`layer-sample-vertices`, `gl/compute/layer-sampling.ts`), per reached instance
  and vertex: the segment's positions (on the GPU once a geometry, since 4b) and transforms (uploaded a layer), the grid sampled by
  the texture's linear filter as the mesh shaders sample it (`texture3dFrom2dLinear`), the result written into the
  storage the layer is drawn from, at the cells' corners as the splat puts them (as the shaders since the half-cell
  fix). `r8` layers pack four elements a `u32`: a thread a word, the words at the range's
  ends read and their other bytes kept.
- **The GPU owns those ranges.** The segment's layer is a `DataArray` with a writer (`DataArray.gpu`, a
  `GpuDataWriter`), its array zeros shared by all such arrays (`getSharedZeros`: never written, never reused by
  `createDataArray`): no CPU copy of the 35 MB. The merged object puts the writer on the layer's part
  (`SegmentedArray.Part.gpu`); the WebGPU render item, after it writes a layer from the CPU, writes on the GPU each part
  whose range or the words around it it wrote (`writeGpuParts`): a part changed, a layout made again (a full write), a
  layer moved without its bytes copied. A layer moved by a GPU copy keeps them.
- **Released** where replaced (`releaseGpuData` in the layers' `create*`) or its segment dropped (`disposeSegment`);
  the storage buffers and the grid through the context's deferred release.
- **CPU readers read back** (`readGpuDataArrays`): the geometry export reads a merged member's GPU-written layers by a
  dispatch into a buffer of its own, as drawn.
- **A new context** (a backend switch, a restore): a writer of another context's kernels makes its object foreign
  (`hasForeignResources`), the merged visual made again on the new one.

## Steps

1. **Baseline** (done): `merging/assembly-gaussian-smooth-layer-{one,all}` (merged) set the layer in one of two colors
   in turn, so that every frame smooths; measured above. They and `merging/assembly-gaussian-smooth-merged` pin
   resolution 1: `auto` quality gives the assembly 3, which nothing smooths at.
2. **The kernel** (done): `layer-sample-vertices`, dispatched for one mesh into a plain buffer, from offsets in its
   positions, transforms and output. Kernels may sample a texture filtered: a linear texture of a format every device
   filters gets a filterable binding and the kernel the sampler `sLinear` (as render shaders). Tested on WebGPU
   (headless) against a CPU model of the mesh shaders' sampling (`texture3dFrom2dLinear`): within 0.76 of a byte step
   for `rgba8` and `r8`, an `r8` range at an unaligned offset keeping its neighbors' bytes, an unreached instance zero.
3. **Merged objects** (done): as the design above. Tested on WebGPU (headless): the image against the per-group
   visuals' 35 pixels apart where the plain layer is 264 (the CPU path was 92: the filter's half cell, and the vertices'
   interpolation of it), the export's read back equal to the bytes drawn, a geometry made again (its segments laid out
   again) sampled again, a device loss drawing the same image on the new device (since the half-cell fix the plain layer
   288 pixels apart, merged still 35). Measured against step 1: on all units
   2.75 to 2.13 s a frame (per group 2.04), on one unit 18.8 to 11.3 ms (per group 12.9); the splat the rest.
4. **The splat on the GPU** (planned, see [The splat on the GPU](#the-splat-on-the-gpu-step-4)): 4a to 4e.

## The splat on the GPU (step 4)

The splat (`calcMeshColorSmoothing`, `splatVertex`) is what is left of a smoothed layer's frame on CPU meshes: about
2 s of the 2.13 on all units, merged and per group alike. The texture meshes splat on the GPU already
(`computeColorSmoothing`: `color-smoothing-accumulate` into fixed-point atomic sums, `color-smoothing-normalize` into the
atlas), from the isosurface's vertex arena; CPU meshes have no vertices on the GPU a kernel can read.

### Design

- **WebGPU with compute only**, as the sampling. WebGL and contexts without compute keep the CPU splat: CPU meshes have
  no position texture for the point-sprite accumulation (one uploaded a geometry would do), decided after 4e's
  numbers.
- **The vertices on the GPU once a geometry** (`GpuMeshVertices`, `geo/geometry/mesh/gpu-vertices.ts`: positions as
  `f32` xyz, groups as `u32` with 4c, storage buffers of the context): kept by the mesh's values (a segment's), uploaded
  at the first GPU smoothing and again for a new version of the position cell (a geometry made again, possibly into
  the same arrays) or another context (the old one's released then), released with the values (`disposeSegment`; the
  per-group visuals' with 4d). The writers resolve them at each write, not at their creation. The splat and the sampling read them, so a
  layer uploads only its transforms (64 bytes an instance), reached flags and colors. Not the render item's vertex
  buffers: the geometry pages are `VERTEX | INDEX`, not `STORAGE` (the usage on every page, and bindings of
  `maxStorageBufferBindingSize`), a slot's attributes may be interleaved, the vertex buffers are written at the first
  draw (after the layer's write) and the transforms are a vertex buffer or instance storage by context.
- **One splat, two vertex sources**: `color-smoothing-accumulate` (the texture meshes' arena, within the WebGL point
  sprite's cells) and `color-smoothing-accumulate-mesh` (the vertices above from a vertex offset, every `stride`-th from
  0 as the CPU: `ceil(vertexCount / stride)` samples), sharing the splat (`splat` in `color-smoothing.wgsl.ts`): the
  weights (`2 - distance` within 2) and the fixed-point sums as they were. Two kernels, not a define: their storage
  bindings differ. The
  colors as the texture meshes' path loads them (a `group` or `groupInstance` texture of the layer or the color theme),
  the instances splatted as `getSplatInstances` picks them (the subset's transforms and rows, as
  `getTextureMeshSplatInstances`).
- **The sums over the region** (`getSplatInstances`' region: the reached instances' cells), not the whole grid: the
  grid spans every instance (the shaders read it without a bounds check), so one unit of the assembly cleared and
  normalized the whole assembly's grid (20 bytes a cell). An offset uniform in both kernels; the atlas cleared, then
  normalized in the region; slabs where the region's sums exceed the storage limit, as now.
- **The grid as the CPU makes it**: `getGrid`'s lattice and `getVolumeTexture2dLayout`'s atlas (slices in rows of
  `width / gridDim[0]`, the normalize's `uTexCols`), the same `gridTransform` and resolution. An `rgba8` texture for
  `r8` layers too (`r8unorm` is not a storage format in core WebGPU): the shaders and the sampling read `.r`; the
  export reads the grid back with a stride by its format (`readTextureBytes`).
- **Who gains**: the per-group visuals (their `volumeInstance` grids, the fragment shaders' sampling) and the merged
  ones (the sampling kernel), layers and the color theme's smoothing (`applyMeshColorSmoothing`). The merged visual's
  color smoothing without a context (`processValues`, per vertex on the CPU) stays.
- **Equal to the CPU splat within the rounding**: the same cells and weights; the CPU sums bytes in `f32` in vertex
  order, the GPU unorm colors in fixed point (1/16384) in any order, so a cell's byte within 1; but where a cell's
  weight is near zero (one vertex at the splat's radius), the fixed point's step a part of it: measured one cell of
  27k off by 4 (its color uniform). Not seen; a larger scale would risk the sums' `u32`.

### Steps

- **4a. Baseline** (done): `merging/assembly-gaussian-smooth-layer-{one,all}` with `-per-group` (the mesh visual
  without its merge builder) and `-texture-mesh` (`tryUseGpu: true`, the default) variants, 4 frames a run (a frame of
  seconds), the render objects by type in their notes. A frame, the layer set each:

  | Layer on              | WebGPU merged | per group | texture meshes | WebGL2 merged | per group | texture meshes |
  | --------------------- | ------------- | --------- | -------------- | ------------- | --------- | -------------- |
  | all units             | 1.88 s        | 1.90 s    | 56 ms          | 2.75 s        | 2.23 s    | 1.06 s         |
  | one unit              | 13.4 ms       | 13.2 ms   | 7.5 ms         | 18.2 ms       | 11.5 ms   | 31 ms          |

  The texture meshes' compute splat 34 times faster than the CPU's on all units: step 4's target. Their WebGL point
  sprites half the CPU's on all units, slower on one (a pass a slice over the whole grid). Paths: on WebGPU 5 of the 7
  groups texture meshes, the 2 smaller CPU meshes (below `SmallSurfaceCells` and `SmallSurfaceAtoms`), merged; on
  WebGL2 all 7 texture meshes. On CPU meshes the frame is the splat (per group nothing else: the fragments sample).
  Tested (headless WebGPU, `webgpu-texture-mesh-smoothed-layers`): a texture mesh's layer drawn near the CPU mesh's
  (134 to 171 pixels from run to run against the plain layer's 283; the surfaces 32 apart without a layer), its grid
  non-zero only around the reached instance of three. Its grid not as the CPU's: another lattice (`getGridDim` on the
  sphere's box, not `getGrid`'s cells), every ninth arena vertex (never indexed) in the order of the kernels' appends;
  4c's kernels on CPU meshes use the CPU's grid and vertices.
- **4b. The vertices on the GPU** (done): `GpuMeshVertices` and its release; the sampling kernel reading them (no
  upload a layer). Tested (headless WebGPU): two layers sampled from one upload, the positions' bytes released with the
  representation (not without `releaseGpuMeshVertices`); step 3's tests (bytes as drawn, a geometry made again, a device
  loss) unchanged. Not measurable on the benchmark: the splat's 2 s hides the upload (all units merged 2.19 to 2.31 s
  against per group 2.17 to 2.25 in the same runs).
- **4c. The kernels** (done): `color-smoothing-accumulate-mesh` (`computeMeshColorSmoothing`), the region
  (`uRegionMin`, `uRegionMax`) in both accumulates and the normalize (the texture meshes': the whole grid, as before),
  `GpuMeshVertices` with the groups; `calcMeshColorSmoothing` splatting on the GPU where its input has the mesh's
  vertices there (`gpuVertices`, of the context's compute): the same grid, instances and region, an rgba texture (a
  single-channel one of the CPU's not reused). Tested (headless WebGPU, `webgpu-mesh-color-smoothing`): three copies, a
  layer on one, `rgba8` and `r8` grids as the CPU's (layout, transform, instances splatted), within a byte but cells of
  near-zero weight (at most one in a thousand, by up to 8), none outside the region; the texture meshes' visual cases
  (17 on WebGPU) and their test unchanged.
- **4d. The meshes' smoothing on the GPU**: `calcMeshColorSmoothing` with a compute context splatting by the kernels,
  per group and merged, colors and layers. Tested: images as step 3's (merged and per group against the CPU splat),
  the export's grids read back, the visual geometry cases on WebGPU.
- **4e. Measured** against 4a, merged and per group, one unit and all; then WebGL decided (the CPU splat kept, or
  positions as a texture for the point sprites).

## Follow-ups

- **The half-cell shift** (done): the splats (CPU and GPU) put a cell's value at its corner (`p = i`), the filter reads
  texel `i` at its center (`p = i + 0.5`); the mesh shaders, WebGL (`texture3dFrom2dLinearCorners`) and WebGPU
  (`texture3dFrom2dLinear_<texture>`), and the kernel sampled half a cell off in x and y (z mixed from the corners
  already). Now at `p + 0.5`: the GPU sampling equal to the CPU path's (`getTrilinearlyInterpolated`) within the
  rounding, tested. The direct volume adds the half itself (its grids are of cell centers).
- **The kernel's pipeline** (done): compiled in the background once a merged surface smooths
  (`prepareLayerSampling`), and by a precompile list with a merged mesh (`Canvas3DPreparer.precompile`,
  `prepareObjects`: whether it smooths is not in its signature); else synchronously at the first layer.
