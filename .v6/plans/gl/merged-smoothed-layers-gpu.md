# Smoothed layers sampled on the GPU (merged surfaces)

Status: steps 1 to 3 done; step 4 and the follow-ups open. Follows [merged-units-core.md](merged-units-core.md) step
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
  and vertex: the segment's positions and transforms (uploaded once a layer, 12 bytes a vertex), the grid sampled by
  the texture's linear filter as the mesh shaders sample it (`texture3dFrom2dLinear`), the result written into the
  storage the layer is drawn from. `r8` layers pack four elements a `u32`: a thread a word, the words at the range's
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
   again) sampled again, a device loss drawing the same image on the new device. Measured against step 1: on all units
   2.75 to 2.13 s a frame (per group 2.04), on one unit 18.8 to 11.3 ms (per group 12.9); the splat the rest.
4. **Then decide on splatting on the GPU** too: the accumulate kernel (`color-smoothing-accumulate`, today on the
   texture meshes' vertex arena) reading merged vertex buffers. Per-group CPU meshes would gain from it as well.

## Follow-ups

- **The half-cell shift**: the mesh shaders (`texture3dFrom2dLinear`) sample a grid half a cell from where the splat
  put its values: the splat writes a cell's value at its corner (`p = i`), the filter reads texel `i` at its center
  (`p = i + 0.5`). Per-group smoothed layers and colors have it, and so the GPU sampling (the CPU path sampled at the
  corners). Check WebGL's `texture3dFrom2dLinear` too; fix in the shaders (and the kernel) by sampling at `p + 0.5`.
- **The kernel's pipeline** is compiled on first use (`createKernel`, synchronously): in no precompile list.
- **Positions and transforms** uploaded a layer (12 bytes a vertex beside the vertex buffer): the kernel could read the
  merged object's vertex buffer and instance data instead.
