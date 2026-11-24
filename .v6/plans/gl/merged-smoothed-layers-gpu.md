# Smoothed layers sampled on the GPU (merged surfaces)

Status: planned. Follows [merged-units-core.md](merged-units-core.md) step 4b.

## Why

A smoothed layer (overpaint, transparency, emissive, substance) on a merged surface is sampled per vertex and instance
on the CPU, where per-group visuals sample their grid in the fragment shader. Measured (60 copies of 7qpd as a root
structure: 7 groups of 60 instances, CPU meshes at resolution 1, WebGPU; `merging/assembly-gaussian-smooth-layer-*`,
the layer set each frame): on all units 2.75 s a frame merged against 2.04 s per group, the same 35 MB of layer data; on
one unit 18.8 against 12.9 ms.

- **Shared by both paths**: the splat into the grid, on the CPU (about 2 s).
- **The merged path's own**: the trilinear sampling at every vertex of every reached instance (8.9 M samples), and the
  35 MB written to the CPU arrays and uploaded.

## Today

`Visual.setOverpaint` (and the others) with a smoothing context and no graphics context: `applyMeshOverpaintSmoothing`,
then `calcMeshColorSmoothing` splats into a grid and samples it (`getTrilinearlyInterpolated`), giving the record's
`tOverpaint` as `vertexInstance` data. The merged object packs the records' arrays into its layer's range of the
storage pages (`DataLayers`, written by range: `TextureChanges`), a segment's data at its `vertexInstance` offset (the
sixth `SegmentTexels` element). The averages (`transparencyAverage`, `transparencyMin`, `emissiveAverage`) come from
the group-instance data before smoothing.

## Design

- **WebGPU with compute only.** WebGL and contexts without compute keep the CPU path. Images as the CPU path's within
  the rounding (f32 against doubles, rgba8).
- **The grid as now**: splatted on the CPU, uploaded as a texture (as the per-group grids, `kind: 'volume'`). The
  segment holds it (not a schema value: `tOverpaintGrid` is a shared cell) until the layer is cleared or written again.
- **A kernel samples the vertices** (`layer-sample-vertices`), per reached instance and vertex:
  - its position from the merged object's vertex buffer (the segment's range);
  - its transform from the instance data;
  - the grid sampled by the texture's linear filter, as `texture3dFrom2dLinear` (the per-group shader's): images as
    the per-group visuals', not bytes as the CPU path's (filtering weights are of the hardware's precision);
  - the result written into the layer's storage page at the segment's `vertexInstance` offset. `r8` layers pack four
    elements a `u32`: a thread a word, not a vertex, so that no two threads write one word.
- **The GPU owns those ranges.** The record's layer marked GPU-written: its CPU array is the group-instance input,
  not the drawn data, and the merged object's packer skips its range. A layout change (segments moved, a part's object
  made again, a page reallocated) dispatches again for the moved segments, from the grids they hold.
- **CPU readers read back.** The geometry export (through `getMergedMembers`) reads a GPU-written range back, so it
  exports the layer as drawn; the export is asynchronous already, so a readback before it reads the members.
- **A new context** (a backend switch, a restore): written again through the update that makes the parts again
  (`writeLayers` on new geometry).
- **Memory**: the per-vertex data on the GPU only. The CPU arrays stay at group-instance size, not the 35 MB.

## Steps

1. **Baseline** (done): `merging/assembly-gaussian-smooth-layer-{one,all}` (merged) set the layer in one of two colors
   in turn, so that every frame smooths; measured above. They and `merging/assembly-gaussian-smooth-merged` pin
   resolution 1: `auto` quality gives the assembly 3, which nothing smooths at.
2. **The kernel** (done): `layer-sample-vertices` (`gl/compute/layer-sampling.ts`), dispatched for one mesh into a
   plain buffer, from offsets in its positions, transforms and output. Kernels may sample a texture filtered: a linear
   texture of a format every device filters gets a filterable binding and the kernel the sampler `sLinear` (as render
   shaders). Tested on WebGPU (headless) against a CPU model of the mesh shaders' sampling (`texture3dFrom2dLinear`):
   within 0.76 of a byte step for `rgba8` and `r8`, an `r8` range at an unaligned offset keeping its neighbors' bytes,
   an unreached instance zero. Not against `getTrilinearlyInterpolated`: it samples at the cells (the splat's points),
   the filter at the texels' centers, half a cell apart, as the per-group shader.
3. **Merged objects**: GPU-written ranges, skipped by the packer, dispatched again after layout changes. Tests: images
   against the per-group visuals' (the same filter), a segment added and removed keeping the layer, a backend switch,
   the geometry export reading the drawn layer back.
   Measured against step 1.
4. **Then decide on splatting on the GPU** too: the accumulate kernel (`color-smoothing-accumulate`, today on the
   texture meshes' vertex arena) reading merged vertex buffers. Per-group CPU meshes would gain from it as well.
