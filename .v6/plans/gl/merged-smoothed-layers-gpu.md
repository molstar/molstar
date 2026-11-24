# Smoothed layers sampled on the GPU (merged surfaces)

Status: planned. Follows [merged-units-core.md](merged-units-core.md) step 4b.

## Why

A smoothed layer (overpaint, transparency, emissive, substance) on a merged surface is sampled per vertex and instance
on the CPU, where per-group visuals sample their grid in the fragment shader. Measured (60 copies of 7qpd as a root
structure: 7 groups of 60 instances, CPU meshes, WebGPU), the layer on all units: 1.50 s a frame merged against 1.03 s
per group, the same 35 MB of layer data. These are medians of frames that set and cleared the layer in turn, so they
understate the cost of setting it (step 1 measures again).

- **Shared by both paths**: the splat into the grid, on the CPU (about 1 s).
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

1. **Baseline**: a benchmark case `merging/assembly-gaussian-smooth-layer-{one,all}` (merged) that sets one of two
   layers in turn, so that every frame smooths; measured; the numbers in merged-units-core.md and the changelog entry
   corrected.
2. **The kernel**, dispatched for one segment into a plain buffer. Tested on WebGPU (headless) against
   `getTrilinearlyInterpolated`: the bytes within a small tolerance (the filter's precision against doubles).
3. **Merged objects**: GPU-written ranges, skipped by the packer, dispatched again after layout changes. Tests: images
   against the per-group visuals' (the same filter), a segment added and removed keeping the layer, a backend switch,
   the geometry export reading the drawn layer back.
   Measured against step 1.
4. **Then decide on splatting on the GPU** too: the accumulate kernel (`color-smoothing-accumulate`, today on the
   texture meshes' vertex arena) reading merged vertex buffers. Per-group CPU meshes would gain from it as well.
