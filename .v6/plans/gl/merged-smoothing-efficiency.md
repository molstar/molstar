# Merged smoothing: what is left of its cost

Status: item 1 measured and dropped, item 2 done; items 3 and 4 open. Follows [merged-smoothed-colors-gpu.md](merged-smoothed-colors-gpu.md) and
[merged-smoothed-layers-gpu.md](merged-smoothed-layers-gpu.md): the items of the efficiency review of that work not done
yet. Done from it: a layer written once a geometry (`layersWritten`), a rebuild keeping the geometries, the regroup only
where a layer may smooth, the layers' kernels alone prepared, a layer's input reused (`DataArray.spare`); not done: a
split's objects made again (measured: no time, `merging/assembly-gaussian-smooth-layer-move` 2.89 against 3.08 ms).

## Items

### 1. GPU objects made and destroyed by every write (dropped: measured)

- **What**: a splat (`smooth`, `gl/compute/color-smoothing.ts`) makes two kernels and three buffers and destroys them; a
  sampling write (`computeLayerSampling`) a kernel and its own compute pass; a layer's data (`createGpuLayerData`) its
  transforms and reached buffers.
- **Measured** (a CPU profile of `merging/assembly-gaussian-smooth-layer-one`, WebGPU, 48 updates): a frame 3.3 ms, of it
  the update 0.74 ms on the CPU and the GPU's 1.81 ms; under `setState` the kernels' and buffers' creation, destruction
  and bind groups about 0.16 ms an update (5% of the frame), the input's upload (`load`) as much, the reached instances
  (`getSplatInstances`) 0.09 ms. Not worth the state that keeping them takes (kernels a writer, transforms by their
  cell, scratch sums).

### 2. A trajectory frame writes the layers twice (done)

- **What**: a new structure (a frame) remakes the geometry; the merged visual smooths its layers again (`resmooth`)
  with its state's loci, of the previous structure; the layers' nodes then set them remapped to the new one
  (`SyncStructureRepresentation3DState`), which are not equal to what was written (`Loci.areEqual`: another structure),
  so written again (`layersWritten` skips only equal ones).
- **Design**: the update smooths the layers remapped to its structure (`Overpaint.remap`, ..., as `setState` remaps
  them) and keeps those as written and as its state: the nodes' set then equal, skipped.
- **Test**: a new structure and geometry with a smoothed layer, then the layer set as a node sets it (its loci of the
  new structure): its cells not written again (without the remap: written twice). **Measured**: the saving one layer
  write a frame (with a layer on all units of 60 copies of 7qpd about 33 ms on WebGPU, seconds on the CPU); no trajectory
  case (no benchmark structure has frames, an alternating structure would measure its surfaces' remaking).

### 3. The grids span every instance

- **What**: a smoothed layer's or color's grid covers the box of all instances (the per-group shaders read it without
  a bounds check), its texture kept by the segment's writer (`createGpuLayerData`) for its writes again; the normalize
  writes zeros to every cell outside the splat's region, into a new texture already zero. Measured: the merged colors'
  textures 99 to 132 MB (60 copies of 7qpd), a layer on one of 60 copies a grid for all.
- **Design**: where sampled per vertex on the GPU (`perVertexOnGpu`: only `layer-sample-vertices` reads it), the grid
  sized to the splat's region (the reached instances and their neighbours), the kernel reading 0 outside it (a bounds
  check, as the CPU's interpolation); the normalize not writing outside the region of a new texture. The per-group
  grids as they are (their shaders unchanged).
- **Test**: the sampling with a grid smaller than the instances' box: the instances outside it zero, those inside as
  before; the bytes of the textures (`getByteCount`) of a layer on one unit, before and after.

### 4. Data per vertex for every instance

- **What**: a smoothed layer's data per vertex and instance is `vertexCount × instanceCount` elements a segment, the
  unreached instances zeros in the data pages (a layer on one of 60 copies: a segment's 60 instances stored, one used,
  about 5 MB of 60 copies of 7qpd).
- **Design**: the data of the reached instances alone, through a base a instance (none: no layer) read by the merged
  shaders where they read `vertexInstance` layers; the colors keep every instance (all reached).
- **Why last**: the merged shaders' layer reads change (WGSL and GLSL, the records' tables), for memory alone; done
  where the memory measured matters (a layer on few copies of many).

## Steps

1. **Item 1** (dropped): measured, 5% of the frame.
2. **Item 2** (done): the remap in the update.
3. **Item 3**: the region-sized grid for the per-vertex sampling; the bytes measured.
4. **Item 4**: decided by the memory items 3 leaves (a layer on few copies of many measured), then designed in detail.
