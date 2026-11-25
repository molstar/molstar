# Plan: Sphere impostors drawn from chunks of culled instances

## Status

- **Paused**: its reason was measured away. The writes that waited inside (8-40 times slower) did so because the
  benchmark ran frames uncapped (no vsync: the GPU's queue always full, 2.5-3.6 s of GPU work in a 4 s segment,
  mostly the sphere raster); paced by the display (`benchmark-particle-tiles.mjs --vsync`) the same path's writes are
  as fast inside as outside (up to 1 MiB 115 µs, uploads 40 ms of an orbit), frames at 60 Hz (p95 19 ms). The empty
  records went by the raster split in the paged prediction (below). Chunks still remove the records left (level 0's
  where the raster does not take a member whole) and make a pass's recording replay while moving; taken up again where
  that shows in a frame at a display's pace.
- **Step 1 measured** (Chrome, the yeast tiled at 16³ with levels, orbiting inside after the dolly, the records of both
  phases read back every 15 frames): of the levels' 22,564 records (4 × 5,641 members) 0-8 hold instances; the
  raster takes almost every instance. Those few are large molecules near the camera (4,000-300,000 items a frame in
  1-8 instances), so their chunks would be full (99-100 % at 64 lanes): no vertices wasted. The ~11,000 records
  encoded a frame are the prediction's (`predictRecords`: a member's sphere in the frustum and a level's range,
  inside the cell every member at the near levels), nearly all empty. The store's bound of chunks (`setChunkCapacity`)
  52-58 M (1.6 of 2 GiB), the frames' need at most about 4,700: a bound that safe cannot be allocated; capacity by
  read-backs with headroom (as the raster's) and a frame of the record path where they overflowed.
- **A cheaper part first**: the paged prediction leaves out no level the raster takes all of (the merged culler's
  `toRaster`: a member's largest sphere at a level's nearest distance within the raster's radius); with it the far
  levels' records go, level 0's stay (its range from 0: inside the cell every member's): about half.
- **The raster split in the paged prediction done**: a level's impostor record left out where its member's largest
  sphere (its raster row times the largest scale uploaded, `PagedCuller.maxScale`) is raster-sized at the level's
  nearest distance (`isRasterSize`, conservative: the member's sphere's nearest). Chrome, the yeast: frames inside much
  shorter (orbiting 205 → 744 frames in 4 s, p95 113 → 8.4 ms; the dolly 233 → 550, p95 43.7 → 16.8 ms), outside
  unchanged. The writes inside still slow (a write up to 1 MiB about 6.8 ms orbiting): not the draws alone, so the
  GPU-process explanation is at least incomplete; the GPU's compute passes (the culling, the raster) not in the timings
  measured so far (render passes only): measured next, before the chunks.

## Why

A tiled particle view of the yeast cell (5,640 members, four levels, two phases) encodes about 11,000 indirect draw
records a frame while the camera moves inside it ([particle-tiles.md](../particle-tiles.md) 9l): 0.9 s of a 4 s orbit
in `drawIndexedIndirect` and `multiDrawIndexedIndirect` on the main thread, and Chrome's GPU process busy decoding
them, so that every `writeBuffer` waits: 8-40 times slower than outside, 1.5 s of the orbit, with the GPU idle (3.5 ms
a frame) and no garbage collection (10a3). The records change every moving frame (`predictRecords`), so the pass is
recorded again every frame and draw lists save nothing. A record is a member at a level because a member's geometry
(its vertex count, first item and stride) is its own; the chunks carry exactly that per instance.

## What there is

- **The paged culler** (`merged-paged-culler.ts`, `paged-cull.ts`): per member and level a region of pose slots
  (`indices`), counted by `counts` from `bases`; per region a draw record (`args`: the member's vertex count at the
  level and first item from `layoutRecords`, its instance count and first instance from `paged_cull_args`). The
  impostor's vertex stage reads a slot's pose and member through `sCulledIndices[instance_index]` and `CulledData`,
  its item, sites and stride from `tSegment`, the level's range and index from the draw block (`uDraw.uLod`,
  `uLodLevel`). No other per-member uniform: members differ in `tSegment` and `CulledData` alone.
- **The sphere raster's chunks** (`sphere-raster.ts`, `sphere-raster.wgsl.ts`): per raster record `ceil(count / 64)`
  blocks, a prefix sum, then the regions kernel appends per instance its chunks of up to 64 items: two `vec4<u32>`,
  `(slot, first item, items, level | stride << 2)` and `(store index, group instance, group offset, level scale)`.
  Paged objects are raster inputs already (`getRasterInput`), but only their raster regions are chunked, into a
  shared scratch buffer sized from read-backs (chunks past it dropped for a frame).
- **The resolve** (`spheres.wgsl.ts`, `resolveSphere`): an impostor's vertex output made from a chunk and a lane,
  its level's range by `uRasterLevels[level]` and its scale from the chunk: the vertex stage a chunk draw needs,
  bound in the resolve's pipeline only.

## Design

1. **Chunks of the impostors' regions**: the blocks, prefix sum and regions kernels run over the levels' records too
   (and the raster's records where the raster does not draw a pass: `getLists(withRaster)`), per phase, into a
   buffer the object owns (not the raster's scratch: it lives from the cull to the frame's last pass, the depth,
   color, pick, marking and transparent passes alike). Each phase's chunks a range of its own, its count written into
   an indirect draw.
2. **One draw an object and phase**: a vertex stage variant (`dSphereChunks`): the chunk `instance_index`, the lane and
   the quad's corner from `vertex_index` (64 items × 6 vertices; lanes past the chunk's items moved out of the clip
   volume), the sphere read as `resolveSphere` reads it (shared code, not a copy), the level's range from an object
   table by the chunk's level (as `uRasterLevels`), not the draw block. Bindings: the chunks in place of
   `sCulledIndices` (the chunk holds the slot), `CulledData` as now (the budget of 8 storage buffers kept).
3. **The record path kept**: for objects whose chunks do not fit (the bound below), WebGL, and what the chunk draw
   does not cover at first (the ghost's pass: its coarsest level's uniforms; solid interiors' caps: their draw blocks
   bound around one draw each, a later step).
4. **Capacity without dropping**: impostors cannot miss a frame's chunks. The bound from the store (`setChunkCapacity`:
   per member and level its instances times its chunks) capped by the device's binding size; where the bound is past
   it, the record path. Measured on the yeast in step 1.
5. **What goes away for those objects**: `predictRecords`, `setListRecords`, `compactRecords` and the record lists:
   a pass's draws fixed (an object's two), so its recording replays while the camera moves.

## Risks

- **Wasted vertices**: a chunk draws 384 vertices however few items it has; the raster's chunks (far instances,
  coarse levels) average about 7 items in the yeast. The impostors' are near, their levels finer: likely fuller, but
  measured first (step 1). Smaller chunks for impostors (16 or 32 lanes) if not.
- **The second phase**: its regions after the first's (`prior`), its chunks after the first's chunks, its draw its own.
- **Clipping**: instance clipping by the chunk's segment (as the resolve), pixel and primitive clipping as now; the
  ghost's pass keeps what the clip objects remove (its chunks, or the record path).
- **Animation**: wiggle and tumble per vertex need the group and segment, both in the chunk; the raster's eligibility
  (opaque, still) does not apply.

## Steps

1. **Measure**: on the yeast tiled at 16³ with levels (the Chrome benchmark's path), per moving frame the impostors'
   chunks if they were made (the regions kernel run over the levels' records into a scratch, its counters read back):
   chunks, items per chunk, the bound against the binding size. Decide the lanes per chunk, or stop.
2. **The chunks of the levels' regions** into an owned buffer, per phase, with their indirect draws; tested by
   reading back: the chunks' slots and items equal the records' instances and items.
3. **The vertex variant and the draw**: one draw an object and phase; tested by images equal to the record path's
   (opaque, pick, marking, transparency, the second phase), the Viewer's visual cases on WebGPU.
4. **Measured** by `benchmark-particle-tiles.mjs`: draws a frame, the uploads' time inside while moving, frames.
5. **Then**: the ghost and the caps by chunks; `merged-culler.ts`'s objects (unpaged) the same way.
