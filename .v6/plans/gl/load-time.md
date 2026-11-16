# Plan: Loading large merged scenes faster (the yeast cell)

## Status

- **Steps 1 to 6 done** (As built; step 3 as the grid's centers shared, not the grids deferred; step 4 the cells
  dense, not compacted; step 6 the targets' bounds alone, not the Morton order); step 7a done, 7b-d done as
  [instance-index.md](instance-index.md) step 3 (the members' grids pending, built at the first CPU cull). Measured (Why). Takes the open part of [particle-pack.md](../particle-pack.md) step 7 (the scene-wide
  instance grid) for what loading needs; the rest of step 7 (the CPU culling's grid) stays there.
- **Stopped here** (the follow-ups under As built done). Where it stands, the yeast cell in Prism (WebGPU, Chrome
  headless, the pack served locally; ms from navigation): the pack fetched by 0.81 s, the cell first drawn at 2.47 s
  (fading in), whole at 2.94 s, its final frame 3.34 s. Its 14.2 s of CPU to the first frame (Why) now about 1.7 s from
  the data to the first frame, the main thread's alone: the pipelines compile meanwhile (precompiled from 1.44 s, their
  promises resolved as the thread frees, ~2.48 s), so only less main-thread work moves the first frame.
- **What is left of those 1.7 s** (ms; none measured further): the commit's uploads, the poses' values as stored
  (503 MB) and the sites (163 MB), ~275 (data: only a smaller encoding helps, e.g. rotations in 6 B); the rest of the
  commit ~150 (`read` 42, `load` 23, the culler's layout: not looked into); the response's body ~190 (712 MB read
  whole: a parse while streaming would overlap it); the hierarchy theme's colors ~160 (`computeDistinctColors`,
  sampling in Lab with RGB round trips: a table of the 8-bit RGB's linear values gives ~23 ms with the same colors, more
  only with others, which sessions would see); the merging ~120 (`createMergedRenderObject` 80, `getMergeKey` 38: loops
  a member, 5,639); the target representations ~260, thin (`createValues` 30, `getSitesMaxRadius` 28,
  `getInstanceRows` 22). Also left: the markers' zeros uploaded (31 MB, ~20 ms with their allocation: lazy markers
  touch marking and picking; not worth it alone).

## Why

- Measured: the yeast cell (31.5 M instances, 5,639 targets, `merging/yeast-pack-uniform` as the benchmark builds it),
  WebGPU, RTX 4070 Ti, one CPU profile of the load to the first frame: 14.2 s of CPU (19.3 s sampled, 5.1 s idle). The
  pack's parse and the particle list take milliseconds; the time is the GPU-side structures. Inclusive, ms:

  | Phase                                                                                            | ms    |
  | ------------------------------------------------------------------------------------------------ | ----- |
  | render objects (`createRO`)                                                                      | 4,952 |
  | — the instances' bounding sphere (`calculateInstancesBoundingSphere`, twice)                    | 2,305 |
  | — the targets' instance grids (`calcInstanceGrid`: bottom 1,391, top 143)                        | 1,615 |
  | the merged culler's creation (`createMergedCuller`: the first layout)                            | 3,281 |
  | the first dispatch (`dispatch`: the second layout and the instance data again)                   | 2,967 |
  | — the scene cells (`buildSceneCells`, both layouts; `writePoints` 1.9 s of it)                   | 3,816 |
  | — the instance data (`writeCulledData`: `writePoses` 904, `writeBuffer` 1,085; both layouts)     | 2,119 |

  A load in Prism (WebGPU) shows the same functions on top (`writePoints` 24.5% self, `buildSceneCells` 27.7% total,
  `writePoses`, `writeBuffer`, `getMaxDistanceSq`, the grids' `_build` and `calcBottomGrid`).
- **The culler lays out twice**: its first dispatch with the sphere raster turns on the raster regions
  (`dispatch`: `rasterRegions = true`, then `sync`), part of the layout key (`getFixedKey`), so the scene cells are
  built again (`syncLayout`) and every instance's data written again (`syncInstances`), though only the regions and
  records change, not the cells. The clip ghost's region the same (`ghostRegion`).
- **The bounding sphere twice**: `createValues` and `updateBoundingSphere` each compute the instances' sphere
  (`calculateBoxTransformBoundingSphere`: `writePoints` and `getMaxDistanceSq` over every instance).
- **The targets' grids where the GPU culls**: coarse already (`getTargetCellSize`, about 32 instances a cell), they serve
  only frames culled on the CPU (kernels compiling, a fallback), yet are built for all 5,639 targets at load.
- **The scene cells on the CPU**: each instance's center by its transform (`writePoints`), binned into a uniform grid
  per size class, a cell's sphere from its instances, uploaded (`cellSpheres`, `instanceCells`): the CPU never reads
  them back, the GPU has the instances' poses already.
- **The instance data converted on the CPU**: the poses written as 32 B a GPU instance from the list's rotations and
  positions (`writePoses`), about 1 GB uploaded (`writeBuffer`).

## Design

1. **One layout**: the regions (raster, ghost) apart from the cells: turning one on lays out the records and regions
   again (a new `setRegions`), not the cells nor the instance data; or the raster regions known at creation where the object
   can be drawn by the raster (`isSphereRasterEligible`, `sphereRaster` on). The first dispatch then reuses the cells.
   Est. -1.9 s (the second `buildSceneCells`) and -1.0 s (the second instance data).
2. **The bounding sphere once**: the instances' sphere computed by `createValues` kept for `updateBoundingSphere` where
   the transforms' version is the same (a cached sphere by version, as `CenterGroupData` caches its own). Est. -1.1 s.
3. **The targets' grids' pass shared**: the grid's centers (it writes every instance's center and their box) give the
   instances' bounding sphere (the box's center, the largest distance from it), the first bounding sphere then without
   a pass of its own. Not the grids deferred until a CPU cull needs them (as planned first, frames waiting for the
   kernels): a grid sets its target's instance order (`aInstance`, the poses' `order`), so one built later reorders
   the instances and the merged object writes their data again (about 1.4 s at the first hover, whose marked lists cull
   on the CPU), and the merged culler reads a member's grid for its levels (`cellSize` <= 1 gridless: drawn at full
   detail at every level).
4. **The scene cells on the GPU**: a compute pass from the instances' poses (already uploaded): each instance's center
   and size class (its segment's invariant sphere), its cell (`atomicAdd` per cell: counts), an exclusive scan of the
   counts (offsets), a scatter (each instance's slot in its cell: `instanceCells`), the cells' spheres from their
   instances (atomic min and max of the boxes as ordered `u32`s, then the sphere). The grids' extents and sides from the
   scene's bounding box (item 2's sphere, the classes' sides from the segments' radii: small, on the CPU). No CPU arrays
   (`cellOf` 8 B an instance); a moved layout (`moveLayout`) a dispatch again or the spheres moved as now. Est. -1.9 s
   (the first `buildSceneCells`), the cells' upload gone. WebGL keeps the CPU's cells (no compute; it does not use
   them: its culling is the targets' grids'). As built the cells dense, not compacted (no scan): their count is the
   grids' (known on the CPU, no read back for the kernels' sizes), an empty cell a thread of the cell test each; the box
   by the members' bounds, not a pass over the centers; WebGL has no merged culler (compute only), so no CPU cells.
5. **The instance data expanded on the GPU**: the list's poses uploaded as stored (snorm16 rotations, float32
   positions and scales: about 16-20 B an instance) and expanded to the culler's 32 B by a kernel. Est. -0.9 s
   (`writePoses`), the upload about halved. As built at the sync, not the dispatch: the draws read the data before the
   culling's kernels are ready (the identity region), so the expansion is submitted there, its kernel created at once.
6. **Precomputed in the pack** (`particle-pack.md`): what holds for the scene as loaded, recomputed where it changes
   (a move, `Particle.transform`): per target its instances' bounding sphere (items 2, 4: no pass over the instances);
   the particles of a target in Morton order of their position (a scene cell's instances adjacent: the GPU binning's
   scatter coherent, the CPU grids' builds cheaper, step 7's runs). Not the cells themselves: they depend on the merge
   (which targets one object holds) and the size classes, the renderer's choices. A format version: older packs
   computed as now. As built the bounds alone, an optional field of a target's row (no version), computed at the first
   use where a pack lacks them; the Morton order not now (it changes the particles' keys, its gain unmeasured).
7. **The targets' grids on demand where the GPU culls** (after step 6: the yeast cell's `calcInstanceGrid` 1.6 s of
   the 6.3 s left): what deferring them failed on in step 3 is the instance order, and a GPU-culled object does not need
   one in its data. Every draw of it reads its instances through the culler's indices (`sCulledIndices`), the CPU
   culling's lists too (ranges of the identity region, `identityBase`, `shiftLists`), so the grid's order can live
   in the identity region alone (the instances' indices in grid order, 4 B an instance), the data in the list's order.
   - At load, a GPU-culled target (`isTargetCulledOnGpu`) builds no grid: its transform data in the particles' order
     (`aInstance` the identity, the poses without `order`), a deferred grid (no cells, the cell size by its bounds,
     `getTargetCellSize`: what the culler's levels read, distinct from gridless `cellSize <= 1`), its bounding sphere
     the box sphere from the pack's bounds (step 6) where at least 1,000 instances (the box sphere either way; a few
     instances' as now, by the boundary helper).
   - At the first CPU cull of the merged object (a view the GPU does not cull: the image pass's own frame, illumination's
     thickness pass, `createFrame` in `tracing.ts`; a ghost not dispatched), its members' grids built from their
     transforms (`calcInstanceGrid`, as now at load) and their orders written into the identity region (the culler's
     new `setIdentityOrder`); the CPU lists then as now. The marked lists' draw positions (`createInstancePositions`,
     the inverse of `aInstance` now) the inverse of that order. Moves (`moveTransformData`, the culler's
     `isLayoutMoved`) of a deferred grid: nothing to move but its placement.
   - Before the culling's kernels are ready (compiling at the first load), no CPU cull: the object not drawn, its
     kernels pending work the frames wait for (as a pipeline compiling, `holdFramesUntilPrepared`); else the first
     frames would build the grids at load again.
   - Not on WebGL (no GPU culler: its culling is the grids'), nor objects culled on the CPU, nor targets that are not a
     pack's range (no bounds: their box needs a pass anyway; as now).
   - What remains of the grid at load where the GPU culls: none; where a CPU cull happens, the build (about 1.4 s for
     the yeast cell) at that cull, once.

## Steps

1. One layout (item 1), measured; the visual cases with the raster and the ghost (`renderer/sphere-raster/merged-*`,
   the clip ghost cases) unchanged.
2. The bounding sphere once (item 2).
3. The grid's centers for the bounding sphere (item 3).
4. The scene cells on the GPU (item 4): the cells' contents against the CPU's (a test on Dawn: the same instances per
   cell, spheres within rounding), the culling's visual cases, the yeast cell's frames.
5. The instance data expanded on the GPU (item 5).
6. The pack's precomputed bounds and order (item 6), with the pack's writer (`cli`) and reader.
7. The grids on demand (item 7), in parts: (a) the culler's identity region from an order (`setIdentityOrder`), tested
   on Dawn (CPU lists over a permuted region draw what they drew); (b) the merged object's members' grids built at its
   first CPU cull, the order written, the marked positions from it; (c) GPU-culled targets without a grid at load
   (`createTargetParticlePoses`, the deferred grid, the bounds' sphere); (d) no CPU cull before the kernels. Each
   measured; the culling's tests (`webgpu-gpu-culling`, the image pass and illumination with a merged object), the
   Prism particle tests, the yeast cell's load and frames.

Each step measured as in Why (the load's CPU profile, the first frame), the frames after unchanged.

## As built

- **Step 1**: the culler's layout split in two, the cells (`syncLayout`: the scene cells, the instances' cells, the
  segments' data, the cells' spheres) and the regions (`layoutRegions`: the regions, the identity region's base and
  indices, the counters, the records); a region turned on at a dispatch (`setRegions`: the raster's, the ghost's) lays
  out the regions alone where the layout is otherwise current, its keys updated, not the cells nor the instances'
  data. The yeast cell's load (the CPU profile as in Why): the first dispatch 2,967 -> 76 ms, `buildSceneCells`
  3,816 -> 1,916 (once), the instances' data 2,121 -> 1,409, the regions 124; the sampled load 19.3 -> 13.4 s. Tested:
  `webgpu-gpu-culling.test.ts` (Dawn: the raster's first dispatch lays out nothing, `MergedCullerCounts.layouts`;
  before the change one more), the visual cases with the raster, merged objects and the ghost.
- **Step 2**: `calculateTransformDataBoundingSphere` (`gl/renderable/util.ts`): the instances' sphere of an object's
  transforms, the box sphere's center and largest distance (1,000 instances and more) kept per transforms' cell by its
  version, the poses' and the invariant sphere's center; a radius changed alone (an element target's sizes, the padding)
  adds to the distance. The spheres' `createValues` and `updateBoundingSphere` and the particle target's placement use
  it. The yeast cell: `updateBoundingSphere` 1,007 -> 409 ms (the rest targets under 1,000 instances: the boundary
  helper's sphere, its radius not additive), the bounding spheres together 2.1 -> 1.9 s (the first pass as before, within
  the profiles' noise): less than the 1.1 s estimated. Tested: `transform-bounding-sphere.test.ts` (a radius changed
  without a pass, the transforms changed with one, both equal to the uncached sphere).
- **Step 3**: the grid's pass shared: `calcBottomGrid` (`core/math/geometry/instance-grid.ts`) also keeps the sphere of
  the centers it wrote (`InstanceCenters`: the invariant center, the box's center, the largest distance, from the
  float32 centers: one pass over them, not the transforms), a moved grid none (`moveInstanceGrid`); the transform
  data seeds step 2's cache with it where a grid is built (`setTransformDataCenters`, poses and matrices). The yeast
  cell: the bounding spheres 1.9 -> 0.65 s (the rest the targets under 1,000 instances), `createRO` 4.5 -> 3.1 s, the
  grids' build as before. Tested: `transform-bounding-sphere.test.ts` (the sphere after a grid's build without a pass,
  within rounding of the uncached).
- **Step 4**: the scene cells built on the GPU (`gl/compute/scene-cells.ts`, `scene-cells.wgsl.ts`): `getSceneGrid` on
  the CPU (the box of the members' bounds, the side by the instances per cell, a member's class by its sphere's
  diameter: per segment, no pass over the instances), `scene-cell-bin` a thread per instance (its sphere from the
  culled data as `instanceSphere`, its cell, the cell's box of centers and largest radius by `atomicMax` on floats'
  ordered bits, the min inverted), `scene-cell-spheres` a thread per cell (an empty one a negative radius, which the
  cell test skips); built at the dispatch after a layout or a move, before the cell test, in the same pass. The cells
  dense, not compacted: the yeast cell 1,127,204 dense where 508,858 hold instances; their count known on the CPU (the
  buffers' and dispatches' sizes, no read back). A cell's radius is its instances' largest own radius, not the
  segment's largest scale's (tighter). A move builds the cells again from the moved instances (`moveLayout`: the
  segments' scales by each one's move), so members need not move alike; the records and segments' data kept. The
  yeast cell (the CPU profile as in Why): the culler's creation 3,443 -> 1,275 ms (`syncLayout` 2,039 -> 84, the
  instances' data the rest), add and commit 3.4 -> 1.4 s, the CPU to the first frame 9.4 -> 7.1 s; the first frame
  1,178 -> 1,218 ms (the build on the GPU within its noise), a still frame as before; storage +42 MB (the boxes, 28 B a
  dense cell, kept for a move's build; the dense cells' spheres and levels). Tested: `webgpu-scene-cells.test.ts` (Dawn:
  each instance in the cell of its center and class, each cell's sphere that of its instances' box and largest radius,
  the empty ones negative, built twice into the same buffers), `scene-cells.test.ts` (the grids: classes, sides, the
  cap), `webgpu-gpu-culling.test.ts` (the GPU culling's images as the CPU's, a move laying out nothing), the visual
  geometry and renderer cases on WebGPU.
- **Step 5**: the poses' slots expanded on the GPU (`gl/compute/culled-poses.ts`, `culled-poses.wgsl.ts`):
  `CulledData.layout` writes the rest (`poses: false`) and lists its pose parts (`poseSlots`, `slotFirst`); a part's
  values go up as stored (a range's slices of the list's arrays as they are, an index list's elements gathered: grid
  or float32 positions, snorm16 or float32 rotations, float32 scales or a uniform one), with a record a part (its values'
  words and kinds, the grid, the center, the similarity, its first slot and thread); `culled-poses` a thread a slot
  (its part by a binary search) writes what `writePoses` writes. At the culler's sync (`writeCulledData`), not its
  dispatch: the draws read the data before the culling's kernels are ready (the CPU culling's lists in the identity
  region, a pick), so it is complete before them; the kernel created at once (not prepared in the background: the first
  submit waits for its compile). Batches of at most 64 MB of values (the binding limit where lower), a part split
  across them, one buffer reused and each batch submitted (`ComputeContext.submit`, new): a write after a submit lands
  after it, the queue's writes going ahead of commands not yet submitted. The yeast cell: `writeCulledData` 1,173 ->
  781 ms (`writePoses` 432 gone; the values' upload 212 ms for about 500 MB, the slots' were 1 GB), add and commit 1.7-2.0
  -> 0.84 s, the first frame 1.4-1.7 -> 1.3-1.4 s, the frames as before; no storage kept (the values' buffer
  destroyed after its submit). Left: the other per-instance arrays the CPU writes (the slot words, `aInstance`, the
  segment indices: 464 ms of the 781, many small writes). Tested: `webgpu-culled-poses.test.ts` (Dawn: the data the
  CPU writes, every kind of positions, rotations and scales, ranges and indices, a similarity, a snorm16 slice at an
  odd pair; in one batch and in batches splitting parts), `webgpu-gpu-culling.test.ts` (members of poses drawn as of
  matrices), the visual geometry and renderer cases on WebGPU.
- **Step 6**: a target's bounds (`ParticlePackTarget.bounds`: its positions' box and largest distance from its center,
  rounded outward to 0.001 Å) written by the pack's writer (`writeParticlePack`, from the positions as written) and
  kept by the reader on the list (`Particle.setTargetBounds`, by position: not kept by a move); a pack without them
  computes them for all targets at the first use (`Particle.getRangeBounds`, `computeRangeBounds`: on a grid in its
  integer steps, 231 ms for the yeast cell where decoding each position took 539). A target's poses take their box for
  the grid's cell size (`createTargetParticlePoses`: the box expanded by the largest scale times the invariant center's
  offset from the placed point), not a pass over the instances. Not the sphere: a few instances' is the boundary
  helper's (tighter than a box sphere, and what `calculateInstancesBoundingSphere` gives on the other paths, as a
  representation's transform: seeding it with the box sphere made two paths' bounds differ, the mesh target's 11 Å in
  `particles.test.ts`, and moved the smoothing grids of an instanced mesh); it is kept instead per transforms' version
  and invariant sphere (`getFewInstancesBoundingSphere`): the values' and the bounds' update compute it once. The
  yeast cell (re-converted): `createRO` 3,250 -> 2,460 ms (the box pass 334 gone, the few instances' spheres 684 ->
  350), the CPU to the first frame 7.3 -> 6.3 s; without bounds in the pack `createRO` about 2,690 (the bounds
  computed). The pack's JSON 2.7 -> 3.3 MB. Tested: `pack.test.ts` (the bounds written hold every particle, read
  without a pass, computed alike without them, a moved list's again, none for a range not a target's),
  `target-common.test.ts` (a pack's target sized without the box's pass, as from the instances' box),
  `transform-bounding-sphere.test.ts` (a few instances' sphere kept, computed again for another radius or transforms).
- **After the steps, the representations' creation** (a profile of the load after `instance-index.md` step 3: the
  sites geometry itself a view of the pack's bytes, 7 ms, `Spheres.fromSites`; the time what is derived from it):
  - **The sites' boundary** (713 -> 265 ms): a chain's box (`getFastBoundary`'s) from the quantized integers'
    extremes decoded once (`getSitesBox`, `sites.ts`; float32's rounding monotonic: the same box, the same sphere to
    the bit), not the sites decoded into arrays first; not stored in the pack (the helper's sphere has up to 196
    extrema a target).
  - **The colors** (420 -> 5 ms): a theme constant over runs of particles says so (`ColorTheme.particleRuns`: the
    entities' colors of a pack, a run a target); a particle target within one run takes one color (`getTargetTheme`:
    `uniform`), not a color a particle (31.5 M); the textures 376 -> 256 MB.
  - **The shared quad indices** (432 ms and 466 MB -> 24 ms and 39 MB): a merged object's draws of items (no elements)
    address their segment's by its first vertex (`drawIndexed`'s `baseVertex`, the first index 0; the culler's indexed
    records too, the sphere raster reading the base vertex), so the indices cover the largest segment's items
    (`drawCountMax`), not the merged geometry's 20.3 M; the culler indexed by that too. WebGPU; WebGL as before.
  The yeast cell: the representations 2.26 -> 1.04 s, the first frame 1.8 -> 1.07 s, the CPU to the first frame 5.5
  -> 3.5 s, the frames as before (4.1 ms). Tested: `target-structure.test.ts` (the sites' sphere equal to the decoded
  positions' fast boundaries), `particle-colors.test.ts` (a target of one run uniform, of two or of arrays a color a
  particle), the culling's and the sphere raster's tests on Dawn (merged impostors drawn by base vertex, direct and
  indirect), the visual geometry and renderer cases.

- **With the viewer's defaults** (Prism: the sizes by `particle-size`, colors by the hierarchy): an element target's
  sizes made first and its values with one uniform size, their largest (`createElementSizes`), not the representation's
  size a particle made (`groupInstance`: 31.5 M) and overwritten, nor the bounds computed twice for another padding;
  the hierarchy's and the compartments' colors by their runs too (`particleRuns`); the sites' largest radius by their
  types read in place (`getSitesMaxRadius`). The yeast cell so: the representations 3.04 -> 1.1-1.25 s, `createRO`
  2.02 -> 0.60 s (the sizes 792 -> 1.5 ms, the few instances' spheres 819 -> 356 ms). Left: those spheres (the
  boundary helper over the invariant sphere's extrema an instance, up to 196, for 2 to 14 instances: every
  representation's), the sites' boundary (246 ms), the culler's creation (611 ms).

- **A Prism profile**: (1) the uniform colors made every member of a merged object its first member's color (a merged
  texture's granularity is its first member's; uniform maps to none): a constant color now filled per instance
  (`createFilledInstanceColor`, `createTargetColors`; one evaluation, 38 ms for the yeast cell), not uniform. (2) The
  bounding spheres by the whole's size: a particle target of a list of 1 M particles or more takes the box sphere at
  any count (`instanceSphereFast`, `LargeParticleListCount`; `calculateTransformDataBoundingSphere`,
  `Visual.setTransform`), seeded from the pack's bounds: no pass (356 -> 1 ms). (3) An app's first frames (Prism's,
  not held for pending work as the benchmark's) culled the merged object on the CPU while the culling's kernels compiled,
  building every pending grid (`ensureIndexes`: 2.3 s, `calcBottomGrid`, `writePoints`): such a frame now draws
  nothing (`MergedCuller.areKernelsReady`; [instance-index.md](instance-index.md) step 3's part d after all). The yeast
  cell: the representations 1.04 -> 0.78 s, the CPU to the first frame 2.9 s; no grid built at load.
- **Distinct extrema**: a bounding sphere's extrema made distinct (`BoundaryHelper.getSphere`, 98 of them) by a sort by x
  and its index array, ~12 us a call: twice a target (its sites' boundary, its spheres' invariant sphere), 240 ms for
  the yeast cell's targets in Prism. Now one pass over a reused open-addressed table of the coordinates' bits, ~1.7 us
  a call.
- **Runs on the GPU**: of the culler's uploads (1.6 GB for the yeast cell's 31 M instances, 4 B an instance each) the
  identity region (written again where the regions grew: three times) and CulledData's slot words, `aInstance` and
  segment indices are runs, `(start + step i) | bits`: written by a kernel (`u32-runs.ts`, `CulledDataLayout.runs`; the CPU
  writes what is not one: a part's slot words in an order of its own, a segment's `aInstance` not serial). The culler's
  own uploads 758 -> 3 MB, `createMergedCuller` 586 -> 476 ms; left: the poses' values as stored (503 MB) and the scan of
  `aInstance` (41 ms). (2) A multi-sampled frame's jittered samples (an idle tick keeps the frame) culled with views of
  their own, a merged object's on the CPU (its GPU culling of the frame unjittered): its first such culling built every
  pending grid (1.7 s, at random after load). The samples now cull as the camera unjittered (`setJitter`,
  `Camera.cullView`): the frame's GPU lists.
- **Object colors**: a target's color constant over its particles was filled per instance (`createFilledInstanceColor`:
  merged, a uniform color is the first member's), 126 MB of the yeast cell's data page uploaded and 42 ms. Now a color
  type of its own, `object` (`ColorTypeIndex` 5): a texel, the texture's define; merged, the members' packed as a table
  and read by the segment (`aSegment.w`; the sphere raster's resolve sets it from the chunk's slot); the merge key holds
  the color type (`object` and `instance` share the define). The data page 321 -> 195 MB, all uploads 996 -> 744 MB.
  Left in it: the markers (31 MB, zeros at load) and the sites (163 MB).
- **Shared serial instances**: a target's `aInstance` in its own order was an array filled serially each (`fillSerial`,
  58 ms for the yeast cell's 31.5 M instances), a merged object's a copy of them all (126 MB), CulledData's runs a scan of
  it (`isSerial`, 41 ms). Now views of one serial array shared by all (`getSharedSerial`, core `util/array`; read only:
  replaced, never written; a grid's order its own array), a merged object's of serial members a view too
  (`rebuildInstances`), known serial without a scan (`isSharedSerial`). The CPU 58 + 41 ms and the copy's loop to the
  shared array's growth (25 ms, once), the JS heap after load 1.52 -> 1.32 GB.
- **Chains' extremes in the pack**: a sites target's boundary took each chain's box from a pass over its quantized
  positions (`getSitesBox`, 113-140 ms for the yeast cell's 20 M sites). A pack keeps each chain's extremes
  (`ParticlePackSites.extremes`, `ParticleSitesQuantized.extremes`; computed where a pack has none): the box decoded
  from them, the same to the bit.
- **Target signatures once per variant**: the precompile signatures of a particle target representation
  (`computeParticleTargetSignatures`) made a visual's values for each target (the yeast cell's 5,639: 302 ms). A
  target's signatures depend on its props, geometry choice, sizes and surface defines, not its data: the props are one
  object per kind and quality (`createParticleTargetPropsProvider` interns them), the signatures memoized by it and the
  rest, each list added once. 302 -> 12.5 ms; the first frame drawing the cell 2.86 -> 2.55 s after navigation (Prism,
  WebGPU, the pack fetched by 0.83 s).
- **Sites' defines in the signatures**: the precompile's signatures come from an empty geometry, so a sites target's
  impostors' `dSites` (quantized sites) and `dSiteColor` (the second color) were not among them: their programs
  compiled after the scene's commit (220-250 ms), the first frame waiting. A sites target declares them
  (`getSitesSignatureDefines`, as the geometry and the colors decide): compiled with the rest while the scene builds.
  The yeast cell first drawn 2.98 -> 2.47 s after navigation, its final frame 3.46 -> 3.34 s.

## Open questions

- Item 4, revisited after step 6: the scene cells compacted (a scan of the occupied cells, their count read back or the
  cell test dispatched indirectly), against the dense cells' empty threads (the yeast cell: 55% of the cells empty)
  and their storage?
- After step 5: the slot words, `aInstance` and the segment indices written on the GPU too (the segments' bases and
  the parts' orders uploaded, no per-instance array of the CPU's but `aInstance`'s), or in fewer, larger writes?
- The Morton order of a target's particles (item 6): the CPU grids' build (`calcBottomGrid`, 1.3 s) and the GPU
  binning more coherent, measured first (sorted at read in a test build); it changes the particles' keys.
