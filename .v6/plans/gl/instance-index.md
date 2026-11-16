# Plan: The instance store and the instance index (the instance grid split)

## Status

- **Steps 1 to 3 done** (As built; step 3 for the GPU-culled merged objects alone: the others keep the grid's order,
  decided); step 4's pose changes done (As built), its adds and removals open (after the editable lists' plan); steps 5 and 6
  (WebGL) deferred (decided; the measurement and the options under Steps). Dynamic, editable
  particle lists: their storage and changes in [../particle-slots.md](../particle-slots.md) (slots in pages, a list's
  changes), what this plan's adds and removals take; what edits them (a tool, its state, undo) is not this one's. Did [load-time.md](load-time.md) step 7b-d (the grids on demand); its step 7a (the culler's identity region in
  an order) is this plan's WebGPU draw order.

## Why

- **The instance grid (`core/math/geometry/instance-grid.ts`) does five things, one of them an index**:
  1. **The instances' order**: built, it reorders the data, `aTransform`, `aInstance` and the poses' `order`
     (`transform-data.ts` `updateTransformData`, `updatePosesData`), a cell's instances a contiguous range; every
     reader of the data reads grid order, culled or not.
  2. **The cells** (`cellOffsets`, `cellSpheres`): the CPU culling's frustum, occlusion, clip and level tests
     (`renderable.ts` `cullSegment`), its draws a range a cell (base instances; WebGL without them attribute offsets).
  3. **The batches** (the top grid): `cullSegment`'s alone (`getCullGrid`, `getGridClipClasses`).
  4. **`cellTransform`**: the matrices copied in grid order (1's storage), read as `aTransform` everywhere.
  5. **Besides**: `cellSize <= 1` the flag of no levels of detail (`merged.ts` `writeSegmentLod`, the merged culler's
     records and segments' data, `markedSegment`, the sphere raster's levels); `centers` the bounding sphere's seed;
     `layout` and `placement` the moves (`moveInstanceGrid`, the culler's `isLayoutMoved`).
- **No dynamic updates but a whole move**: with the order in the data, an instance added or removed, or a pose
  changed across cells, builds the grid again, writes all the data again and lays the merged object out again
  (`rebuildInstances`); only a move of all instances alike keeps it (`moveTransformData`, the particle targets alone).
- **WebGPU's culled objects use none of 2-4 a frame** (the merged culler bins its own scene cells,
  `scene-cells.ts`): they pay the build at load (the yeast cell: `calcInstanceGrid` 1.6 s of the 6.3 s left) and keep
  their data in grid order for the frames the CPU culls (a view of its own: the image pass, illumination's thickness
  pass; the kernels compiling).
- **The levels hang on the grid**: an object has levels of detail where it has a grid, not by a property of its own.
- **Suspected**: the position iterators (`mesh.ts` `createPositionIterator`, and cylinders', lines', text's,
  texture-mesh's, direct-volume's, `center-group.ts`, `sphere-sites.ts`) read transform row `k` for instance `k`,
  where row `k` holds instance `aInstance[k]`: per vertex and instance colors and sizes
  (`createVertexInstanceColor`, `createVertexInstanceSize`: a slot by the instance, read by `aInstance` in the
  shaders) at another instance's positions wherever a grid reorders; the geo exporters (`extensions/geo-export`) pair
  rows and instances the same way. Not verified (step 1).

## Design

1. **The instance store**: an object's instances in a stable order (their ids), never reordered by an index:
   transforms or poses, `aInstance` the identity. Added: appended. Removed: the last moved into the gap (decided: a
   swap, compact; its id kept by `aInstance`, so picks and per-instance values follow it). A pose changed: in place, the changed ranges uploaded (the GPU's
   expansion writes slots by id, `culled-poses.ts`). The bounds kept with it: a box grown by adds and moves, computed
   again after removals (deferred: a conservative box meanwhile).
2. **The instance index**: derived from the store, versioned with it: cells of instance indices (a list a cell, not a
   range of the data), their spheres, the batches over them (the CPU culling's, as now); built where a CPU cull needs
   it, not at load. Draws read the instances by index:
   - **WebGPU**: culled objects already do (`sCulledIndices`; the CPU culling's lists ranges of the identity region in
     the index's order, `setIdentityOrder`). Objects culled on the CPU keep their data in the grid's order (decided:
     their grid is needed at the first frame anyway, small ones the most; an update rebuilds their grid and reorders
     their data), no read by index a vertex.
   - **WebGL2**: merged objects read the instances' data from texture by index already (`tInstanceData`,
     `uInstanceBase`): a draw slot's instance by one more fetch.
   - **WebGL1**: the data copied in the index's order (as now), made where drawn.
3. **The levels of detail by a property of their own**: whether an object has levels and the cell size their scales
   take, not `cellSize <= 1`.
4. **Updates, local and scene-wide**:
   - **An object's index**: a pose changed within its cell grows the cell's sphere (refit later), across cells moves
     the index between the lists; an add or a removal one list; the lists kept with slack or counted again (O(n): a
     million instances in milliseconds). Changes batched a frame.
   - **WebGPU scene-wide**: the merged culler's scene cells built again on the GPU a batch of changes (milliseconds,
     `load-time.md` step 4); no incremental cells.
   - **WebGL scene-wide**: one CPU index over all of a merged object's instances ([particle-pack.md](../particle-pack.md)
     step 7's grid), lists of indices into the store, in place of the members' grids; the pack's Morton order then an
     optimization (coherent lists), not what the ranges need.
5. **Moves**: a move as a whole moves the store's transforms or the poses' matrix and the index's spheres (as
   `moveInstanceGrid`), no index built; the culler's `isLayoutMoved` by the index's layout.

## Steps

1. The suspected position-iterator pairing (Why): a test of per vertex and instance colors on an instanced mesh whose
   grid reorders (a volume or a positional theme), and of an export; fixed where it holds.
2. The levels of detail by a property (item 3), the `cellSize <= 1` reads replaced.
3. The store and the index on WebGPU (items 1, 2): the data in id order, draws by index, the index built where a CPU
   cull needs it (`load-time.md` step 7's gain: the grids' build off the load).
4. The updates (item 4), designed for both, pose changes first (decided), then adds and removals: a benchmark first
   (a particle list's poses changed a frame, a part of them and all: what a frame costs today, where), then pose
   changes in place (the store's changed ranges, the culled data's slots, the scene cells built again on the GPU, a
   CPU index refit), then adds and removals (a swap), the merged object's segments without a re-merge.
5. WebGL2 by `tInstanceData`, WebGL1 the copy in the index's order.
6. The scene-wide CPU index for WebGL (item 4).

Steps 5 and 6 deferred (decided). Measured: the yeast cell rotating on WebGL2 (no GPU culler: the members' grids, about
20 M cells, built at creation) about 180 ms a frame, all `cullSegment`: about 150 ms the cells' tests, about 35 ms the
lists' emits, an entry a member's run of a cell, about one a visible instance for members spread over the scene. The
options weighed: the data copied in scene-cell order (attributes and multi-draw as now: the tests about 1 M cells, the
emits as they are: est. 180 -> 45-60 ms); a compaction on the CPU as the GPU culler's (the scene cells tested, a
visible instance's row written into its member's region a level, an instanced order attribute, the data read from
the instance data texture by row, a multi-draw a page of 4 M instances, ANGLE's ~1.2 GB textures: a draw a member and
level; the most work and gain); the members' grids coarser on WebGL (as where the GPU culls; MG_1184 culled on the CPU
measured mixed before: zoomed 6.3 -> 18.6 ms, whole 15.7 -> 5.2).

Each step: the culling tests (`cull-segment`, `webgpu-gpu-culling`, `transform-data`, the particle targets'), the
visual suites (geometry, renderer), MG_1184 and the yeast cell's load and frames.

## As built

- **Step 1**: confirmed and fixed: the position iterators (mesh, cylinders, lines, text, texture-mesh, direct-volume,
  `center-group.ts`, `sphere-sites.ts`) read an instance's row (`getInstanceRows`: the inverse of `aInstance`, kept
  per its version, none in the instances' order), the geo exports their matrices in the instances' order
  (`getInstanceMatricesById`; color smoothing keeps `getInstanceMatrices`: it pairs rows with `aInstance`). Also in
  5.13 (the root `CHANGELOG.md`). Tested: `instance-positions.test.ts` (a mesh by matrices and by poses, lines; the
  exports' matrices), failing before the fix.
- **Step 2**: the transform data's `instanceCellSize` and `instanceBatchSize` (`BaseSchema`), set by
  `updateTransformData` (what the grid is built with, and what a rebuild takes, `Visual.setTransform`);
  `hasInstanceCells` (`transform-parts.ts`: a cell size above 1) in place of the grid's `cellSize` reads
  (`cullSegment`, `markedSegment`, the renderable's `cull` and marked lists, the merged object's segments' levels and
  its culler's records and segments' data, the sphere raster's levels; `getCullFacts` by `instanceCellSize > 0`,
  watched); a merged object's own 2 (its members' cells culled). The grid's `cellSize` is now only what it was built
  with. Tested: `transform-data.test.ts` (the sizes as created and rebuilt, matrices and poses), `cull-segment.test.ts`.
- **Step 3, GPU-culled merged objects**: a grid pending (`InstanceGrid.pending`, `createPendingInstanceGrid`: no
  cells, the cells' invariant sphere kept), asked by `createPosesTransform`'s `pendingGrid` (the particle targets
  culled on the GPU, `createTargetParticlePoses`) and kept pending through rebuilds (`updatePosesData`) and moves
  (`moveTransformData` moves its placement: the culler's moves as before); the instances in their own order (`aInstance`
  the identity, the poses without `order`). The merged object (`merged.ts`) builds its pending members' grids at its
  first CPU cull (`ensureIndexes`, by `cullVersion`: `buildPendingGrid`, `calcPendingInstanceGrid`, the data not
  reordered) into cells of its own for `cullSegment`, and writes their orders into the culler's identity region
  (`setIdentityOrder`); a change of the members drops them (`syncIndexes`: the region in the instances' order again);
  the marked lists' draw slots the built order's (else the identity). Without a reader by index the pending grid is
  built and applied as before (`applyPendingGrid`: an object not merged at its renderable's creation and update, a
  merged one without a culler). A large target's bounding sphere seeded from the pack's bounds (step 6: the box's center,
  the distance plus the invariant center's offset), no pass. No CPU cull at load: the frames before the culling's
  kernels hold already (`load-time.md` 7d not needed). The yeast cell: `calcInstanceGrid` 1.6 s -> 0, `createRO`
  2,460 -> 1,007 ms, the CPU to the first frame 6.3 -> 4.7 s, add and commit 0.6-0.9 s, the frames as before (4.0 ms).
  Tested: `webgpu-gpu-culling.test.ts` (Dawn: pending members drawn as built ones on the GPU and in a CPU-culled view,
  their data untouched; unmerged ones applied), `target-common.test.ts` (a GPU-culled target's grid pending, no pass).
- **Step 4, pose changes**: measured first: `particles/pose-frames-{all,part}` (the benchmark's on-demand cases: the
  PetWorld pack culled on the GPU, a new list a frame with every particle or a tenth moved): 52-57 ms of update a frame,
  a tenth moved as costly as all (everything built again: the colors, the culled data, the culler's layout, the targets'
  bounds and spheres). A pose change is the list's lineage (decided): `ParticleList.posesFrom` (the first list of a line
  of copies with other poses, as `movedFrom` records moves), made by `Particle.withPoses` (the coordinates, the
  rotations; what is cached by position not kept), read by `Particle.isSamePoses`. The particle target treats such a
  list as a move does (`posed`): its transforms built again, its colors and sizes kept unless positional, its pending
  grid and cell size kept (`createTargetParticleTransform`'s `posed`; `lodLevels` bumped only where the grid changed);
  the transform data keeps a pending grid and `aInstance` where they hold (`updatePosesData`; and a bug fixed: a
  pending grid asked for again was dropped, `createPosesTransform`); the culler lays out nothing where the members'
  pending grids are those split (`isPosesChange`, `posesLayout`: the cells built again on the GPU, the poses' slots
  written again, not the regions, the records, the words, `aInstance`, the segments; `MergedCullerCounts.poses`). A
  frame's update 52-57 -> 24 ms; what is left (per frame): the poses' values uploaded and expanded (all of them, about
  8 ms), a few instances' spheres (the boundary helper, 7.5), the targets' bounds (3). Tested:
  `particle-list.test.ts` (`withPoses`: one line, the caches by position its own), `target-common.test.ts` (a posed
  target's pending grid and cell size kept), `webgpu-gpu-culling.test.ts` (Dawn: pending members posed otherwise, no
  layout, drawn as built at the new poses).

## Open questions

- Changed ranges with the lineage (only the moved particles' values uploaded, their targets' spheres and bounds):
  for partial edits, with the editable lists' plan?
