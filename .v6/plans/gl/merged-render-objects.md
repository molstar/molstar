# Plan: Merged render objects, the path for many units

## Status

- Steps 1 (benchmarks, `auto`), 2 (the update costs), 3 (geometry without copies), 4 (layers without copies), 5 (the
  costs per member: location iterators, WebGL's segment draws, the merged `sync`, props kept as one object, the
  allocations of a recolor, the color themes' per-location work) and 6 (GPU culling, its occlusion test and `auto`;
  levels per instance and view-independent draw lists measured and dropped) done. Left: the visuals' own update work.
  The GPU culling's instance copy on weak GPUs (about 4 ms of a 55 ms PetWorld frame on an integrated GPU) measured
  and not taken up: a smaller instance transform would save about 2% and change `aTransform` in every program. Replaces batching
  ([batched-draws.md](batched-draws.md), dropped) for scenes of many objects: those come as many units of few
  representations, which CPU-side merging serves.

## Why

A representation of many units (assemblies, unit cells, many chains or ligands) makes a render object per unit and
visual. Merged (`renderMerged`), they are one render object: the renderer walks and culls one row, the GPU draws once
per level and pass. Measured (WebGPU, per frame): `units-1000-spacefill` merged 0.43 ms, `-unmerged` 2.4–2.8 ms (the
same with batching, 2.43 ms). Merging works on WebGL as well.

## Today

- **Where**: units representations (`units-representation.ts`) and particles, per geometry kind
  (`mergeRenderObjectsByType`): meshes, spheres, cylinders, lines, points (`MergedTypes`); not text, images, texture
  meshes, direct volumes.
- **When**: `renderMerged`: on, off, or `auto`, more than 30 units (`AutoRenderMergedThreshold`, 100 before step 1);
  members must share a
  merge key (`getMergeKey`: the defines, instance granularity, `uLod`, level distances) and supported data types
  (`canMergeValues`: no palette; color and size per uniform, instance, group, group instance or vertex; the layers per
  instance or group instance).
- **How** (`createMergedValues`): one render object whose vertex attributes, elements and theme layers are the
  members' arrays as parts (`SegmentedArray`: each at its vertex, element or texel base, indices biased by the vertex
  base at upload; a layer's image sized as a copy would be), instances concatenated with per-instance segment
  data (`aSegment`: group and data offsets; `aSegmentSphere`; `aSegmentLod`: the members' level scales); `segments`
  holds per member its instance, vertex, draw and element ranges. Members keep their own render objects (not drawn).
- **Updates** (`sync`, on the merged renderable's `update`): every member's value versions compared, a changed
  member's geometry and layer parts advanced (the render item writes those parts only), level scales copied,
  aggregates and bounds recomputed; a change of the layout (any
  member's sizes) rebuilds the whole. A new set of members (`key` of member ids) creates a new merged object.
- **Drawing**: culling and levels per segment on the CPU (`createSegmentedMdbList`, `cullSegment`): per level a list
  of draws (ranges and instance ranges); picking resolves the instance to member and local instance
  (`pickMergedLoci`).

## Questions to answer first

- Where merging pays: units-N merged against unmerged for N from 10 to 1000, per frame (CPU and GPU), setup, and
  updates (hover marking, a color theme change, one unit's geometry change such as a trajectory frame). Whether `auto`
  should merge from fewer units, or always where the key allows.
- What updates cost merged: the comparisons over all members and the copies per changed block; a hover over one unit of
  a merged object, a trajectory frame of all units (each member's geometry copied).
- Memory: members' data and the merged copies (twice the geometry and layers on the CPU, the copy on the GPU).
- Coverage: which common representations don't merge today and why (the merge key, the data types: e.g. volume
  colors, palettes), and which matter.
- Transparency: a merged object is sorted as one (blended transparency without order-independent passes).

## Design directions

- **Members' data as parts, not copies** (steps 3 and 4): the merged render item writes each member's arrays into its
  own GPU ranges at the member's base; no store needed, as members never upload on their own. Indices could stay
  relative (WebGPU `baseVertex` per segment draw) instead of biased at upload: saves
  the bias on updates, measured first.
- **Updates per member**: a member's change (marking, colors, geometry of a frame) updates its ranges and its segment
  data; a member's size change moves its ranges (the others stay); no full rebuild, no comparisons over all members
  where the change is known (the member's update tells the merged object).
- **When**: from the measurements, a threshold or always (where the key allows); membership changes (units added or
  removed by a selection) cheap with ranges.
- **GPU culling and levels of detail**: a merged object is one object with many instances; its block is read with a
  uniform index in every draw. A kernel tests the segments' instances (bounding spheres, levels from the distances)
  and writes per level the visible instances (a compacted list) and indirect draws; the per-instance data (transform,
  segment data) then read by the compacted index from storage instead of instanced attributes. Its cost measured first
  and properly (the same image as the CPU path; `draw-overhead/*`, merged objects of many instances).

## GPU culling (design, step 6)

For merged render objects on WebGPU with compute and `indirectInstanceRanges` (`BaseGeometry.Params.gpuCulling`,
`canCullOnGpu`); CPU culling (`cullSegment`) stays for WebGL and devices without them. A plain instanced object is
culled on the GPU as a merged object of one member (`mergeRenderObjectsByType`).

- **Cells, not instances**: a thread per cell of the members' instance grids (a member without a grid one cell) tests
  its sphere against the frustum and the object's distance range and picks its levels, as `cullSegment` does per cell;
  finer cells cost threads, not draws. Testing each instance cost more than it saved.
- **Compaction by copy**: a thread per instance of a drawn cell copies its instanced attributes (transform, instance
  id, segment data, about 112 bytes) into compacted buffers, a slot per level in its segment's region (`atomicAdd` on a
  counter per level and segment); the compacted buffers are the render item's instanced vertex buffers. The programs
  stay the same: no storage reads per vertex. Reading the instance data from storage through compacted indices instead
  (4 bytes per slot) costs every vertex an index and seven vec4 reads (measured: the opaque pass 12% slower for sphere
  impostors, worse for meshes).
- **Draws**: one indirect record per level and segment (`drawIndirect` or `drawIndexedIndirect`), its instance count
  written per frame, the other words once per layout. Not a record per cell: indirect draws replay whether they draw or
  not, and without multi-draw indirect with a count (not in WebGPU) the GPU can't merge the visible cells' ranges as the
  CPU culling does (measured: 51k cells, 21 ms per frame against 1.2 ms).
- **Fallback**: the compacted buffers end with all instances in their own order (the identity region); views not culled
  on the GPU (another view of the frame, a scaled model, the kernels compiling) draw the CPU culling's lists shifted to
  it, with the same programs and buffers.
- **When**: before the passes (`Renderer.prepareCull`, called by the draw pass after its `update`; passes update the
  renderer with their render pass open, a compute pass would end it), once per frame for its first view.
- **Draw lists**: the draw key of indirect lists is their buffers (the GPU writes the counts): recordings stay valid
  while the camera moves, were GPU-culled objects not view dependent (step 5 below, dropped).
- **Cost**: levels times segments indirect draws every frame, about 0.28 µs each (as the CPU culling's draws per
  segment); the copy about 112 bytes per drawn instance; memory about 112 bytes per instance per level plus one.
- **Occlusion**: per cell against the previous frame's Hi-Z levels (GPU textures, `textureLoad`, r32f), below.
- **When on** (`auto`, the default): from `AutoGpuCullingMinInstances` (20k) instances and `AutoGpuCullingThreshold`
  (30) per draw (segment and level). A GPU-culled object costs about 0.1 ms per frame of its own (its dispatches and
  the Hi-Z levels it asks for), the CPU's culling about 10 ns per instance in frames whose view changed (less where a
  batch is emitted whole). Measured (WebGPU, RTX, rotating, CPU submit, CPU / GPU culling): one unit of 100, 1000 and
  10k copies 0.17 / 0.33, 0.15 / 0.25, 0.11 / 0.23 ms; 100 units of 100 copies 0.24 / 0.29; of 1000 copies 1.12 /
  0.30 (`merging/instances-*`); the 60-copy assembly of 7qpd (`render/assembly-60-*`) stays on the
  CPU: culled on the GPU by 30 per draw alone, its frames drew 27 calls for 16 and submitted 0.07–0.12 ms more, the
  opaque pass unchanged. On the integrated GPU (AMD RDNA2, 2 CUs, `--force_low_power_gpu`), frames bound by the GPU:
  one unit of 100 to 10k copies CPU / GPU culled 5.6 / 5.7, 10.6 / 11.2, 26.2 / 26.2 ms; 100k instances 25.7 / 26.3
  rotating, 26.0 / 26.1 still (the culling pass 0.55 ms and the Hi-Z levels 0.18 ms on the GPU, the CPU's culling
  hidden behind it): nothing lost where the GPU bounds the frame. `postsynstubby_2` with pixel levels
  (`merging/postsyn-lod-pixel3*`, the frame, GPU culling and its occlusion / CPU culling / GPU culling without
  occlusion): whole 54.7 / 60.6 / 60.1 ms, zoomed 370 / 401 / 404, rotating 60.6 / 60.3 / 67.0 and 290 / 312 / 320;
  the culling pass 3.5–4.6 ms there, paid for by the occlusion. Re-measured with draw lists on (by default): the frame
  57.8 ms, zoomed 372, the opaque pass 55.1 and 365 of it (drawing the impostors); less data per instance copied (a
  3x4 transform for the 4x4, a quarter less) would save about 1 ms of 55 (2%) and change `aTransform` in every program
  of both backends: not taken up.
- **Not on the CPU**: which instances were drawn (the culled counts of timing stats), and the order of instances within
  a region, which varies between frames: it matters only for equal depths.

## GPU occlusion culling (plan, step 6.6)

Cells of GPU-culled merged objects also tested against a depth pyramid on the GPU, in `cell-cull` after the frustum:
for densely packed scenes (mesoscale), where most of what is in view is hidden behind the front.

- **Why**: `postsynstubby_2` (PetWorld, mesoscale preset, `merging/postsyn-*`) is bound by its opaque pass (4.7 ms
  whole, 20–23 ms zoomed to 0.3 of the distance), the same with CPU or GPU culling; only drawing less helps. The CPU
  Hi-Z (`hiZ.enabled`) cut the zoomed pass by 12% (23.3 to 20.6 ms) with the same draws: its readback lags (in
  back-to-back frames it is rarely ready) and the 200 Å cells are coarse, few whole cells are hidden. Not an estimate
  of what a GPU test can do.
- **The pyramid**: `HiZPass` builds it after the opaque pass from `depthTextureOpaque` (`r32f` levels, those from
  `minLevel` side by side in one atlas texture, read back for the CPU test). On the GPU the atlas is read directly
  (`textureLoad`): no readback, no lag beyond the frame.
- **One phase first**: the next frame's cells tested against the previous frame's pyramid with its view and projection
  (the CPU test's projection of a padded sphere, its level by the screen size, four texels compared, ported to WGSL).
  A cell that comes into view behind what moved shows a frame late (a frame of popping while the camera moves; the CPU
  test, with its readback lag, pops longer). Still frames are stable (the skip stays valid: a pyramid of the same view).
- **Two phases later, if the popping shows**: draw the cells visible last frame, build the pyramid from that depth,
  test the rest against it and draw what is newly visible (a second set of records per level and segment, the opaque
  pass split in two). Exact, at the cost of a pyramid build within the frame.
- **What occludes**: opaque depth only (transparent objects are occluded, never occlude). Off where a pass needs what
  is hidden: illumination (path tracing reads hidden geometry for shadows and reflections; the CPU test is off there
  too), and views without a pyramid of their own (they draw the CPU lists, as now).
- **Finer cells**: occlusion of whole cells needs cells small against the occluders; with the GPU testing cells
  cheaply, a smaller `cellSize` for GPU-culled objects (or a grid sized from the instance spheres) is part of it.
- **Draw lists**: `isChunkCurrent` refuses unplanned replays while a CPU occlusion test is set; GPU-occluded objects
  keep their recordings (their counts on the GPU), with step 5's view independence.

Steps:

1. **Measure the potential** (done): the GPU test in `cell-cull` against the previous frame's atlas
   (`OcclusionPyramid`, `HiZPass.gpuPyramid`, given by `Renderer.setOcclusionPyramid` and `prepareCull`; the pyramid
   rebuilt every frame, also while a read for the CPU test is pending; off with illumination as the CPU test). Found
   on the way, both of the existing Hi-Z: the levels went wrong from the second build after a prepare (the first
   level's `uInvSize` set in place into the vector of the level the cell held last: Canvas3D's prepares plan two
   levels, so the second; checked by `webgpu-hi-z`), and the test read the 2x2 texels around the rectangle's center,
   which miss a side where it is off a texel's center (now the texels it covers). With them, occlusion drew visible
   cells away (`postsynstubby_2`, zoomed: 4% of the pixels changed, the CPU test as well); fixed, no pixel differs
   (checked against the same frames without occlusion, zoomed 1, 0.3 and 0.1). Measured (`merging/postsyn-*`,
   WebGPU): the opaque pass 4.78 to 4.28 ms whole, zoomed 23.2 to 20.6 ms still and 20.4 to 18.1 rotating; cells of
   50 or 20 Å (default 200) 19.5 and 19.4 ms zoomed; the pyramid 0.02 ms, the culling pass unchanged. So for this
   model at most about 15%: most of the zoomed pass is the visible molecules, which culling whole cells doesn't
   reduce. They cost by the spheres drawn, not the pixels (`merging/postsyn-gpu-hiz-zoom-*`): the opaque pass 15.5,
   18.0 and 29.6 ms at 0.23, 0.92 and 3.7 Mpx (about 14.6 ms with no pixels, 4 ms per Mpx); by LOD preset at 0.92 Mpx
   5.3 ms `performance`, 18.0 `quality` (the default), 130 `ultra`. So the lever is the LOD levels and the cost per
   sphere impostor, not occlusion. The impostors (zoomed, 0.92 Mpx, 20.4 ms): `ignoreLight` 18.3, `approximate`
   23.2; measured by hand, the fragment shading made trivial 18.1, no fragments at all (the quads collapsed) 13.4, a
   triangle per sphere instead of two (3 vertices) 18.2. About 0.45 ns per sphere drawn, whole (11–12 M spheres,
   4.5 ms) and zoomed (44 M, levels 0 to 3: 4.1, 16.8, 20.3 and 2.9 M): about 48 spheres a pixel zoomed. The levels
   chosen by distance in Å, not by the size on screen, are the cost; the impostor's own work a tenth or two.
2. **Hi-Z for GPU culling without the readback** (done): `hiZ.gpu` (on by default) builds the levels whenever an
   object was culled on the GPU in the frame (`Renderer.gpuCulledCount`, `HiZPass.setGpuCulling`), also with
   `hiZ.enabled` off; the readback only for the CPU test. Measured (`merging/postsyn-lod-pixel*`, pixel levels,
   WebGPU): the opaque pass 1.23 to 1.09 ms whole, 6.25 to 5.70 ms zoomed (the same with Hi-Z off), the levels 0.03 ms.
3. **Grid size** (dropped): cells of 50 or 20 Å against 200 gained about 6% zoomed (step 1: 20.6 to 19.4 ms, distance
   levels), far less than the pixel levels; not pursued.
4. **Draw lists** (dropped with step 6.5).
5. **Two phases** (dropped): rotating 2° a frame (`-rotating`), 3 of 921,600 pixels differ from the frames without
   occlusion, zoomed; no popping shows.
6. **WebGL2, the CPU test** (measured, dropped): the benchmark's frames ran in one task, where WebGL's fences never
   signal (a sync's status changes between tasks only): the CPU test never had levels (every read timed out after
   `maxFrameLag`); on WebGPU only a read landing in the warm-up served the still frames (step 1's 12%). The Hi-Z cases
   now yield a task per frame (`yieldFrames`): reads land every 4 frames or so. Then (`merging/postsyn-lod-pixel3-zoom*`,
   WebGL2, RTX), against no occlusion (opaque pass 6.72 ms, CPU 3.83 ms, 11.6k draws): the cells within
   `2 * batchSize` of the camera tested (as now) 6.69 ms, CPU 5.62 ms, 12.5k draws; every cell tested 6.46 ms, CPU 17.1
   ms, 26.9k draws. A hidden cell splits the ranges `emit` joins, so the draws grow even with a test for free (e.g.
   visibility per cell computed on the GPU, read back): their CPU costs more than the opaque pass saves. WebGPU's test
   compacts what is drawn and its draws stay; the CPU test stays as it is (`hiZ.enabled` off by default).

## Steps

1. **Benchmarks** (done): `merging/units-N-<spacefill|cartoon>-<merged|unmerged>[-hover|-recolor]` (on demand; N 10
   to 1000 chains of 10 atoms, one representation; the update before each frame timed as a stage of its own, after a
   first update outside the timing so that its programs are compiled), with setup times and GPU bytes. Per frame,
   merged / unmerged, WebGPU (NVIDIA, D3D12): spacefill still 0.27 / 0.27 ms (N 10), 0.57 / 0.58 (100),
   0.70 / 1.59 (300), 1.68 / 2.53 (1000, submit 0.48 / 2.40); hover 1.84 / 4.54 (1000); recolor 13.1 / 29.2 (1000,
   almost all of it the representation's update on the CPU); setup 162 / 257 ms (1000), 19 / 7 (10); cartoon still
   3.41 / 3.79 (1000), hover 1.63 / 4.42, recolor 13.2 / 28.4. WebGL 2 at 1000: spacefill 1.87 / 7.99, hover
   4.72 / 16.4, recolor 16.4 / 36.2; cartoon 1.88 / 5.13. GPU bytes about equal (members have no GPU items of their
   own). Findings: merging pays from a few hundred units and is neutral below (within the noise), never clearly worse
   from 30 units on: `auto` from more than 30 units. A merged object issues as many draws as the units (a draw per
   segment: each unit's geometry is its own, a draw spans one segment's vertices for its instances): what merging saves
   is the CPU (one object walked, nothing rebound between the segments' draws), not draws. Updates dominate at scale:
   a color change of 1000 units 13 ms merged, mostly per-unit visual and theme work. Not measured yet: a geometry
   update (a trajectory frame), CPU memory.
2. **Profile the updates** (done; hover and recolor at 1000 units). The recolor's CPU went to checks repeated per
   member and update: `hasForeignResources` (every texture and buffer of every render object, each update; now once
   per object and context, `render-object.ts`), `getMergeKey` (the key string built per member and update; now cached
   per values, rebuilt when a key cell's value changed: updates bump versions of unchanged values), and the merged
   `sync` taking a snapshot of all members' versions after each sync (now updated where compared). Recolor at 1000
   units, merged / unmerged: WebGPU spacefill 5.11 / 22.2 ms (13.1 / 29.2 before), cartoon 5.32 / 21.1; WebGL 2
   spacefill 6.69 / 32.0 (16.4 / 36.2), cartoon 6.85 / 23.8. Hover merged 1.25 ms (WebGPU), 5.75 (WebGL 2). What
   remains of a merged recolor (about 4.7 ms) is the units' own theme and visual work (location iteration, the color
   textures) and the copies into the merged layers: the per-member updates of step 5. Later, the render items'
   `update` compared every value's version by name in one object of all the values (a dictionary): now arrays per
   list (defines, attributes, textures), and WebGPU's theme layers return early where none changed; the update stage
   of unmerged recolor at 1000 units 6.4-6.8 to 5.0-5.2 ms (WebGL 2), 4.9-5.2 to 3.7-3.8 (WebGPU).
3. **Geometry without copies** (done): the merged vertex attributes and elements are `SegmentedArray`s of the members'
   arrays (`segmented-array.ts`), written by part into the render item's buffers (WebGPU geometry pages, WebGL
   buffers), an update writing the parts whose member changed; indices biased by the member's vertex base in a
   temporary while written. The geometry store ([shared-draw-resources.md](shared-draw-resources.md) step 2) is not
   needed for it and deferred (measured there). The geo exporter exports a merged object's members (their own values;
   the merged values are segmented). At 1000 cartoon units (WebGPU): the merged CPU geometry copy 51.8 MB to none,
   frames 1.45 ms merged / 2.23 unmerged, add and commit 26 / 169 ms, GPU bytes unchanged; recolor 4.89 ms.
4. **Layers without copies** (done): the merged theme layers are images of the members' layers as parts, written by
   whole words (WebGPU data pages) or rows (WebGL, a few MB at a time) composed from the parts (`SegmentedArray.read`:
   a word or row may hold texels of two members), an update writing the changed parts' words or rows. Layers grow
   with groups times instances (per group and instance colors, overpaint, transparency): measured with
   `merging/assembly-60-*` (60 copies of 7qpd, its chains the units, `element-index` colors), the merged copy of the
   spacefill layers 3.5 MB to none (a 1M-atom unit with 60 copies: 240 MB of colors alone); frames and recolor
   unchanged within the noise (the recolor, 68 ms, is the color theme). Checked on Dawn: a merged object's markers
   written by part draw what a whole upload draws (the test fails with the part writes left out); WebGL's partial row
   writes are not checked on a GPU yet (no native GL on the development machine).
5. **The costs per member left** (in part): measured at 1000 units. Location iterators, created per visual and update
   (a theme change, one per unit), as a class instead of an object of closures and accessors
   (`location-iterator.ts`): recolor merged 4.57 to 3.60 ms (spacefill), 4.88 to 3.73 (cartoon), unmerged 19.2 to
   18.5; the color loops' `move` faster as well: the 60-copy assembly recolor 70.7 to 57.7 ms. The benchmark's warm-up
   now draws frames until a few in a row draw without skipping (an update's programs are created by the frames drawing
   them and compile in the background; the first run measured frames without the updated objects, and its draw
   count, the one printed, showed 9 calls). Then, profiled (`merging/units-1000-spacefill-merged-recolor`, WebGPU,
   3.32 ms): the merged `sync` 12%, its layout key a string over all members built every update, and about 30 named
   reads of each member's values (objects of a hundred keys: slow named reads); now the cells resolved once per member
   into objects of one shape, the versions seen in arrays, the layout compared as numbers (`sync` 12% to 5%). And the
   clip objects allocated and their cells bumped per visual and update (`Clip.getClip`): now only where the props' clip
   changed. The recolor 3.32 to 2.49 ms (cartoon 3.56 to 2.84; WebGL2 2.32 and 2.48). Left: the visuals' own updates
   (values, `prepareUpdate`'s comparisons, about 40%), the merge checks (`getMergeKey`, 9%: 26 version reads a member,
   scattered cells, no change seen; only notifications would avoid them), garbage (10%); no single large item. Then
   the props kept as one object while none changes (units and complex representations, `shallowMerge2`; the update
   made a new one every time): the visuals skip the props' comparisons and the values from the props (`updateValues`,
   `updateRenderableState`, some 20 cells a visual) where neither the props nor the geometry, transform or sizes
   changed. Checked by `props-identity` (fails with either left out). Recolor at 1000 units (WebGPU / WebGL2): merged
   spacefill 2.49 / 2.32 to 1.75 / 1.50 ms, cartoon 2.84 to 1.98; unmerged spacefill 8.58 to 6.76, its frame 10.9 to
   6.4 (fewer cells bumped, fewer uploads). Complex visuals merged their props once more (a new object every update):
   now against their own last props (`shallowMerge2`; callers may change theirs in place, as the particle
   representations do). Not extended to particle, volume and shape representations: a visual per representation, the
   props' work a few µs per update against the data's (a mesoscale recolor: per type, millions of particles).
   Allocations of the recolor (`pnpm benchmark --heap-profile`, sampled with the objects collected since): 453 MB over
   the run under the visuals' updates, a fifth of it `Object.keys` of the themes' `deepEqual`, compared per visual for
   the same pair of themes (now cached per pair, `ColorTheme.areEqual`, `SizeTheme.areEqual`), a third the update
   function's promise per visual (`createOrUpdateVisual`, now a promise only where a visual's update is one, one
   context object per update) and a quarter a new location iterator per visual and recolor (now kept while the group,
   structure and props stay; the color and size writers reset it): 166 MB. Recolor at 1000 units (WebGPU): merged
   spacefill 1.75 to 1.41 ms, cartoon 1.98 to 1.55. Left: the visuals' update closures (`withPreparedGeometry`), small.
   Mesoscale (`merging/postsyn-lod-pixel3-recolor[-entity]`, particle representations, 517k particles): a uniform
   color per type 0.85 ms an update; the entity colors (by name and by function in turn, a color per particle) 15.8
   ms, the palette interpolated per particle (`palette.color`): now a color per entity or compartment computed once
   (`particle-entity`, `particle-compartment`), 12.1 ms (WebGL2 14.4 to 11.0). Then the instance colors written by
   index, each instance's location without moving the iterator (`LocationIterator.getInstanceLocation`; `move` was a
   third): 7.0 ms (WebGL2 6.1). And the particles' location by their index set's kind resolved once (an interval or
   an array, not `OrderedSet.getAt` per instance): 5.9 ms (WebGL2 4.8). Large assemblies
   (`merging/assembly-60-*-merged-recolor`, 60 copies of 7qpd, `element-index` colors per group and instance): 58 ms,
   a third a binary search per location for the element's index in its unit; now by subtraction where the unit's
   elements are consecutive, the unit's index cached for consecutive locations: 25 ms (cartoon 8.4 to 4.0); the rest
   the palette's interpolation per element and copy. Measured and dropped: the group and vertex color writers by index
   as the instance one (no gain there; the instance writer itself slower, 7.5 to 9.5 ms, by the JIT).
   On WebGL the segments' draws without base instances re-pointed the instanced attributes per draw (`changeOffset`,
   about 40% of a hover frame); now (WebGL2 without the base instance extensions) merged render items read their
   instance data (transforms, instance ids, segment data: 7 texels an instance) from a texture at the draw's first
   instance (`uInstanceBase`, set per draw as `uDrawId`) plus `gl_InstanceID` (`dInstanceData`, macros in
   `common-vert-params` in place of the attributes, which stay declared and unread; programs offset only the
   attributes they read), filled by the WebGL render item from the instance values (`getInstanceDataSchema`). At 1000
   units (WebGL 2): spacefill still 1.59 to 0.47 ms, hover 4.39 to 1.49; cartoon still 1.59 to 0.71, hover 4.40 to
   1.89, recolor 4.82 to 3.85. WebGL1 keeps the offsets. Profiled zoomed (`merging/postsyn-lod-pixel3-zoom-nohiz`,
   WebGL2 on D3D11, no base instances: 11.6k draws, CPU 3.7 ms): culling 61% (`cullSegment`, `emit`, the frustum
   tests), drawing 35% (a draw per range, `offsetAttributes` 8%, the uniforms 7%). The offsets skipped for programs that
   read no instanced attribute (`Program.readsAnyAttribute`, the instance data from the texture): CPU 3.72 to 3.47 ms.
   Measured and dropped: multi-draws (`WEBGL_multi_draw`, a list a call) with each draw's first instance read by
   `gl_DrawID`: CPU 3.47 to 2.6 ms, but the opaque pass 6.69 to 7.7 ms zoomed and 1.32 to 1.46 whole from a texture
   (read once a vertex: one more dependent read to the seven of the instance data), 8.5 and 1.59 from a uniform array
   (ANGLE's dynamic indexing on D3D); `gl_DrawID` itself costs nothing. The frames there are bound by the GPU.
   Named reads per frame (values objects of a hundred keys), audited: culling read about ten cells by name per object
   (and per member of a merged object); now resolved once per renderable and member (`CullValues`). WebGL draws without
   base instances re-pointed the instanced attributes per draw and back to 0 after, also where every draw's first
   instance is 0 (an instance per unit); now only where it differs. The draw paths (the render items, the uniform
   writers) were array-based already. Measured (submit, WebGPU / WebGL2): `merging/units-1000-spacefill-merged` 0.60 /
   0.34 to 0.31 / 0.31 ms, `-unmerged` 2.35 / 7.48 to 2.31 / 5.64 (the hover frame on WebGL2 15.4 to 13.1),
   `merging/instances-1000-spacefill-merged-rotating` 0.92 / 0.90 to 0.89 / 0.86. Left on WebGL2's unmerged frames:
   texture binds per object (each texture bound, no cache of the units' bindings: every direct bind would have to
   invalidate it).
6. **GPU culling and levels** for merged objects ([draw-table.md](draw-table.md) step 7), WebGPU; CPU culling stays the
   fallback. Fewer draws than segments would need a segment per vertex (the draws span segments), not instances.
   Measured first (`merging/*-rotating`, `merging/instances-<copies>-spacefill-*-rotating`: 100 units of 100 or 1000
   copies, culled per frame by their instance grids): at 1000 units of one instance the CPU submits a rotating frame in
   0.32 ms (culling a small part); at 100k instances (100 units of 1000 copies) culling was about 90% of the frame,
   2.29 ms merged (2.70 unmerged), through `Frustum3D`, `Plane3D` and `Sphere3D` per batch and cell. With the planes as
   coefficients and the spheres read from their arrays (`cullSegment`, `renderable.ts`): 1.24 ms merged (WebGL 2 2.71
   to 1.15), 1.52 unmerged; what remains is about 80% of the frame, the frustum tests and the draw ranges' merging
   (`emit`), about 1 ms per 100k instances per frame, linear in the cells tested. GPU culling would take that off the
   CPU: worth it from some 100k instances in view-changing frames (large assemblies), not below. Designed above (GPU
   culling), in parts:
   1. Measured (done): `draw-overhead/10k-shared-first-instance-indirect*-bundle` (bindings set once, a draw per object
      by its first instance, NVIDIA, D3D12): 10k indexed indirect draws replay in 2.76 ms per frame, 0.28 µs each,
      against 0.99 ms direct; draws with an instance count of 0 cost the same (half 2.74 ms, all 2.80 ms). GPU culling
      replays segments times levels draws in every frame (still ones too), about 0.28 µs each; CPU culling costs about
      10 ns per instance in frames whose view changed (`merging/instances-*`), nothing in still ones. Break-even about
      30 instances per segment and level (100 units of 1000 copies, one level: 28 µs of draws against 1 ms of
      culling; 1000 units of one instance stay on the CPU): the threshold of step 7, `indirect-first-instance` now
      requested where available (`OptionalFeatures`).
   2. Indirect draw lists (done): `MultiDrawBaseData.indirect` (a storage buffer and byte offset), a record per entry
      of the list (`IndirectRecordBytes`: indexed 20 bytes, else 16) in place of its arrays; the WebGPU render item
      draws them by `drawIndexedIndirect` or `drawIndirect`, WebGL throws (none made for it); the draw key of an
      indirect list is its buffer and offset (the GPU writes the arguments: recordings stay valid). The capability
      `indirectInstanceRanges` (WebGPU with `indirect-first-instance`, false on WebGL). Checked on Dawn
      (`webgpu-indirect-draws`): merged spheres and meshes, colored by chain, drawn from indirect lists of their
      segments (first instances other than 0) draw what the same direct lists and their own draws draw; the test fails
      with the first instances left out or the byte offset ignored.
   3. GPU culling of merged objects (done): `gpuCulling` (off by default at first; a merged object of the other kind
      is a new one). The culler (`merged-culler.ts`) writes the cells, each instance's cell and segment, the instance data
      (`packInstanceData`), the records and the identity region when the layout or the instances change; per frame
      `cell-cull`, `instance-copy` and `cull-args` (`gl/compute/cell-cull.ts`) after zeroing the counters. Storage
      buffers can be vertex buffers (the render item binds a storage buffer of an attribute as it is). Checked on Dawn
      (`webgpu-gpu-culling`): merged spheres with cells of a few instances, GPU-culled, draw the instances the CPU
      culling draws (the records' counts read back) and the same image; the test fails with a transform column not
      copied; `Renderer.update` never culls (a compute pass would end the pass it is called in); a canvas with draw
      lists draws and picks. A dispatch is skipped while the planes, the distance range, the layout and the instances are those of the last (still frames cull nothing again). The instances of a workgroup (consecutive, mostly of one segment) reserve their slots by one global atomic per level (a counter per segment for 100k instances contended: the copy 0.08 to 0.02 ms). Measured (`merging/instances-1000-spacefill-*`, 100k instances, WebGPU): rotating, the frame 1.14 to 0.79 ms, CPU submit 1.05 to 0.15 ms, the culling passes 0.02 ms on the GPU, the opaque pass 0.55 ms against 0.59; still, 1.32 to 0.76 ms (the CPU culls again every frame, the GPU not). Measured and dropped first: a test per instance, read from storage through compacted indices (the opaque
      pass 0.65 against 0.58 ms, the per-vertex reads); cells as draws (51k records, 21 ms per frame). The kernels
      (`MergedCullKernels`, the same for every object) compile once a culler is created, and with the programs where
      a precompile list names a GPU-culled object (`ObjectSignature.gpuCulled`: recorded lists, and lists derived from
      params where `gpuCulling` is on; 'auto' depends on the instances, derived as culled on the CPU); checked on Dawn (`webgpu-gpu-culling`: a context with nothing in its scene
      has them after precompiling a recorded list; fails without). Hi-Z: step 6.
   4. Levels per instance (dropped): a cell across a level boundary draws its instances at both levels' ranges; about
      4% of the vertices of the cell path are such duplicates (`postsyn`, pixel levels), not worth a test per instance.
   5. Draw lists, GPU-culled objects not view dependent (dropped): draw lists are off by default, and a GPU-culled
      object's submit is 0.15–0.25 ms without them.
   6. Occlusion on the GPU (done): per cell, above (GPU occlusion culling, steps 1 and 2).
   7. The threshold (done): `gpuCulling: 'auto'` by default (above, When on); an object alone is merged to be culled
      on the GPU (`canMergeRenderObjects` with one member), WebGL keeps it as it is. Checked by `merged.test` (which
      objects merge, and how) and on Dawn (`webgpu-gpu-culling`: one object merged alone, and not with 40 instances
      by default).
   8. CPU culling by batch (done; WebGL, and WebGPU without `gpuCulling`): a batch inside the frustum, within the
      distance range and one level alone, and too far for the cells' occlusion tests, emits its instances as one
      range (`cullSegment`; the cells' instances are consecutive, the draws the same as cell by cell, tested by
      `cull-segment`). Measured (`merging/postsyn-lod-pixel3*`, WebGL2, RTX): the culling 3.6 to 0.2-0.3 ms whole (162k
      cells, 6k batches), 3.6 to 2.3 ms zoomed (batches across level boundaries and the view's edge still by cell);
      the frame whole 4.09 to 1.44 ms, zoomed GPU-bound (7.1 to 7.0 ms, its CPU 6.6 to 4.0 ms);
      `merging/instances-1000-spacefill-merged` 1.48 to 0.95 ms. Measured and dropped: joining a level's draws across
      gaps of up to 16 instances (to the same image: the gap's cells culled, its spheres dropped by distance or
      clipped): 11.7k to 2.8k draws zoomed, CPU 5.6 to 3.4 ms, but the opaque pass 6.6 to 7.8 ms (the gaps are mostly
      instances at a coarser level, drawn with the finer level's spheres).

## Open questions

- Whether units of different representations (same kind and key, different materials) should merge as well (a
  material per segment): only if measurements show many representations of few units each.
- Picking and the object ids of members after updates per member.
