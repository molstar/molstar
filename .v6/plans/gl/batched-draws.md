# Plan: Batched draws, many objects per draw

## Status

- **Dropped** (the code removed; this plan kept as the record). Batching served scenes of many small objects in
  separate representations, which Mol* doesn't have beyond about 1k objects; many objects come as many units of few
  representations, which CPU-side merging serves better: `units-1000-spacefill` merged 0.43 ms per frame, unmerged and
  batched 2.43 ms, unbatched 2.38 ms (one object walked and drawn against a thousand). At 1k separate objects the
  unbatched path with shared object groups and draw lists is enough. Against it: two programs per set of defines (solo
  and merged), uniform array windows, member bookkeeping, material arrays, about 1.7k lines in the WebGPU backend, and
  large sphere objects about 7% slower batched (cause not found). Merging is planned instead:
  [merged-render-objects.md](merged-render-objects.md).
- Kept from the work: shared uniform pages and object groups ([shared-draw-resources.md](shared-draw-resources.md)
  step 4), sphere and point positions as a data layer, clip objects as a frame and a shape (smaller material blocks,
  cheaper clip tests), the draw table's CPU work (chunks replayed after view changes, visibility from the table), the
  `draw-overhead/*` measurements (indirect and merged draws).
- Lesson: experiments that read the object block at a constant index read another object's block (the index wasn't 0)
  and drew other data; their speed-ups were not valid. An experiment must draw the same image as what it is compared
  to (check it, or keep the data correct by construction).

## Why

Scenes of 10k–100k small objects (user data: shapes, annotations, many small structures) draw one object per draw. With
shared geometry pages, theme layers in storage and shared object groups ([shared-draw-resources.md](shared-draw-resources.md))
a draw still binds its object group at its dynamic offsets, its vertex buffers and index buffer at its ranges. Measured
with `draw-overhead/*` (raw WebGPU, 10k tiny cubes in one bundle, NVIDIA, D3D12), per frame:

| How the draws are issued                                               | ms   |
| ---------------------------------------------------------------------- | ---- |
| 10k draws, a shared bind group set per draw at a dynamic offset        | 1.1  |
| the same, as indexed indirect draws                                    | 4.1  |
| the same, half of them with an instance count of 0 (culled)            | 4.1  |
| 10k draws, nothing set between them (objects by `firstInstance`)       | 1.0  |
| one draw of all, instanced                                             | 0.16 |
| one draw of all, merged (rebased indices, an object id per vertex)     | 0.15 |
| merged, half the objects culled at random: 2469 runs, nothing rebound  | 0.29 |

Mol\* itself: `draw-table/boxes-10k-drawlists` 4.75 ms. Rebinding between draws costs, issuing draws barely does; an
indirect draw costs its replay whether it draws or not. So culling (CPU or GPU) gains GPU time only where it means fewer
or cheaper draws: runs of consecutive objects in a merged draw.

## Design

### Batches

- A **batch** holds objects that draw with one pipeline and one set of bindings: the same program (shader, defines),
  material layout (material blocks read by index), object group (pages) and geometry page. The usual scene is a few
  representations, each of many render objects, each of many instances: objects of any instance count batch, their
  instances in a storage array of the batch read by `instance_index`.
- Its geometry is **merged**: the objects' vertex ranges and index ranges allocated in the batch's region of a
  geometry page, indices rebased to the region (`index + vertex base`), and an object index per vertex (u32,
  `aObject`). A draw of a range of the batch's indices draws the objects in it, whatever its length.
- **Object data by index**: the object blocks of the batch in a storage array with the block's stride; batched program
  variants read `uObject` as `objects[aObject]` (vertex) and pass the index flat to the fragment stage. The theme
  layers' offsets are object uniforms already, so layers need nothing else. Material blocks likewise, from a storage
  array per material layout, by an index in the object block.
- **Instances by index**: the batch's instances (transform, id) in a storage array, an object's start in its block; the
  shader reads instance `start + instance_index` in place of the instanced attributes. A member draws its index range
  with its instance count (and the culled or level ranges with their base instances); consecutive single-instance
  members draw as one range (step 4).
- **Order**: a batch draws in its index order, not the table's: for passes where order does not matter (opaque,
  depth, picking, marking, the transparent passes of order-independent transparency); blended passes in sorted order
  keep per-object draws.

### Drawing

- The renderer walks a pass as before (the table's order, the row checks); the members of a batch drawn one after
  another bind nothing between them (shared groups, buffers, draw block), and the pass state merges their draws while
  each continues the index range of the one before (`PassState.drawIndexedRun`): runs of consecutive objects, one
  `drawIndexed(count, 1, firstIndex)` per run. All visible: one draw; a hover marking one object: one draw in the
  marking pass. The order of the draws is kept, so every pass merges (also blended ones); runs form where the table's
  order follows the batch's index order (objects added together).
- Draw lists record the merged runs per chunk; a change of what a pass draws records the chunks it touches.
- **GPU culling** then per batch: a kernel tests the batch's rows (`DrawTableBuffer`, since removed), compacts the indices of the
  visible objects into a per-pass index range and writes one indexed indirect draw: the indirect cost per batch, not
  per object; recordings stay valid while the camera moves.

### Memory and updates

- Per vertex 4 bytes more (`aObject`); the rebased indices replace the object's own (an object in a batch keeps no
  separate index range).
- An object's geometry changing size gets a new range in the batch (the old freed): runs may split, the batch stays
  drawable. Removing objects leaves holes, compacted when the batch's fill is low.
- The batch regions come from the geometry store ([shared-draw-resources.md](shared-draw-resources.md) step 2): the
  builders write into the batch's region where the batch is known at build time, otherwise one copy on adding.

### WebGL

- The same merged ranges draw in WebGL 2 (`drawElements` with an offset per run, object data from a data texture by
  index); WebGL keeps per-object draws until measured.

## Steps

1. **Upper bound** (done): `draw-overhead/10k-merged-bundle` and `-merged-culled-bundle`, see Why.
2. **Object blocks by index** (done, meshes): with `dBatched` (`batched-objects.ts`) the object block is declared as
   `sObjects: array<Object>` (read-only storage, bound at the batch's first block with a dynamic offset) and a private
   index `gObject`, every `uObject` of the program reads `sObjects[gObject]`; the mesh shader sets the index first in
   both stages (`assign_object_vertex_wgsl`: from `aObject`, passed on flat as `vObject` at location 13). The buffers
   of uniform blocks (ring, arenas, uniform pages) have storage usage and an alignment valid for both bindings. Checked
   with `WebGPUContextOptions.batchObjects` (each mesh a batch of its own, `aObject` zeros): the `webgpu-batched`
   visual profile (all cases pass) and `webgpu-batched-objects.test.ts` on Dawn (images equal to the uniform blocks';
   fails when the batched stage discards or binds another block). No cost measured on `boxes-10k-drawlists` and
   `units-1356-gaussian-drawlists-rotating`. Out-of-range indices read the binding's last block (robust access).
3. **Merged geometry** (done, meshes; with `batchObjects`): `object-batch.ts`, a batch per key (program name, schema,
   defines, material layout, transparency, filtered textures, data page of the theme layers) holding its members'
   object blocks as one storage array (a range of the uniform pages, bound whole at its start; `aObject` the member's
   index), their instances as another (a range of the data pages, five `vec4` per instance, `uBatchInstanceStart`) and
   their geometry merged: a range of the geometry pages per vertex buffer slot and one of indices rebased to the merged
   vertices. Items join while they can draw from one (indexed, no segments, indirect draws or object or material
   textures used by their defines), leave when not, and draw their own index range of it with their instances (or
   their culled and level ranges at their base instances); their own geometry is released while in a batch. Growth
   replaces the batch's ranges and writes all members again from their CPU data (queue writes only: a GPU copy would
   run after the writes of its submission). Object blocks of members are written in place when changed, not
   snapshotted per submission. Checked: `webgpu-batched` visual profile (all pass, also with draw lists;
   `geometry/mesh/instanced` fails with the instances not indexed), Dawn (`webgpu-batched-objects.test.ts`: a cartoon
   of two chains in one batch, two cartoons of other materials in another, fewer index ranges bound, images equal;
   fails with the indices not rebased or the material index wrong), `object-batch.test.ts`. Measured (each member
   still a draw of its own, nothing rebound between them but the draw block's offset): `draw-table/boxes-10k-drawlists`
   5.08 to 2.84 ms per frame, `-drawlists-sweep` 9.4 to 5.7 ms, `boxes-10k-sweep` (direct) 120 to 89 ms,
   `units-1356-gaussian-drawlists-rotating` about equal (`--backends webgpu,webgpu-batched`).
   3b. **Material blocks by index** (done): the material blocks of a layout in a storage array (`getMaterialArray`,
   a range of the uniform pages, an index per material on first use, written when changed), `uBatchMaterial` in the
   object block, `gMaterial` set with the object; the material id out of the batch key (the material layout in it), so
   that objects of different representations share batches (`boxes-10k`).
4. **Merged draws** (done): `PassState.drawIndexedRun` merges the one-instance draws of consecutive members continuing
   each other's index range, issued with the next change of the pass state or `flush` (the context flushes before
   bundles, clears, the blend constant, the pass's end, a recording's end; render items before drawing directly).
   Nothing is bound between members: the draw block holds what the pass sets only (`uLod` in the object block,
   `hasObjectLod`), shared by the members of a batch (`getBatchDrawBlocks`, an arena of the batch). Members of several
   instances draw one after another, nothing bound between them either. Checked as step 3 (all visual profiles pass,
   also with draw lists). Measured (`--backends webgpu,webgpu-batched`, per frame): `draw-table/boxes-10k-drawlists`
   5.17 to 1.66 ms, `-drawlists-rotating` 7.0 to 4.0 ms (CPU 3.8 ms of it), `-drawlists-sweep` 11.9 to 4.4 ms (CPU
   4.2 ms), `boxes-10k` (direct) 64 to 38 ms, `boxes-10k-sweep` 143 to 79 ms, `units-1356-gaussian-drawlists-rotating`
   about equal. The GPU waits 0.1–0.2 ms where the passes replay; what remains is the CPU of the passes' walks
   (`draw-table.md`).
   Then the CPU side (profiled with `--profile`): chunks replay without walking after view changes where their rows'
   visibility is unchanged (`draw-table.md`, Chunks), the scene's visibility from the table's `visibleVersion` (not a
   hash over all objects per hover), the last pipeline of a program memoized (`WebGPUProgram.lastPipeline`), pending
   ranges of `BufferPages` released once per submission (setup was quadratic in the frees), batch members not written
   again on their first draw. Per frame with `webgpu-batched`: `boxes-10k-drawlists-rotating` 4.0 to 1.5 ms (CPU 3.8
   to 0.42 ms), `-drawlists-sweep` 4.4 to 2.8 ms, `boxes-10k-sweep` (direct) 79 to 66 ms. What remains of the sweep:
   the marking passes record the chunks of the hovered objects (two per frame), each draw of a recording captured
   again for its digest (`captureDraw`, about a sixth), the marked objects' updates.
   4b. **The default** (done): `batchObjects` on unless false (`webgpu-unbatched` visual and benchmark profiles keep the
   unbatched path checked). One program per set of defines, whatever the object's size, so that precompiled programs
   are the ones drawn: program items predict the batch from the defines (`isBatchedByDefines`). A size limit (large
   objects unbatched) was tried and dropped: a second program per defines, unknowable to precompiling. Instead the
   blocks moved from storage to uniform arrays bound a window at a time (`getBlockWindow`: 144 object blocks of meshes,
   24 material blocks, both windows bound at their offset, the index modulo the window): blocks read from storage made
   `render/globule-300k-spacefill-mesh` 11.1 to 13.1 ms per frame (material 1.1 ms of it, object 1.3 ms; copying the
   blocks into private variables once per invocation 14.6 ms), uniform windows 11.9 against 11.4 ms unbatched; the
   instances stay in storage (as fast as instanced attributes). The draw block is bound only where the code reads it.
   Runs break at windows: `boxes-10k-drawlists` (a material per box) 1.2–1.4 ms with storage arrays, 1.85 ms with
   windows (unbatched 5.4 ms), `-drawlists-sweep` 3.4 ms (9.5 ms); the material block was large (2400 bytes, the clip objects' arrays of 16 per material). The clip objects as a frame
   and a shape each (`Clip.Objects.frame`, `shape`; 80 instead of 144 bytes per object, one matrix per test instead
   of a matrix and a quaternion) make it 1376 bytes, 40 per window: `boxes-10k-drawlists` 1.19 ms again. Clip sets out
   of the material block (an index per material, the sets in a uniform array of their own: 680 materials per window)
   remain possible, for scenes of many distinct materials.
5. **The other geometry kinds**: cylinders and lines (done; indexed, no object textures: the mesh's hooks in their
   shaders, `BatchedPrograms`; checked by their visual cases blank with the batched fragment stage discarding). Batches
   count their ranges as resources (`stats.resourceCounts`). Open:
   - spheres and points, in three steps: (a) `tPositionGroup` (positions and groups, read per item from
     `vertex_index`) a data layer in the data pages, `rgba32f`, so that they have no object texture; (b) unindexed
     members: a member is a range of the batch's vertex indices (six per item), `aObject` per vertex as for meshes, its
     first vertex in its object block, the item read at `(vertex_index - first) / 6`; runs of consecutive members
     drawn unindexed (`draw(count, 1, first)`); (c) levels of detail, below. (a) and (b) done: spheres and points in
     `BatchedPrograms`, unindexed members with `PassState.drawRun`, `vertex_index` relative to the member's first vertex
     (`uBatchVertexStart`, also for indexed members: batched meshes, cylinders and lines had read their per-vertex data
     at the rebased index); spheres keep `uLod` in the draw block (`DrawLodPrograms`) until (c), so that their runs
     don't merge yet. Checked: all sphere and point visual cases draw batched (blank with the batched stage
     discarding). Measured: `units-1000-spacefill-unmerged` (1000 sphere objects) 3.15 to 2.78 ms, but large sphere
     objects slower: `globule-300k-spacefill` 4.40 to 4.70 ms, `-large` 14.4 to 15.6 ms (step 5d);
   - 5d. **Solo and merged programs** (done): two programs per set of defines (`dBatchSolo`, a define, so solo and
     merged members are in batches of their own): solo, the object index from the draw block (`uBatchObject`, uniform
     over the draw), for members drawing alone (several instances, levels of detail, more than `MaxMergedVertices`, 512);
     merged, `aObject`, for runs of small one-instance members (their draw blocks shared without `uBatchObject`).
     Precompiling prepares both (program items of batched objects: the size is not in the signature). Measured
     (`globule-300k-*`, one large object each): meshes 12.0 to 11.4 ms solo (unbatched 11.5 ms); spheres still about 7%
     slower than unbatched solo (4.9 against 4.6 ms, `-large` 15.9 against 14.7 ms), cause open: not the per-vertex
     index. Experiments reading the object block at a constant index read another member's block (the globule's index
     isn't 0) and drew other data: their faster times are not valid;
   - levels of detail (spheres): an object draws a level as a prefix of its items (the items ordered so that a prefix
     subsamples: `lodLevels`, the level's count, distances, overlap and radius scale), the levels' prefixes overlapping.
     Their parameters are per object and change with the props only: a level table in the batch's storage (per member
     its levels, `uLod` of each), read by level. The level of a draw can't come from `vertex_index` (the prefixes
     overlap): it comes from `instance_index`, the draw's `firstInstance` holding level and first instance
     (`level * levelStride + baseInstance`); the shader reads level and instance from it. Which levels and instance
     ranges a frame draws is decided on the CPU today (`cull`, `MultiDrawBaseData`), direct draws, any
     `firstInstance`; GPU culling (step 6) writes the same ranges as indirect draws, where a nonzero `firstInstance`
     needs `indirect-first-instance` (optional feature; without it, a level's draws rebind a per-level offset). The
     level table is data the culling kernel reads as well (distances per level);
   - text and image: object textures (the font atlas, the image), not batched.
6. **GPU culling** per batch ([draw-table.md](draw-table.md) step 7; done instead per merged object,
   [merged-render-objects.md](merged-render-objects.md) step 6): a kernel tests the members' rows and instances,
   picks their levels of detail from the level tables and writes indirect draws (level and first instance in
   `firstInstance`, see step 5).

## Follow-ups

- **Merged render objects** (`renderMerged`, CPU-side: many units or visuals in one render object, their members
  segments, `dSegmented`) are one object: their block is read uniformly in every draw, also in GPU culling's indirect
  draws and per level (the level per draw); their per-instance data (`aSegment`, `aSegmentSphere`, `aSegmentLod`) read
  by instance from storage like the batches' instances (their cost to be measured: see 5d). Not batched today (`isBatchedNow`): batching them means
  the segment attributes in the batch's instance records. Analyse and plan making merging cheaper and more general
  (building, updating a member without rebuilding the whole, more geometry kinds) as an alternative to batching for
  scenes of many small objects of few representations; with the geometry store
  ([shared-draw-resources.md](shared-draw-resources.md) step 2: members contiguous in a page drawn as ranges, no copy
  into merged values).

## Open questions

- Batch keys where facts differ by pass (a marked object drawn in the marking pass only): runs handle membership, the
  pipeline's render state must not differ within a batch.
- Picking ids and per-object lod in a batch (the object block holds them; lod as a cull before the runs).
- How large a batch: one per program and material and page, or bounded (many objects changing at once rewrite runs of
  one batch).
- Clip objects and solid interior caps of batched objects.
