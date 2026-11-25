# Plan: many objects culled on the GPU (the mesoscale explorer's split)

## Status

- **Current**: revised (below, "Revised: the pool removed, fewer objects instead"): steps 5 and 6 (the culling and
  the raster over all objects, the culling pool) removed (step 8), steps 3 and 4 kept; option 1 as a budget done
  (step 13: the split back to about the `noraster` level); next steps 9-12. Their measurements below stay as the
  reference of what the split costs and what a pool gains. Benchmark cases
  `mesoscale-split/postsyn-{whole,zoom0.3}[-rotating]-{noraster,off,auto,auto-objects{0,4,16,all}}` (on demand,
  `postsynstubby_2-no_bonds.bcif` in `.benchmark/`): PetWorld as the mesoscale explorer builds it, a spacefill
  representation per model (150 objects, `StructureFromPetworld`: a model's units merged, instanced by its operator
  group), its `quality` parameters, a uniform color each, Hi-Z on. `noraster` as on a device without `atomicMinMax64`,
  `auto` as on one with it, `off` the objects of `auto` with the raster off.
  - **Step 5, the culling over all objects**: the cullers' buffers are ranges of pages (`merged-cull-pool.ts`, a page
    per record stride and binding limit), an object a row of parameters (`CullObjectVec4s`), the three kernels over a
    page's cells, instances and records by global index (`cellObjects`, `recordObjects`, the instance cells as
    `vec4<u32>` with their object: the instance kernel binds 8 storage buffers, the core limit); one counter clear and
    3 dispatches a phase per page. The records' first instances, `identityBase` and the raster's bases are the page's;
    the culler's attributes are cells of its own following the page's buffers. Writes land before the commands of their
    submission (`StorageBuffer.write`): the rows go to one of 16 slots per run (`uRowBase`), a growth's copies and a
    freed range's zeroing are submitted first (`flush`). Per moving frame on the split: culling dispatches 1012 -> 7.7,
    clears 169 -> 1.3; whole rotating `off` 8.53 -> 5.11 ms (`noraster` 4.04), zoomed rotating `off` 15.9 -> 10.7 ms
    (`noraster` 11.8), `auto` 24.4 -> 19.9 and 30.0 -> 20.3 ms (the raster's per-object dispatches left); MG as one
    object unchanged (rotating whole 2.31 -> 2.07, zoomed 4.00 -> 3.99 ms).
  - **Step 6, the raster over all objects**: per object a row of what the kernels read (`SphereRasterRowVec4s`, its
    clip objects after the rows: the raster kernel keeps 8 storage buffers), uploaded per phase into one of 8 slots;
    the objects grouped by the buffers they bind (a merged group: the culling page and the same data pages; an
    unmerged object a group of its own), a group's rows consecutive. Per group one regions dispatch (a thread's object
    by a binary search over the rows' first threads; once per workgroup where all its threads are of one object, none
    for a group of one), one args dispatch for all groups (a thread per group: each object's first chunk and single in
    its group's dispatch, the group's indirect arguments and totals), one raster and one singles dispatch per group (a
    workgroup's chunk object by a binary search over the first chunks). Chunks stay in their objects' ranges, so the
    resolve is unchanged; the second phase uploads the rows again with the objects whose culler ran one flagged. On the
    split, per moving frame: dispatches 2378 -> 24, compute passes 340 -> 5, no bind groups created; whole still
    10.4 -> 3.1-4.6 ms, whole rotating 19.9 -> 4.4-6.3 ms (`noraster` 2.6-4.4), zoomed rotating 20.3 -> 5.5-6.8 ms
    (`noraster` 11.3-12.5; the runs' spread wider than before, the frames now CPU-bound: about 0.3-0.4 ms waiting on
    the GPU). Cost for one object: the threads read their object's parameters from the rows instead of uniforms
    (regions: 6 vec4s per thread, one per slot), MG's raster pass whole 0.52 -> 0.79 ms, zoomed 1.67 -> 1.87 ms,
    PetWorld in the Viewer zoomed 0.93 -> 0.96 ms (GPU, `--gpu-timing`); a workgroup sharing the row (a barrier and a
    copy per thread) was slower there (0.93 ms).
  - **Step 7, the CPU left** (profiled, whole rotating `auto`, 6.2 ms a frame without the profiler): by samples (the
    time deltas put a blocked main thread on the next sampled function), the frames' JS about 2.3 ms a frame: the
    opaque pass's two parts 45% of it (the resolves, a draw per object and part, 17%; the draw lists' replay 20%), the
    raster's preparation 19% (`getCapacity` and `readChunkCounts` 9%: a string key and a map per object, twice a
    frame), the culling's two phases 26% (the per-object JS: `isDispatched`, the row, `request`). Garbage collection
    28% of the frames' window, mostly the setup's (2 GB allocated in the run, 1 MB a frame by the frames:
    `getRasterInput`'s object per call, `setPlane` and `setHiZValues` per object in `request` where once per run would
    do, `setDispatched` and `setGpuView` building arrays per plane, the render items' bind group keys). Candidates, each
    small: the view's planes and pyramid set once per run, not per object; `getRasterInput` kept per dispatch; the
    raster's capacities by object (no string keys); the resolves of the objects without raster spheres in view left
    out (their chunk counts, read back a frame late).
  - **Step 13, option 1 as a budget per frame**: `gpuCulling: 'auto'` marks an object merged and culled on the GPU by
    its spheres, not its instances (`rasterOnly`, `Renderable.rasterSpheres`); the renderer culls such an object on the
    GPU only where the raster takes it (`Renderer.prepareCull`): never with `sphereRaster` 'off' or a view the raster
    does not draw, all with 'on' and 'only', the `sphereRasterObjects` (default 8) largest visible ones with 'auto'
    (hidden ones not at all). The others draw the CPU culling's lists, as while the kernels compile: no rebuild when
    the choice changes. Objects culled on the GPU by their instances (the membrane) and asked for (`gpuCulling: true`)
    are not counted. Split, median ms, same session (`auto` before -> after; `noraster`; `auto-objectsall`, the old
    behavior): whole still 12.1 -> 3.10 (2.36; 12.2), whole rotating 26.0 -> 4.20 (2.91; 19.5), zoomed still - -> 5.03
    (7.81; 9.96), zoomed rotating 31.0 -> 8.92 (8.63; 20.8); `off` now without GPU culling for the raster's sake: whole
    rotating 8.79 -> 3.66, zoomed rotating 18.0 -> 10.5. Budgets (whole rotating / zoomed rotating): 0: 3.97 / 11.1, 4:
    4.13 / 9.12, 8: 4.20 / 8.92, 16: 4.70 / 7.34. MG as one object unchanged (whole 0.95, rotating 2.47, zoomed large
    7.96 ms).
  - **Step 8, the pool removed**: the code as before step 5 (`CullBatch` and the bind group cache kept). Removed
    against the pool, the same session (median ms, the runs' level about 10% above the earlier ones): split whole
    rotating `auto` 6.22 -> 26.0, `off` 5.17 -> 8.79, `noraster` 4.36 -> 3.94; whole still `auto` 4.26 -> 12.1; zoomed
    rotating `auto` 7.36 -> 31.0, `off` 10.7 -> 18.0, `noraster` 12.0 -> 13.1. MG as one object faster again: whole
    `auto` 1.41 -> 0.94, rotating 2.98 -> 2.02, whole large 3.79 -> 3.32, zoomed large 7.76 -> 7.55. Tests, the
    native suite (but the two failing before) and the WebGPU visual profile (228) pass.
  - **Step 1, where the time outside the passes goes**: not the GPU. A Chrome trace (the GPU process' threads) shows
    its main thread decoding and recording the WebGPU commands (`WebGPUDecoderImpl::HandleDawnCommands`, `Queue::Submit`)
    5.2 s against 0.8 s `noraster` over the run, about 40 against 6 ms a frame. The commands per moving frame (counted
    in the page), `noraster` / `off` / `auto`: dispatches 8 / 1012 / 2378, compute passes 2.6 / 337 / 340, buffer clears
    1.3 / 169 / 172, bind groups created 58 / 64 / 373, draws 134 / 78 / 380.
  - **Step 4, bind groups of compute kernels**: a kernel keeps the bind groups of its last 4 binding sets, keyed by the
    resources' ids and versions compared as numbers (`createKernel`, `webgpu/compute.ts`): the raster's regions kernel,
    alternating between the records of the two culling phases, created 319 bind groups a frame, now 2. Whole rotating
    `auto` 30.3 -> 28.8 ms.
  - **Step 3, one compute pass per culling phase**: the objects' culling recorded into a `CullBatch`
    (`Renderable.prepareCull`, `prepareCullSecond`, `MergedCuller.dispatch`): their counter clears, then one pass of
    all dispatches (`Renderer.prepareCull`). 340 -> 5 compute passes a frame; whole rotating `auto` 28.8 -> 24.4 ms,
    `off` 9.36 -> 8.53 ms, the others within the runs' spread. What is left is per dispatch: about 2400 a frame.
  - **Step 2, the raster threshold** (`AutoGpuCullingMinRasterSpheres`, tried at 10M and 100M): 136 of the 150 objects
    have 100k spheres or more (95 at least 1M, 42 at least 10M; median 295 instances). Per frame, 100k / 10M / 100M /
    `noraster`: whole still 13.1 / 4.52 / 3.54 / 2.88 ms, whole rotating 24.4 / 8.91 / 4.00 / 4.04, zoomed still - /
    3.98 / 5.95 / 8.10, zoomed rotating 30.0 / 8.44 / 8.99 / 11.8. The raster pays for the large objects zoomed in,
    the many medium ones cost more than they gain; one object of a few million spheres gains (3J3Q, PetWorld in the
    Viewer), so the count of objects, not an object's spheres, is what decides: a budget per scene, the threshold
    unchanged.

## Where things stand

- **The trigger**: `gpuCulling: 'auto'` culls an object on the GPU from 50k instances, or, where the device has the
  sphere raster, from 100k spheres (`AutoGpuCullingMinRasterSpheres`, `resolveGpuCulling`, `repr/merged.ts`). The
  explorer's per-model objects have tens to thousands of instances of hundreds to thousands of atoms: without the raster
  one of the 150 is culled on the GPU (the membrane, 242k instances), with it most of them (149 draw calls still), and
  each also goes to the sphere raster (`sphereRaster: 'auto'` takes the objects culled on the GPU).
- **Per object culled on the GPU, per moving frame**: its culling in compute passes of its own (`MergedCuller.dispatch`:
  a counter clear, a pass of 3 dispatches; the second phase another pass of 3), the raster's 4 dispatches in each phase
  (`SphereRaster.draw`: about 40 uniform values set per object, its kernels' own bind groups) and a read of its chunk
  counts, its resolve draws in both parts of the opaque pass (scissored to its bounding sphere: the whole cell for an
  entity spread over it). A still frame dispatches the raster again (the culling only where the view changed).
- **Measured** (RTX 4070 Ti, WebGPU, 1280x720; ms per frame, CPU submit in brackets; without GPU timing, which adds a
  query per pass and distorts hundreds of small passes):

  | case                   | noraster    | off         | auto        |
  | ---------------------- | ----------- | ----------- | ----------- |
  | whole, still           | 2.88 (2.64) |             | 11.7 (6.29) |
  | whole, rotating        | 4.04 (3.78) | 9.36 (8.01) | 30.3 (15.5) |
  | zoom 0.3, rotating     | 11.8 (11.3) | 15.8 (10.5) | 30.6 (18.1) |

  With GPU timing (whole, rotating, `auto`): culling 1.59 + 1.58 ms (the two phases summed over the objects' passes),
  raster 2.08 + 1.85 ms, opaque 0.36 + 0.33 ms, 7.98 ms in all against 1.81 ms `noraster`; the frame's GPU time without
  timing is about 30 ms: most of it outside the passes (the passes' boundaries, the clears and uploads between them, the
  raster's reads), to attribute with a GPU trace.
- **CPU** (profiled, whole rotating `auto`, about 19 ms a frame): the kernels' dispatch preparation 37% (`getBindGroup`
  builds its key, an array and a joined string, per dispatch: 9% alone), the raster's `prepare` 25% (uniform values,
  `getCapacity`, `readChunkCounts`), the two culling phases 45% (with the raster's second phase inside), the resolves 9%.
- **The same model in the Viewer** (`mmcif_particles`) is one particle target object: one culling and one raster
  dispatch set per phase. The split, not the model, is what costs.

## Options

1. **The threshold for the raster's sake by objects, not spheres alone**: 100k spheres made the raster pay for one
   object (3J3Q, PetWorld as one object); a culled object costs a fixed 0.1-0.3 ms a frame. Require many spheres per
   object and per instance draw, or count the objects culled for the raster in a scene and stop taking more past a
   budget (a scene property, not a representation's: needs the decision where the scene is known). Quick relief, keeps
   the raster for the large objects; measure which of the 150 objects pay.
2. **One compute pass per culling phase for all objects**: the renderer opens the pass in `prepareCull` /
   `prepareCullSecond` and the cullers record into it; the counter clears before it (a clear is not allowed inside a
   compute pass: one clear of a shared counter range, or a clearing kernel). 2N passes to 2; the GPU's per-pass cost and
   the barriers between them go.
3. **Kernel dispatch CPU**: bind groups cached by the bound resources' versions without a string key per dispatch (an
   array compare, the entries built only on a change), uniform writes skipped where nothing changed. Generic: every
   compute kernel gains, not only these.
4. **Culling and raster dispatches over all objects at once**: per-object parameters in a table, the objects' buffers in
   shared pages (`shared-draw-resources.md`'s arenas): 3 + 4 dispatches per phase in all. The largest change, and the
   largest gain where options 2 and 3 leave the per-dispatch cost.
5. **Fewer objects in the explorer**: one representation for the entities of all models (colored and toggled per entity
   by a theme and by visibility per instance, not by representation), as the Viewer's particle target draws PetWorld as
   one object. An app change; the explorer's per-entity UI (visibility, colors, groups) would map to themes.
6. **The resolve per object**: its scissor is the object's bounding sphere, the whole cell for a mesoscale entity, so
   every object reads the visibility buffer over most of the screen; cheap here (0.36 ms both parts), to watch at large
   drawing buffers (a resolve per object restricted to the tiles holding its pixels).

**Decided first**: option 4, dispatches over all objects (many objects made cheap rather than avoided); the budget of
option 1 not for now. **Revised** after steps 5-7: see "Revised" below. Option 1 then done as a budget per frame in
the renderer (step 13).

## Design: culling and raster over all objects (option 4)

What a moving frame dispatches today per object culled on the GPU: the culling 3 + 3 (`cell-cull`, `instance-copy`,
`cull-args` in each phase), the raster 4 + 4 (`sphere-raster-regions`, `sphere-raster-args`, `sphere-raster` and its
singles), a counter clear; with 136 objects about 2400 dispatches and 170 clears. The target: 3 + 3 and 4 + 4 in all,
one clear, whatever the object count. The objects' kernels already share their pipelines (no defines in
`CellCullSchema`, `InstanceCopySchema`, `CullArgsSchema`, the raster's schemas): what differs per object is its buffers
and its uniforms.

- **A culling pool per context** (`MergedCullPool`, owned by the renderer or the compute context): one storage buffer
  per kind of the cullers' buffers, each object a range in it: cell spheres and levels, instance cells, counters, the
  args and second-phase records, the compacted attributes (`aTransform`, `aInstance`, `aSegment`) with their identity
  regions, the segment data. A range allocator with growth (a new buffer, the ranges copied, the objects' bindings
  updated), ranges aligned to 256 bytes where bound with an offset. The draws keep binding their attributes per object
  (`setVertexBuffer` with the range's offset, as the geometry pages) and their indirect records (the record's offset):
  no change to the draw programs. Past `maxStorageBufferBindingSize` a second pool (a dispatch set per pool).
- **An object table** (storage, one row per object culled this frame): its ranges (cells, instances, segments, records,
  identity base, segment data, clip base and count), its level ranges and count, its distance range, its raster
  parameters, its phase. Written per frame for the objects dispatched (one `writeBuffer`); the frustum, the camera
  plane and the pyramid are uniforms shared by all. The objects whose view did not change (`isDispatched`) are not in
  the table and keep their records (`compactRecords`).
- **The kernels by global index**: `cell-cull` over the cells of the table's objects (a cell's object by a binary
  search over the rows' cell bases, or a cell-to-row index written with the table), `instance-copy` over their
  instances (the instance's cell gives the row), `cull-args` over their records; each thread reads its row instead of
  the uniforms. The counters cleared as one range before the pass (scratch per dispatch: objects not dispatched read
  none of theirs).
- **The raster over all objects**: `sphere-raster-regions` over the records of the table's objects into the shared
  chunks buffer (it is shared already, an object a capacity range: `chunkBase`, `chunkCapacity`), the chunk holding its
  row; `sphere-raster-args` one set of indirect args for all chunks (a single range filled by an atomic, or the ranges
  kept and the args kernel summing them); `sphere-raster` and its singles one indirect dispatch each, a chunk's row
  giving its positions and sizes (the data pages: a page index and offset per row, the pages bound as a fixed set of
  bindings, a dispatch per page set where more), its size parameters, levels, clipping and transforms. The per-object
  `uRaster*` values of the resolve stay.
- **What stays per object**: the CPU's decisions (`predictRecords`, `setListRecords`, the level lists), the resolve
  draws (a pipeline and the object's uniforms each), the impostors' draws; unmerged objects the raster takes (`on`,
  `only`) keep their kernels.

Expected: about 2400 to 14 dispatches and 170 to 1 clears a moving frame on the split; the GPU process' command
handling with them (step 1: most of the frame). Unknown: the CPU left per object (`predictRecords`, the list building,
the resolves) once the dispatches go; the cost of the binary search or index per thread (small against the cell and
sphere tests).

## Revised: the pool removed, fewer objects instead

The culling pool (steps 5 and 6) made 150 objects cost about what one does, but at a price the graphics layer keeps
paying for every scene: pages that grow by copying (twice the memory meanwhile, a submission per growth and per freed
range) and never shrink, a pool with no owner that destroys it, objects no longer isolated (a growth rebinds every
object of a page, invalidating their draw lists; a stray write corrupts the others), runs over a page's whole capacity,
a hard limit of 16 runs a submission, and the parameters read from rows by every thread (MG as one object +0.27 ms of
raster). The many objects are the explorer's choice, not the data's: the same model in the Viewer is one object. So the
pool goes, and the explorer's scene gets fewer objects, by the app or by merging.

**Removal** (step 8): `merged-cull-pool.ts`, the kernels by row (`cell-cull.wgsl.ts`, `sphere-raster.wgsl.ts`) and
their schemas back to uniforms per object, `MergedCuller` and `SphereRaster` back to their own buffers and per-object
dispatches; `RangeAllocator.grow` and its test go with them (no other user). Kept: the bind group cache of compute
kernels (step 4) and the culling of all objects in one compute pass per phase (`CullBatch`, step 3), both independent of
the pool; the `mesoscale-split/*` cases. Expected after: as before step 5 (whole rotating `auto` about 24 ms on the
split), MG as one object back to its raster times.

**What the explorer does per entity** (inventory first, step 9): a representation per model (membrane and other
entities apart), its group (`MesoscaleGroup`) holding color (uniform per entity, generated per group), visibility
(`isHidden`, toggled per entity and per group), the graphics mode (levels, `lodPixelSize`, approximate), clipping and
focus. Each of these must keep working, and stay cheap to change, with fewer objects.

**Direction A, the app: one representation for many entities** (as the Viewer's particle path), with:

- **Color per entity** by a theme over instances (`instanceGranularity` is on already): an entity's color a value of
  the theme's per-instance layer, changed by updating the layer's range for the entity's instances, not by rebuilding
  the representation.
- **Visibility per entity inside a render object**, which Mol\* lacks: today a fully transparent group or instance is
  discarded per fragment (`check_transparency`), its vertices, culling and raster still paid. Two ways:
  1. **The transparency layer's 1 as "hidden"**, honored by culling: CPU culling (per instance, `cullSegment`), the
     GPU culler (`instance-copy`), the raster's chunk pass, the vertex stage (degenerate primitives). No new layer;
     but a transparency layer makes the object transparent (`dTransparency`, its passes and variants), so "all values
     0 or 1" must count as opaque with hidden instances (`transparencyAverage`, the draw facts), and the meaning mixes
     two things (a fade to 1 hides at its end).
  2. **A visibility layer of its own** (a byte, or a bit, per instance or group, as the marker and clipping layers),
     read by the same places; no change to transparency's meaning or passes, one more layer (memory: a byte per
     instance, 0.5 MB for PetWorld's instances).

  Hidden instances should also leave the instance grid's cells and the GPU culler's cells when a whole cell is hidden
  (a cell flag), so that hiding most of a model costs nothing per frame.
- **Clipping and focus** unchanged (they are per object today, and the one object takes them).
- What the UI loses: per-representation params (a group with another graphics mode than the rest); to check against
  the explorer's actual uses (step 9).

**Direction B, graphics: merging render objects across representations**: merged objects merge a representation's
units (`renderMerged`, members with their own layers and segments); merging the compatible objects of several
representations (the same geometry kind, program variant and material, by a key) into one merged object would serve
the explorer unchanged and any app with many representations. Each representation keeps its own render objects as the
members (their values, layers, updates), the merged object draws them; a member hidden is a segment left out (the
culling already skips segments: per segment records, `predictRecords`), a member's color its own layer. Costs: who owns
the merged object (the scene or a merging service, not a representation), membership changing as representations come
and go (the merged object rebuilt or patched: its layouts per member are sizes, so appending is cheap, removing needs
compaction), per-member uniforms that differ (the color uniform of a uniform theme: as a per-member value, like the
segment data), picking and marking by member. The dropped batched-draws work (`batched-draws.md`) did this for unmerged
objects and lost to merging; this merges at the level that won.

**Which first**: A's visibility layer is useful on its own (hiding parts of one large object: the Viewer's particle
targets, any merged object), and the explorer change then removes the many objects at their source; B is the more
general fix and the larger graphics change. Proposed order: the visibility layer (A2) measured on one object, the
explorer prototype on it (one representation for the non-membrane entities, color and visibility per entity),
measured against the split; B only if apps other than the explorer need it.

## Steps

1. GPU trace of whole rotating `auto`: done (the GPU process' command handling, see Status).
2. Option 1: measured (Status); the budget done later (step 13).
3. Option 2: done (Status).
4. Option 3: done for the bind groups (Status); uniform writes not yet.
5. Option 4, the culling (done, Status): the pool and its allocator (with tests: ranges, growth, release), the table, the three
   kernels by global index (GLSL not needed: GPU culling is WebGPU only), `MergedCuller` on pool ranges with its
   interface kept; the GPU culling tests (`webgpu-gpu-culling`, `webgpu-illumination-occlusion`, two-phase culling) and
   the visual cases unchanged; `mesoscale-split/*` and the single-object cases (`sphere-raster-mg/*`, the PetWorld
   viewer cases) measured, the call counts per frame again.
6. Option 4, the raster (done, Status): the regions, args and raster kernels by row; the raster's visual and headless tests; measured
   as step 5.
7. A CPU profile of what is left per object: done (Status).
8. Remove the culling pool and the kernels by row (steps 5 and 6), keep steps 3 and 4: done (Status).
9. Inventory of the explorer's per-entity operations and representation params (colors, visibility, groups, graphics
   mode, clipping, focus, picking and labels): what one representation for many entities must support.
10. Instance visibility in render objects: transparency's 1 as hidden (A1) or a layer of its own (A2), decided by step
    9 and a prototype; honored by CPU and GPU culling, the raster, the vertex stage; hidden cells skipped; tests (a
    hidden instance neither drawn, nor picked, nor rastered) and a benchmark (most of a model hidden).
11. The explorer on one representation per kind (membrane, others): color and visibility per entity by layers;
    measured against the split (`mesoscale-split/*` gains a case for it).
12. Merging across representations (B): designed if step 11 leaves needs it does not meet.
13. Option 1 as a budget: objects culled on the GPU for the raster only culled there where the raster takes them
    (`sphereRasterObjects`), the `off` setting honored: done (Status).

## Open questions

- Whether the explorer should keep a representation per model at all (option 5): its UI's per-entity controls (step 9).
- Visibility as transparency's 1 or a layer of its own (step 10).
- Merging across representations: worth it beyond the explorer (step 12).
- Draws of one merged object of many targets, measured: the MG cellpack model as particle targets (2.1M particles,
  848 targets, one merged object culled on the GPU) issues about 1,840 draw calls a canvas frame while rotating (about
  4 ms of CPU headless on Dawn), 13 once still (the pass replayed). The 15k first seen were the headless image's
  multi-sampling (`sampleLevel` 4, the scene drawn per sample), not the canvas: no question here unless a browser
  profile shows the rotating frames' draws.
