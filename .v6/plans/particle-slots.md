# Plan: Particle lists in slot pages

## Status

- Steps 1 to 5, 6a-6c and 6d1-6d4 done (below; 4 the basics, decided), 6d5 in progress, 6e and 7 open. The goal: a tile added or let go under 10 ms (the yeast cell in 8³ tiles, 30 in view: 676k particles,
  3,345 targets; 290-390 ms now, [particle-tiles.md](particle-tiles.md) Status, step 7). The same storage serves lists
  edited or made as they are used (procedural, editable): what [gl/instance-index.md](gl/instance-index.md) step 4's
  adds and removals wait for. What edits a list (a tool, its state, undo) is not this plan's: it takes what changed.
- Speed first (decided): the hot loops as fast as over flat arrays, an update's cost by what changed.
- **Step 1 done, the floor** (the yeast's 30 inner tiles in 8³, 3.4M particles, 4,880 targets; a tile of 140k
  particles of 2,640 targets added and let go; headless WebGPU, CPU): 510-560 ms an update. Per update: the view 162
  ms (its 3.4M particles concatenated, 47 ns each); the touched targets' colors 164 ms (their 3.4M instances, 48 ns
  each: the index theme a call an instance), their transforms and bounds 57 ms (17 ns an instance: the poses path);
  the merged object made again 82 ms; the signatures computed again 17 ms (every target's, an update); the rest
  (garbage, the state) about 50 ms. By the tile alone: the view nothing (its pages are the tile's own arrays: no copy),
  transforms and bounds about 2.3 ms, colors 6.7 ms at 48 ns (1.4 ms at 10 ns), the signatures nothing (a change of
  kinds only), the upload of 140k instances (about 10 MB: transforms, colors) 1-3 ms. **A tile touches thousands of
  targets** (2,640 of 4,880): a cost a target of a few µs alone is the goal; so a tile's instances one block for all
  its targets in the merged object (the GPU-culled one: its culler gathers a member's instances from the scene's cells
  a frame, so their storage order is free), added and let go as one, not a range a target (design 7 revised). With it
  the floor about 4-8 ms: the goal possible, with the colors at about 10 ns an instance (a theme's colors by target
  where they are, an index theme's by key in a pass over the tile, or computed as the tile loads) and the upload of a
  block. WebGL (no GPU culling): a range a target, its floor by the targets touched.
- **Decided**: `count` the slots (`[0, count)`), the pages optional: a flat list as it is, every consumer of one
  unchanged (no break), a slot list read only by the consumers that read pages (adopted one by one, step 4); pages of
  any size (a page a slot range with a flat list of its own: a tile's page its list, no copy, no waste; a slot's page
  found by a search of the starts, the hot loops a page at a time). A tile view's slots its set's particles, a tile's
  page at its base: a particle's slot its key.
- **Step 2 done**: `particle-slots.ts` (`ParticleSlots`): `ParticleList.pages` (the pages by start, their starts,
  the live count) and `pagesFrom` (the version before's pages and the pages added and removed, by identity: one step);
  `create` (a slot list of `count` slots, its pages checked: within the slots, none overlapping, each flat; what
  they share the info's: a catalog's), `update` (pages added and removed, the next version), `getChanges`,
  `forEachPage` (a flat list one page), `findPage`, `locate` (a slot's list and index, none for a gap),
  `getLiveCount`, `flatten` (a flat list by slot, its keys the slots: for what reads flat lists only). A slot list's
  own per-particle fields empty, its labels by slot. The concatenation of lists' particles shared with the tile view
  (`particle-concat.ts`). Tests: `particle-slots.test.ts`.
- **Step 3 done**: `ParticleTileView.createSlots` (a slot view: a tile's list its page at its base as it is, a
  particle's slot its key; the next version of a previous one of the same set, a page kept where its tile's list is the
  same; the targets the source's loaded ones, a loader's); a change's `posed` pages (a page replaced by its particles
  posed otherwise, `Particle.isSamePoses`: a frame step, neither added nor removed); a slot view no composition
  (`isSamePart` false: its consumers by its changes). Tests: `particle-tile-view.test.ts` (a page a tile, its slots
  its keys, the particles the concatenated view's; tiles added and let go, the others' pages kept; a frame step all
  posed), `particle-slots.test.ts` (posed). Measured on the yeast's 30 inner tiles (3.4M particles): a tile added or
  let go 0.26 ms (the most of it the targets' map copied, 4,880), the concatenated view 188 ms. The plugin's view still
  the concatenated one until the representations read pages (step 4).
- **Step 4 done, the basics** (decided: the target path natively in steps 5-6, no interim one): `findPage` the page
  found last tried first (slots in order without a search); `Particle.position` (a slot's by its page; a gap the
  origin), `getBoundingSphere` (a loci of slots by their pages' particles), `getBoundary` (a slot list's of its
  pages' boundaries, each cached on its list: a version's by its pages), `remapLoci` (a loci's slots kept across
  versions, a let-go page's dropped), labels by slot; the five particle color themes on slot lists (graphics
  `theme/color/particle-pages.ts`: a location's page, a picker a page list made once; each theme's choice of it once,
  so a flat list's colors as before): entity (the palette by the shared entity info, a page's index its own), index (a
  slot of the slots: a tile view's key of the set's count), attribute (a page's own column, the first page's range),
  compartment (by the shared compartment info, else the pages' in their order), hierarchy (its compartments and their
  entities those the pages hold, a page's index by them; the structure collected over lists and the index a list's,
  split). Since 6d5 both themes color in the catalog's order, flat and slot lists alike: compartments by
  `compartmentInfo`'s, entities by `entityInfo`'s, the others after them as met; the hierarchy's compartments and their
  entities the catalog's where its info holds them (`ParticleCompartmentInfo.entities` and `withoutEntity`, filled
  from the whole list by `ParticleTileCatalog.fromList`), so a view of some tiles, flat or of slots, colors as the whole
  (MG_1184 in tiles: the views' colors had differed). Spacefill, orientation, fibers and targets, and the move and attach-targets transforms not applicable to a
  slot list yet. Tests: `particle-slots.test.ts` (a location, bounds, a loci across versions), graphics
  `theme/color/_test/particle-slots.test.ts` (the themes applicable as to the flat equal; entity, attribute and index
  colors a slot's as the flat view's particle's; compartment and hierarchy by what a particle is, compartment the same
  in a view of fewer tiles; a gap no particle). Measured, flat lists (the yeast, 31.5M): the hierarchy theme's build
  183-196 -> 174-198 ms, 4.5M of its colors 109-118 -> 101-113 ms, of the entity theme's 100 -> 94 ms: no slower.
- **Step 5 done, colors by the tile**: `ColorTheme.fillParticleColors` (a flat list's colors, a slot list's page or
  the theme's own list, as RGBA all at once); graphics `particle-pages.ts`: a painter a page (its colors and its fill,
  made once a page list), `fillByIndex` (a theme's color index: runs a run at a time, the texel doubled; an array a
  word a particle), `fillParticleColors` (a theme's fill, else a location each), `fillByLocation`. The five themes
  fill (entity, compartment and hierarchy by their index, attribute and index a particle each); the index theme's
  interpolated palette a table of 65,536 steps over more indices than that (a slot list's millions: made once for the
  palette and count, RGBA words; a channel within 0.02 of the exact); the hierarchy's group colors made once for their
  arguments (distinct colors in Lab: 100 ms a theme for the yeast's, flat lists too), its context hash a slot list's
  structure (another, the colors all again). Tests: graphics `theme/color/_test/particle-slots.test.ts` (a page's fill
  its slots' colors, a flat list's its particles', the table's fill its locations'). Measured, the yeast's tile of
  140k in its 31 inner tiles: entity 0.29 ms (2.0 ns a particle; its page's index first 1.36 ms), compartment 0.13 ms,
  hierarchy 0.22 ms, index 0.49 ms (3.5 ns; the table first 1.9 ms), against 6-7 ms a location each; the hierarchy
  theme made again 169 -> 9 ms (its structure collected over the pages: a version's; the open questions).
- **Step 6a done, the store**: graphics `target/tile-store.ts` (`ParticleTileStore`): the members a pool (a target its
  member while it has instances; a freed member's index a new target's), a tile a block of the instances by member,
  an instance its member (`segments`; `DeadMember` in a gap), its particle's slot and its pose (`PoseFloats`:
  `t = p - s R c` and `s`, the unit quaternion: `CulledData`'s slots); a block let go a gap, the last's the store
  shorter, compacted past a share (the blocks moved down in order); the arrays grown by half again, what is there
  copied; `takeChanges` (the ranges written, the members touched and freed, grown, compacted). The poses written by a
  function of its own (every input a parameter: the JIT's monomorphic loop; the positions read inline); a slot list's
  (`gpuPoses`, 6d5) composed on the GPU from a block's runs instead. Tests: graphics
  `representation/_test/tile-store.test.ts` (a block by member, its slots and poses placing its targets; a gap, a
  member freed and its index reused, the last block's; a compaction; a growth). Measured, the yeast's 31 inner tiles
  (3.4M instances, 4,860 members): the 30 in 81 ms; the tile of 140k (2,640 targets) added 2.6 ms (18.6 ns an instance:
  its arrays cold, a run's lookups of 2,640), let go 0.08 ms. Trimmed later where the goal needs it: a tile's grouping
  by target made once (its runs and their members' offsets).
- **Step 6b done, the regions apart**: `MergedSegments.regionBases` (a segment's first slot in each culling region,
  apart from where its data is: a store's member's instances in many blocks, its region a range of its live count;
  `instanceBases` without), read by the culler's segment data and records (`regionBaseOf`); the merged objects as
  before (the WebGPU culling suites: culled poses, GPU culling, two-phase, sphere raster and its split, 32 tests). The
  gaps (proposed, 6d): a segment of their own (the dead instances' member), its draw count 0 and its sphere empty and
  far, culled with its cells before any copy: no kernel changed; the regions sized by all instances (`instanceCount`
  the store's used, gaps with), the members' bases a prefix sum of their live counts: the gaps the regions' room.
- **Step 6c begun, the culled data by capacity**: `CulledStore` (`culled-data.ts`): a store's `CulledData` laid out
  by its instance and segment capacities (the header's count the capacity: only the sections' offsets read it, the
  kernels their count by uniform), so a block added or let go is its ranges written (`writeInstances`: slot words, an
  instance its own two pose vec4s, `aInstance` its index, its segment its member's or the dead segment, its poses
  copied as the store holds them), the segments' data by range (`writeSegments`), the header (`writeHeader`); the
  store grown, all again. Tests: `gl/_test/culled-data.test.ts` (the readers' rows, instance and segment of each; a
  block let go, its range alone written, its instances the dead segment's).
- **Step 6c done, the layers**: a store's per-instance layers (colors, sizes, markers) plain arrays of its capacity, a
  block's texels `[first, first + count)` written and logged (`TextureChanges.update`: the render item writes the
  ranges logged alone, a layer of a fixed size kept in its range), all again where the store grows; `fillBlockColors`
  (`tile-store.ts`: the theme's fill of the page at the block's first texel where the block is in its page's order,
  `TileBlock.inPageOrder`: a pack's tile, a target a run; else into a scratch and gathered by slot). Not grown by
  copying (decided): a buffer's write runs before the submission's commands, a copy in their order, so a copy after a
  growth would overwrite what was written; and a store grown moves its layout's sections, all written again anyway: the
  culler's `ensure` (made again, not copied) as it is, the store's ×1.5 the room. Tests: `tile-store.test.ts` (a
  block's colors in its page order and gathered, an instance's its particle's, the others' texels untouched).
  Measured: the yeast's tile of 140k in page order, its colors 0.48 ms (the index theme; its table first 2.0 ms).
- **Step 6d1 done, the merged hooks**: `MergedInstanceSource` (`merged.ts`; `createMergedValues`'s last argument):
  a member's instances its `instanceCount` (the store's live count; the gaps' member the gaps'), so `computeSegments`
  lays the segments out as the store's (a region its live count's, the instance bases their prefix sums and the culler's
  region bases with them, the gaps' after all, the instances the store's used); its layers of instance granularity
  its cells (the merged object's own: their ranges logged reach the render item), not merged from the members;
  `rebuildSourceInstances` the segment table alone (no transforms, poses nor `aInstance` of the members: `aInstance`
  serial, the culler's data the store's); an instance count changed alone not a new layout (`instanceCountsChanged`:
  the segments, the table and the counts again, `relayoutSourceInstances`, the geometry and layers kept). Tests:
  `gl/_test/merged-source.test.ts` (the segments by the counts, the layers the source's cells, no instance data; a count
  changed alone the segments again and the geometry's arrays the same; without a source merged again). The graphics
  suite and the WebGPU culling suites as before.
- **Step 6d2 done, the culler's hooks**: `MergedInstanceSource.culled` (its capacity, a data version, `write`: all or
  what changed, the gaps' segment the last member; per segment its largest scale), `Merged.source`; the culler's
  `syncSourceInstances` (`merged-culler.ts`: laid out by the source's capacity and the segments, `CulledStore`; the
  header and the segments' data where the layout changed; the instances all where the buffer was made again or the
  capacity changed, else the source's ranges), the source's version in the instances' versions, neither a move nor
  poses alone with a source (its members' grids and transforms not read), the segments' scales the source's;
  `TileStoreCulled` (`tile-store.ts`: the store's changes noted, written by range at the culler's next sync, all where
  it grew or compacted, nothing past its used). Tests: `tile-store.test.ts` (all at first, then the noted ranges, a
  last block let go none, grown all). The culler's own path as before (the graphics suite, the WebGPU culling suites:
  35 tests); the source's path drawn in 6d3 (its test there: a headless render of a tile view).
- **Step 6d3 done, the representation's object**: `target/slots.ts` (`createTargetSlotsObject`; the target
  representation's for a slot list where it culls on the GPU, else the list flattened, picked by its keys, not marked:
  6e). The store's pages by the list's changes (`ParticleSlots.getChanges`: removed, posed again, added; all again where
  not the next version, the props, the size theme or the targets changed, the targets by content: a slot view's map is
  its own); a member a target while it has instances (`createTargetMember`, `common.ts`: its geometry and sizes as a
  visual's, one placeholder instance, colors of instance granularity; kept by target while the props hold), a freed
  member's values kept (no instances); the gaps' member last (`createEmptyMergeMember`, `merged.ts`: no geometry, its
  bounds' cells the first's); the object made again where the members are others (`createMergedRenderObject`'s
  `source`), else its counts set. A member's bounds the store's box of its particles (`memberBoxes`, grown, never
  shrunk) padded by its geometry's reach from the target's center. The colors and markers layers of the store's
  capacity, the added blocks' colors by range (`fillBlockColors`), all where grown, compacted or the theme changed;
  an added block's markers cleared. Picking: an instance's id its store index (`aInstance`), its slot
  `store.slots`, a gap none; marking by slot (a block's offsets of its page's particles where not in page order), the
  layer's markers by `Visual.mark`, their aggregates the source's (`MergedInstanceSource.markers`; no marked lists).
  Only impostor targets (a structure's or sites' spheres) held (`TileStorePlacement.accepts`; one geometry type an
  object), the others left out. `ParticleTilesView`'s `slots` makes the view a slot list (`createSlots`). A source's
  members by their levels on the GPU and the CPU alike (since 6d5): the culler's records a level's items
  (`layoutRegions`) and the segments' strides and scales (`writeSegmentLod`, `getSegmentData`) whatever their grids,
  a member's grid of cells (`instanceCellSize`, `hasInstanceCells`; `getTargetCellSize` by the store's box, a
  quarter's hysteresis) for a CPU cull; a stride over the whole count had read past a member's spheres (stray spheres,
  blobs: MG_1184 in tiles), and drawing all whole cost the GPU 300-1,100 ms a frame (the yeast at 4³).
  Test: `plugin/headless/_test/webgpu-particle-slots.test.ts` (headless WebGPU: a slot view of all tiles the flat
  view's image, tiles let go and added again the same with the object kept, a pick a slot, marked and cleared). Without
  levels (`lodPixelSize` 0): an instance's level is its cell's, and the store's one grid and the members' grids
  differ (an impostor drawn at another level, slightly other).
- **Step 6d4 done, the CPU-culled views and the gaps**: `MergedInstanceSource.getMember` (a member's instances' indices
  in the source and their poses, `InstancePoses` with `t` the position and the center 0); a source member's grid
  pending (`createTargetMember`: `createPendingInstanceGrid` of the props' cell and batch sizes), built from those
  poses where a CPU cull first needs its cells after a change (`ensureIndexes`, `merged.ts`: `cullVersion` and the
  source's data version), the identity region a member's slots' source indices (in its grid's order; the gaps' member
  in its own, no grid: `createEmptyMergeMember` a grid of none). The store's side: `ParticleTileStore.gatherMembers`
  (a pass over the blocks' member tables, the gaps by a pass over the instances), gathered once a change
  (`slots.ts`). The gaps: a record's spheres its index count (`sphere-raster.wgsl`), the gaps' member's 0: neither
  the impostors nor the raster draw them (their instances culled and compacted all the same: a cost while gaps last).
  Tests: `tile-store.test.ts` (gathered by member, the gaps last, a pose placing a point as the store's);
  `webgpu-particle-slots.test.ts` (targets of two shapes, a tile of eight let go and kept as gaps: the GPU's and the
  CPU's cull, the latter from a view of its own, each the flat view's image; the raster alone the same; the serial
  identity region, as before, fails it). Open: a member's cell size fixed (the props'), not by its count and extent
  (`getTargetCellSize`, the flat path's).
- **Step 6d5 in progress, measured and the object kept**: the yeast cell (31.5M) at 4³, 8³, 12³ and 16³ tiles, a
  central region of about one volume (0.7-1.0M particles, ~4,000 targets), four tiles next to it added and let go
  twice (headless WebGPU; a tile's molecules fill 3.4-6.2% of its box at every size, as the cell's 5%). First every
  add made the object again (a slot view's target map has new ids as tiles come: taken as other targets, all made
  again, 300-1,100 ms; then the pool's members by their targets: the object made again, 174-253 MB uploaded, frames
  110-1,000 ms and Dawn's `submit` stalls after). Now:
  - **A pool of members** (`slots.ts`): records (`MergeRecord`: the cells a fill writes their own, the rest the
    empty record's; the shared values once), filled with a target's values (`fillMergeRecord`), the room `PoolHeadroom` past the store's members and at least a member a
    target of the list's catalog (`targetInfo`: a tile set's, the object made once; made again past it); the store's members by target (`ParticleTileStore.memberCapacity`: a new target a member never used
    while there is room, else the longest free; a target let go and taken again its own member back, its values
    as they are).
  - **A source's geometry by allocation** (`merged.ts`, `Merged.replaceGeometry`): a member filled with another
    target's geometry allocated at the arrays' ends, the others not moved; the arrays extended
    (`SegmentedArray.extend`: the same parts first; the uploaders write the parts added, `writeChanged`,
    `changedRanges`), longer where the room is used up (`sourceRoom`: ×1.5, at least 16,384 elements more; the
    parts where they were); laid out whole with room where the unused draw more than the used.
  - **A count change alone**: the segment table's rows of the members replaced alone, the bounds a box of the
    members' spheres (no boundary helper).
  - **The gaps out of the cells**: the scene cells' bin (`scene_cell_bin_wgsl`) and the instance cull skip the gaps'
    segment (`uSkipSegment`), the scene grid's box without it: their stale poses had clamped into border cells, a
    cell's atomics all theirs.
  - **The passes over the members**: the culler's layout key a source's `cullVersion` (not two passes a sync), its
    layout's reads of the members by cells resolved once (`layoutCellsOf`; a values object's named reads slow), no
    move bookkeeping for a source, the segment data once; the merged object's counts and spheres by `memberCells`,
    its aggregates again on a geometry change alone; the pool's counts' cells and bounds' spheres kept.
  - **The members changed alone**: a member's cells cold (scattered over the heap: a pass over ~6k members 1-3 ms,
    three times its warm cost), so an update reads those it changed, not all. The source tells which
    (`MergedInstanceSource.takeMembers`: the store's changed members, the gaps'); the merged sync reads those alone
    (its layout check, version loop and snapshot), the segments and bounds from typed copies (`layout`, the spheres),
    and counts a member's changes (`Merged.memberChanges`); the culler keeps per segment what its layout reads (the
    segment data's rows, draws and levels' counts, bounds; `refreshSegments`), read again where the count moved, the
    first slots and records from them. The pool's counts and bounds of the store's changed members alone, a bounds'
    cell where its sphere changed (a box not shrunk as particles go), a target's acceptance (`accepts`: its geometry's
    choice) kept with `byTarget`, not made again an update (2.7-12 ms of an add at 8³). Fixed region bases (a count
    change no relayout) not taken: the records and segment data were ~2 ms of it, and a member's headroom used up
    every few adds a relayout whole again.
  - **The counts and bounds as copies**: a tile touches about a third of the targets (~2k of yeast's 5.6k), so the
    changed members alone still cost three passes over cold cells (the pool's, the merged sync's, the culler's) and
    the records' prediction a fourth over all (`predictRecords`). The pool keeps a member's count and bounds in
    typed arrays (`MergedInstanceSource.counts`: its cells written with them, the CPU cull and the marking read
    those; a paged source's counts in the copies alone since, the merged object and its CPU cull reading those:
    [particle-tiles.md](particle-tiles.md) 10a4), its reach and grid size read once a fill; the merged sync reads the members whose count or bounds alone
    changed (`counts.take`) from the copies, their cells not (not counted in `memberChanges`: the culler's rows
    stay), the culler its cells' and the prediction's bounds from them (no transient transform).
  - **The poses composed on the GPU**: the store's write ~20 ns a particle in V8 (the decode, the normalization and
    the pose's composition; a bare loop as slow), so a slot list's store (`TileStoreOptions.gpuPoses`) writes an
    instance's member, slot and its member's box alone, and keeps a block's runs (`TileBlock.poses`: a target's
    particles, their first slot, center and scale; those scaled by radius by particle apart). The culler's write of
    the store (`MergedInstanceSource.culled.write`'s `expand`) composes the blocks added (all after growth or a
    compaction) by the flat path's kernel (`culled-poses.ts`, `expandBlocks`): a block's values uploaded as stored
    (a pack's 20 B a particle where a slot is 32 B), a record a run written in place (thousands a tile: an array of
    numbers cost what the composition did). `gatherMembers` (the CPU cull) composes a block's poses on the CPU, once a
    change. The store's write 2-3.5 ms at 8³ (of 3.5-5.5): the rest the runs' own cost (~500 ns a run: the
    placement's lookups and maps, a run a target in a pack's tile).
  - **Fewer submits, fewer reads**: a submit ~0.4 ms of CPU (Dawn), so the poses' last batch (`expandBlocks`) and the
    identity region's run (`U32RunWriter.write`'s `submit` false) go with the next submit (a frame's; a buffer's
    destroy deferred to it already), a CPU write of the identity region submitting first; a member's center kept by
    the store (asked once a member, not once a run); a count change alone the segments' instance bases, counts and
    group offsets from `layout` (`computeSourceCounts`), not the allocations' records again.
  - Measured (CPU of the state update, the merged and the culler's syncs in it; the first frame; the GPU's by
    `onSubmittedWorkDone`): a let go 4.5-9 ms at 8³ and 5-7 ms at 16³, an add of known targets 8-14 ms and 5.5-7 ms,
    a first add (new targets' geometry made) 22-85 ms and 13-41 ms (4³: 140-620 ms), frames 12-18 ms on the CPU and
    11-23 ms on the GPU (the flat view's 5-48 ms). The stalls of 250-1,300 ms: the GPU drawing all whole (above), not
    the layouts. An add at 8³ (a profile of four): the slot object's update ~4 ms (the store's add ~2.5: its runs
    and `memberFor`; the colors ~0.4), the culler's sync ~3.3 ms (the poses' records and the instances' words ~1.4,
    the regions' records ~0.75, the segments' data ~0.5, the cells ~0.5), the merged sync and the state ~1-2 ms. A
    first add's: the new members' geometry and its allocation (the pool's fill, `appendSourceGeometry`, the
    culler's layout). Tests: `merged-source.test.ts` (the members a
    source tells of read alone, their changes counted; none read nothing, undefined all; the counts and bounds of
    those it counted from its copies, their cells not read), `segmented-array.test.ts` (an extension's parts written
    alone, a line of them, another array all), `tile-store.test.ts` (a member by its target; with `gpuPoses` the
    same members, slots and boxes, a block's runs, the poses gathered as the CPU composes them, the culled write's
    blocks from where they are), `webgpu-culled-poses.test.ts` (a store's blocks' slots as its CPU composes them, in
    batches: a run split across them),
    `webgpu-particle-slots.test.ts` (a new target filling a pool member in an object drawn before: the extension's
    parts uploaded, the image the flat view's; skipping them fails it).

## Why

- **An update costs by the whole view, not by the tile**: a view concatenates its tiles by target, so a tile added
  moves the indices of what follows it and the list is made again (45 ms); the representation visits every target
  (3,345); a target the tile has particles of builds all its instances again (95 ms); the one merged object of all
  targets is made again where a target comes or goes (105 ms). A pass over the view alone at 10 ns a particle is 7 ms.
- **What a particle is kept by is its index**: selections, marks, picking, labels, a theme's colors read the index;
  indices that move with every load make each of them a remap (step 7's kept parts: a target's indices taken again).
- **Lists that change as they are used** (particles added, removed, made procedurally) need the same: an index kept
  while the particle lives, storage that grows and shrinks without moving the rest.

## Design

1. **Slots**: a particle's index is its slot, kept while it lives; a list's slots `[0, capacity)`, those in use its
   live ranges (an `OrderedSet` of slots, few ranges: a tile's a range). `count` the live slots; a gap is no particle.
   A flat list (every list today) is one live range `[0, count)`: nothing changes for it.
2. **Pages**: per-particle data in pages of `2^k` slots (one size a list, e.g. 2^16): slot `s` in page `s >> k` at
   `s & mask`; a page a typed array a field (positions as now, float32 or on one grid; rotations float32 or snorm16;
   values arrays or runs). A flat list's page its whole array (`k` its size): its arrays as they are, no copy.
   Pages only where slots live: memory by what is loaded, never by the whole (the yeast's 31.5M slots: its pages for
   the tiles loaded).
3. **Loops by ranges**: what reads all particles walks the live ranges a page at a time (`Particle.forEachRange(list,
   (start, end, page, offset) => ...)`: contiguous typed arrays, a hot loop as over a flat array); a lookup a particle
   (`page >> k`) only where one is wanted (picking, a label, a location). No per-particle indirection in a hot loop.
4. **Tiles in slots**: a tile set's slots are its tiles' ranges, each aligned to pages (its base rounded up: a tile
   its own pages, freed whole); a view the set's slots with the ready tiles live. A tile loaded fills its pages, let go
   frees them; no other particle moves. Its slot is not its key: the key (its index in the whole, `ParticleList.keys`)
   is what it is, the slot where it is stored.
5. **Changes**: a list's next version says what changed since the one before (slot ranges added, removed, posed; the
   targets they are of): what its consumers update by, no diff of the whole (replaces the view's composition check,
   step 7). A version's line (`changesFrom`) as `posesFrom`: two versions of one line, the changes between them.
6. **Editing**: allocate (a free list of slots and pages), free (a gap left), compact (rare: the slots renumbered, a
   remap with it for what keeps indices). The editable and procedural lists' storage; their tools, state and undo are
   theirs.
7. **Representations by changes**: a tile's instances one block for all its targets (step 1: a tile touches
   thousands of targets, a cost a target alone would be the goal): in the GPU-culled merged object, whose culler
   gathers a member's instances from the scene's cells a frame, a block of the tile's instances (each its member, its
   transform, color and size) added as one pass over the tile and one upload, let go as a gap the culling skips
   (compacted when many); a member's bounds and counts by the blocks it has instances in. Members (targets) added where
   a tile brings a new one, their geometry as now. Without GPU culling (WebGL, a few targets): a target's instances by
   slot range (the instance store of [gl/instance-index.md](gl/instance-index.md)), merged objects edited in place,
   buffers with room to grow (×1.5).
8. **What the whole gives**: statistics read the set's or the catalog's (a palette's count the set's, an attribute's
   range the whole's), not the live particles', so they do not change with a load.
9. **Gaps never picked or kept**: picking, loci, labels and marks by slot; a change's removed ranges dropped from
   what is kept (a selection's particles gone with their tile).

## Steps

Each measured against the goal (the yeast cell, 8³ tiles, 30 in view, a tile added and let go: headless WebGPU,
`particle-tiles.md` step 7's case).

1. **The floor** (done): the current update taken apart for one tile (the view, the visit of every target, a touched
   target's instances, the merge), and what each would cost by the tile alone: what the goal leaves for each.
2. **The model** (done): slots, pages, live ranges (`ParticleList` with them; a flat list one page), `forEachRange`, a
   lookup, changes and their line. Tests: a flat list read as before; a list of pages read as its flat equal; a tile
   added and let go, the others' slots kept. Not breaking (decided): `count` the slots, the pages
   optional, a flat list as it is.
3. **Tile views in slots** (done): the view a slot list, a tile's pages filled and freed; its changes. Measured: a view's
   update by the tile (the goal: about 1 ms).
4. **The consumers on ranges** (the basics done; the representations natively in 5-6): the particle representations (spacefill, orientation, fibers, targets), themes, loci,
   labels, export (the pack writer compacts). Measured: no slower on flat lists (the particle benchmarks, the yeast's
   load and first frame).
5. **Colors by the tile** (done): a theme's colors by target where they are (one a member), else a pass over the tile at
   about 10 ns an instance (the index theme by key), or computed as the tile loads. Measured: 140k instances.
6. **A tile-blocked target object** (WebGPU, GPU-culled; decided: a merged object is of one representation, so the
   target representation builds its own for a slot list, not a merge of a render object a target). What there is
   (surveyed): the culler's kernels go an instance at a time, its member read from the culled data (the cells, the copy
   into regions: no member's instances need be contiguous; the workgroup fast path where neighbours share one); the CPU
   side by member ranges (`instanceBases` a member's data, its region slots and its picked ids; a count changed
   rebuilds all: `updateLayout`, `rebuildInstances`, the culler's `syncLayout`); buffers grown to the size, made
   again without their contents (`ensure`); a member's textures parts of one layer. Its parts:
   - **6a, the store** (done; CPU, no GPU): the members a pool (a target its member while loaded: its geometry, levels and
     bounds once); a tile a block `[first, count)` of the instances, by member within it (the fast path), an instance
     its member, pose and slot; a block let go a gap (its instances a dead member, drawn by none), compacted past a
     share; the arrays with room to grow (×1.5), what changed their ranges (the uploads'). Tests and measured: a tile
     of 140k added and let go (the goal: about 1 ms).
   - **6b, the regions and records** (done: the bases apart): a member's region base apart from where its data is (a prefix sum of its live
     counts with room), the indirect records of the members touched; the cells built again on the GPU (as a pose
     change's). Measured: 2,600 members touched.
   - **6c, the layers and buffers** (done): colors, sizes and markers one layer each in block order, by sub-range
     (`fillParticleColors` a block), the culled data's and layers' buffers grown by copying (`copyTo`), not made
     again; the dead members' cells never drawn.
   - **6d, the object**: the target representation's own for a slot list on WebGPU (its members' geometry as now, its
     instances the store's), updated by the list's changes; picking by the block table (an instance its member and
     slot). Measured: the update whole against the goal (the yeast's 31 inner tiles, a tile added and let go). How
     (surveyed; decided: a source hook, not a renderable of its own, which would repeat about 600 lines of the lists',
     the culling's and the raster's): `createMergedValues` and `createMergedRenderObject` given a
     `MergedInstanceSource` (GPU-culled only): the segments from it (a member's count its live count, the region bases
     their prefix sums, the instances the store's used, the dead segment the gaps'), `rebuildInstances` its segment
     table alone, `updateLayout` and `sync` by its layout and data versions (not the members' instance cells), the
     per-instance layers its arrays (by range), the marker aggregates its; the culler's instances by `CulledStore` by
     range, its versions in the keys, the move and pose shortcuts not taken, the segments' scales from it; a member
     pool (the store's segment capacity and the dead one; a freed member drawn of none, `drawCount` 0, in its place;
     past it the object made again, ×2), the representation's own object, not `mergeRenderObjectsByType`'s; picking an
     instance's member and slot (`segments`, `slots`; a gap none). The CPU-culled passes (decided): an instance order
     grouped by member from the store (`setIdentityOrder`), a member's grid from the store's poses where a CPU cull
     first needs it (paid once a change, the GPU views fast). Its parts:
     - **6d1, the merged hooks** (done; `merged.ts`): the source, tested without a GPU (a fake source).
     - **6d2, the culler's** (done; `merged-culler.ts`): tested natively (WebGPU).
     - **6d3, the representation's object** (done; `target/slots.ts`): the member pool, a member's geometry, the
       colors a block, the updates by the list's changes, picking and marking.
     - **6d4, the CPU passes and the gaps** (done; `merged.ts`, `tile-store.ts`): the grouped order and the lazy
       grids; the dead segment kept out of the sphere raster.
     - **6d5, measured**: the yeast's update whole against the goal.
   - **6e, the fallbacks**: the CPU-culled views (the image pass, illumination's thickness pass, kernels compiling):
     their structures built when needed (slow, measured, then by the store if they matter); WebGL: a slot list drawn as
     its flat copy (`ParticleSlots.flatten`: correct, slow), its own path later (`gl/instance-index.md` step 4).
   - **6f, a paged GPU store** (prototyped; the design toward it, not begun). Why: 6d5's update is spread over every
     layer (the store's runs and its per-instance arrays, the pose records, the instances' words, many small uploads:
     `writeBuffer` an eighth of a browser profile, the colors a particle at a time on the CPU), so an add stays at
     11-16 ms at 8³ and a view's worth of tiles blocks for a third of a second; no layer of it is left in a store made
     for streaming tiles.
     - **The prototype** (raw WebGPU, local, the yeast on Dawn): pages of 2^16 particles in two buffers, a tile's grid
       words and snorm16 rotations uploaded as stored (a write a page each, no pass over the particles); a page's record
       (the grid's origin and step, its runs), a runs table (first, count, target, radius); one kernel a particle: its
       run by a binary search of its page's, decoded, culled by the frustum at its target's radius, its pose composed
       (`t = p - s R c`, the unit quaternion) into its target's region (a counter a target; the regions the prefix of
       the targets' live counts, on the CPU). Measured: an add 0.3-1.4 ms at 8³ (80-160k particles, 10 writes) and
       0.04-0.8 ms at 16³, a let go 0.03-0.4 ms (its pages freed); the cull and compose 0.9 ms for 1M visible and
       1.4 ms for 6.4M (the first frames up to 26 ms: the GPU warming), the poses the CPU's within 1e-3 Å. In Chrome
       (the same page): an add 0.8-2.1 ms at 8³ and 0.2-0.6 ms at 16³, a let go at most 0.3 ms; the cull and compose
       0.2 ms for 1M visible and 0.7-0.8 ms for 2.5M of 6M, but up to 10 ms with all 6M visible: the counters'
       contention (an atomic an instance on its target's; Dawn alike, 1.8-10 ms). A workgroup's counts first (its
       particles consecutive in one page: a few runs, counted in workgroup memory by their index from its first
       thread's, one global atomic a run, more than 64 runs the global one an instance): all visible 1.6 ms at 8³
       (6.4M) and 0.8 ms at 16³ (5.6M) in Chrome, 0.7 ms on Dawn; 2.5M visible 0.6-0.9 ms.
     - **The design** (decided):
       - **Store pages by slot**: a slot list's page (a tile: its slots `[start, start + count)`) takes store pages of
         2^13 particles (not the prototype's 2^16: a tile its own pages, a 16³ yeast tile 11-20k particles, so half a
         page lost a tile and its threads idle), each its first slot (`slotBase`): an instance's identity its store
         index (its page's and its index; its slot by the page's record: picking; colors and markers by store index). A store's values one format (positions float32 or a grid's words, rotations none, float32 or
         snorm16; the first tile's, another converted), uploaded as stored; a scale a particle where a target scales by
         its radii and they are not by run (a third buffer). The buffers grown by a GPU copy (not uploaded again); a
         page let go free-listed (fixed size: no compaction). A page's record: its grid's origin and step, its runs'
         first and count (or its members a particle), its count, its slot base, its box (its positions' extent, for the
         object's bounds; no member boxes).
       - **Runs**: a page's runs a block of a runs pool (first, count, member, scale; a member a pool member as now:
         a target's geometry and levels); a page of more runs than an eighth of its particles (a flat list's targets
         alternating) a member a particle instead (u32).
       - **The culler of a source** (a module of its own the merged culler uses for sources; the lists, the raster
         and the passes as they are): a particle a thread, three passes: count (the frustum, the level by distance, the
         occlusion pyramid, the raster's split; a member and level's count, a workgroup's counts first), scan (the
         regions and the indirect records' counts, members x levels, on the GPU), write (a visible instance's pose into
         its region). The CPU lays out nothing of the instances: a tile's change its pages, runs and records alone.
       - **A composed pose**: 32 B, `t` and the scale f32, the quaternion snorm16 (a pack's as stored; float32
         rotations ~1e-4 rad), the slot u32; the culled programs read instance `i` of a member and level's draw at
         `poses[base + i]` (objects not of a source as now).
       - **Gone for sources**: the tile store's instance arrays and `TileStoreCulled`, `CulledStore`,
         `expandBlocks`, the identity region, the culler's `syncLayout` (cells, records, `refreshSegments`), most
         of `MergedSourceCounts`. Kept: the member pool (`fillMergeRecord`, `replaceGeometry`), the records'
         compaction (empty draws cost), the lists, the raster, the passes.
       - **Steps**: 1. the page store and runs (CPU; tested against the tile store's poses); 2. the source culler's
         kernels (count, scan, write; levels, the frustum); 3. the draws by pose region (a slot view's image the flat
         view's); 4. colors, markers and picking by slot; 5. the occlusion, the raster, the ghost, clipping; 6. the
         replaced paths removed; 7. the CPU-culled views a dispatch with their camera.
       - **Step 1 done**: graphics `target/page-store.ts` (`ParticlePageStore`: pages of `StorePageSize`, their records
         `PageRecordWords`, the runs pool by `RangeAllocator` (now with `grow`), a page of many runs its members a
         particle, a run of a target scaled by radii given a particle `ScaleByParticle`, the first list's format and
         others converted, `compose` the CPU's poses by slot); `target/member-pool.ts` (`TargetMemberPool`: the
         members, their counts and centers, a freed one kept for its target; the tile store's now). Tests:
         `page-store.test.ts` (the tile store's members, counts and poses by slot, with radii a particle and by run,
         a target scaled, one not accepted; a list over two pages; members a particle; a let go's pages, runs and
         members freed and reused; growth; a grid store given a float list), `webgpu-buffer-pages.test.ts` (`grow`).
       - **Step 2 done**: graphics `gl/compute/paged-cull.ts` (`createPagedCuller`: `sync` a store's changes, `cull` a
         view; `PagedCullSource`, the store one) and `gl/webgpu/shader/paged-cull.wgsl.ts` (`paged-cull-count`,
         `paged-cull-write`). Its buffers within a stage's eight: the tables one buffer (the records, the runs, the
         members' centers and reaches, the active pages; written whole where a section's capacity changed, else by what
         changed), the positions and rotations as stored, the extras (a page's member and scale a particle: the store's
         extras pool, the record its first, no longer \`PageValues\`'; now the tables' last section, grown in place: a binding freed for the nested lists' cull, particle-nesting.md step 5), the counts (the write's cursors after the scan),
         the bases (the prefix sum's), the poses (32 B: \`t\`, the scale, the snorm16 quaternion, the slot, the member).
         A level the first whose distance range holds the instance's center (none: culled); a workgroup's counts by its
         runs' index times the levels and the level, one global atomic a key. The values' buffers grown by a copy,
         submitted before the pages written then (a write runs ahead of the commands not submitted: the copy over it).
         Tests: \`webgpu-paged-cull.test.ts\` (the regions' counts and instances as the CPU classifies the composed
         poses, their poses within 1e-3 and the quaternion within 1e-4; two pages, a float list converted, members a
         particle, a target scaled by radii, one left out, a let go; a second sync growing the buffers; a level always
         0, a key's region without its level, a particle between runs drawn, and no submit before the writes each fail
         it).
       - **Step 3 done**: graphics `gl/renderable/merged-paged-culler.ts` (`createPagedMergedCuller`, a `MergedCuller` the
         merged object takes where its source is a page store's, `MergedInstanceSource.paged`: the store, its changes, its
         members' reaches): its sync the store's changes, the segments' data (`sCulledData` of no instances) and the
         draw records' words a member changed (`Merged.memberChanges`); its dispatch the paged cull and
         `paged-cull-args` (a record's instance count and first instance from its region's count and base), none where
         the view and the data are the same; its lists a list a level, their records predicted by the members' spheres
         (the boxes of the pages their runs are in, padded by their reach: a draw a level and member costs drawn or not)
         and compacted by a read of the records where the view is still; views not culled on the GPU draw none (the
         image pass among them: step 7). The culled programs' `dCulledPaged` (`culled.ts`: a pose at `instance_index`,
         its store index `aInstance`, its member's segment); precompiling leaves it out (an object of params is not a
         page store's). `target/paged.ts` (`createTargetPagedObject`, the slot list's object of a page store:
         `TargetSlotsConfig.paged` in the target representation, off by default till step 6): its pool as the tile
         store's object's, a member's bounds none of its own (the last's the object's: the pages' boxes and the largest
         reach), the layers by store index. Measured (the yeast, Dawn): a frame's CPU 7-8 ms at 8³ and 7 ms at 16³ (the
         tile store's 11-12 and 10-11), a still frame 1.8-3.6 ms, an add of known targets 5.3-7.2 ms and 4.6-6.1 ms
         (7.4-9.7 and 5.4-6.1), a let go 3.6-5.4 ms and 4.1-5 ms; the GPU's frame as before (10-21 ms). Without the
         prediction a frame drew a record a level and member (12-24k draws, 21 ms a still frame). Tests:
         `webgpu-particle-paged.test.ts` (a slot view's frames the flat view's, no pixel differing on Dawn, tiles let go
         and added again, the object kept; a pick's store index its slot; a first instance always 0 and the colors by
         member each fail it), `webgpu-paged-cull.test.ts` (a pose's store index).
       - **Step 4 done**: the colors a block's list's into its store pages' texels (step 3), the markers by store index
         (`mark`: a slot's page and index), a store page written again its markers cleared (another tile's before),
         picking a store index's slot by its page's record. Tests: `webgpu-particle-paged.test.ts` (a pick through the
         canvas a particle of the list; a highlight's markers, cleared where its tile is let go and added again: no
         clearing fails it). A highlight is not drawn in a headless frame (the tile store's object's alike): its
         markers checked, not its image.
       - **Step 5, the occlusion done**: the paged kernels test an
         instance against the pyramid (`hiz_occlusion_wgsl`, `HiZSchema`): the first phase the last frame's, the
         occluded marked (`occludedBits`, a bit a store index); the second (`cullSecond`, `uPhase` 1) the marked alone
         against the frame's pyramid so far, its counts and bases its own, its poses in the same buffer after the
         first's room (`uPoseOffset`; its records' first instances by `uFirstOffset`: a list's uniforms set values, not
         a storage buffer of its own). The merged culler's second phase (`dispatchSecond`, `getSecondLists`), its lists
         drawn with the first's after it (a still view's next frames too); their records predicted and compacted
         alike (without: a record a level and member again, a CPU frame of 27-30 ms). A view culled without a pyramid
         culled again with one (an app's first frames). Measured (the yeast, Dawn): the GPU's frame 5-13 ms (without
         the occlusion 10-21, the tile store's 9-17), the CPU's 10.5-12.5 ms (the two phases), a still frame 1.8-4 ms,
         an add of known targets 4.8-7.1 ms, a let go 3.7-5.2 ms. Tests: `webgpu-paged-cull.test.ts` (a pyramid of a
         texel: the near one occluding the spheres in view, the second phase against a far one finding exactly those;
         no marking fails it), `webgpu-particle-paged.test.ts` (a jump's frames with two phases the flat view's).
       - **Step 5, the ghost and clipping done**: the clip objects in a section of the
         tables after the active pages (`uClipBase`, `Clip.MaxObjectCount` of them: a stage's storage buffers all
         taken; `clip_storage_wgsl` reads a row by an expression), set by `setClip` and kept through the tables written
         whole. An instance they clip whole is not drawn (`clipSphereClass` 2), with instance clipping one clipped by
         its center; the ghost's (`uGhost`) an instance in view they reach (with instance clipping its center
         clipped), into a region a member after the levels' and the raster's (`uGhostBase`), in the first phase alone,
         not by the pyramid: the first phase's pose room twice the instances with it. The merged culler's ghost records
         the last block at the coarsest level's draws (laid out from the first dispatch with the ghost), its list
         the coarsest level's uniforms at every distance, predicted (a member in view) and compacted alike; a change of
         the clip objects culls again. Tests: `webgpu-paged-cull.test.ts` (each instance in the levels' regions and the
         ghost's as `Clip.classifySphere` decides, by pixel and by instance; no clip objects in the kernel fails it),
         `webgpu-particle-paged.test.ts` (the frames with the ghost and without the flat view's, by pixel and by
         instance; no ghost lists fail it).
       - **Step 5 done, the raster's split**: the poses a page store's `CulledData` (`CulledPagedFlag` in the header's
         last word, no instance arrays: the segments from 1, pose `i` at the slots' first plus `2 i`, its
         `aInstance` and segment its last words, the quaternion snorm16), the culler's `setPosesFront` (written again
         where the buffer grows), so the impostors (`dCulledPaged`, the shared readers now), the raster and the resolve
         read the poses where they are: no per-instance arrays (12 B an instance more, 260 MB at 6M with the poses).
         The kernels' `rasterized` as the merged culler's (a member's raster row in the tables: its padding at each
         level's scale, times the instance's scale); a raster instance into its level's raster region (level-major,
         `uRasterBase`; a workgroup key's low bit), the records' blocks the levels', the raster's (from the first
         dispatch with the split), the ghost's. The raster's regions kernel finds a page store's slot's level by the
         records (the first region of each), the slot the pose's index; its threads the poses' room (a slot's base not
         on the CPU); the raster records of members past the regions first past all. The lists with the raster's blocks
         where the raster did not draw (`withRaster`); the raster's rows predicted only with the split on. A culler's
         first sync uploads every page (a second culler of a store, a scene of its own, drew none: the changes taken by
         the first). The raster rows read again where a member changed (all of them 2.7 ms an add). Measured (the yeast at 8³, Dawn): an add of known targets 5.4-8.1 ms (before 5.5-7.7), its frame's GPU 12-16 ms (10-18); the raster at `auto` no faster in the test's view (11-17 ms, off 13-18). Tests: `webgpu-particle-paged.test.ts` (every particle into the raster's regions at a large
         radius, the frame the flat view's and the impostors'; the regions kernel reading the indices fails it; a later
         scene's culler drawing the store).
     - **Step 7 done, the views of their own** (before step 6: those views drew nothing): `Renderable.prepareCullView`
       (`Renderer.prepareCullView`, no pass open) culls a page store's object for a view of its own, no pyramid, the
       ghost kept: the thickness pass (illumination, a frame of its own) and a ray's pick (`RayHelper`). Into the same
       buffers (decided: no second set, 200 MB at 6M; a cull more, ~2 ms of GPU at 6M): its `cull` draws them where of
       the same camera, whatever the frame (the thickness pass's without the occlusion a superset, the frame's passes
       after it drawing that); the next `prepareCull` culls the frame's view again (the dispatch's key). A second
       camera of a frame (an XR eye) culled by `prepareCull` too (a page store's, by its camera plane: a jittered sample
       the same). Tests: `webgpu-particle-paged.test.ts` (a frame and camera of its own nothing without, its image
       with; the main view's again; the same camera's view after the frame's; an XR eye's: no view dispatch, no eye's
       each fail it).
     - **Step 6 done, the replaced paths removed**: the page store's object the slot list's (`TargetSlotsConfig` gone;
       WebGL draws a slot list flat, `ParticleSlots.flatten`, as before). Gone: `target/slots.ts` and
       `target/tile-store.ts` (`ParticleTileStore`, `TileStoreCulled`; the placement now `StorePlacement`, the
       object's interface `paged.ts`'s), `CulledStore`, `CulledPoseExpander.expandBlocks`, the merged culler's source
       paths (its data by capacity, the source's layout key, the counts' copies, the gaps' segment skipped:
       `uSkipSegment`), `MergedInstanceSource.culled` and `getMember` (its `paged` required), the CPU culling's
       source grids. Tests: `page-store.test.ts` against the poses computed (`t = p - s R c`), not the tile store's;
       `webgpu-particle-slots.test.ts` the paged object's (its CPU culling's test gone with it).
     - **Measured in Chrome** (the Viewer, WebGPU, 1600x900, the yeast in 8³ tiles, a slot view by the camera
       following it, its budget 150k particles a step; Playwright with a CPU profile): moving (a zoom to a quarter of
       the distance, a half orbit, 6 s each, 45-105 view updates) 60 fps (p50 17.5 ms, p99 20-30 ms, no long task but
       one of 60 ms), the CPU 24% busy (a frame 2.5 ms, an update about 7 ms: `setMemberSpheres` 3.3 ms of it, a pass
       over the store's pages); still 17.4 ms; all 512 tiles (31M) drawn 37-42 ms a frame moving (the GPU's). The
       first fill (290 tiles, 25.3M particles, from a whole-cell view): 11.2 s in 203 steps of 57 ms of CPU each
       (frames p50 30 ms, p90 72 ms, steps of 60-210 ms as long tasks): a tile's particles gathered from the list
       (`ParticleTiles`' `subset`: `gatherValues`, `gatherRotations`, `compact`) 22%, the object's update and
       colors 11%, the culler's uploads 11%, the layers' 7%, the frames between 10%; the tile set made first 1.9 s
       (`fromFrames`). Found there: a grown buffer's size not a multiple of 4 (1.5 times: WebGPU's error, the store
       not drawn), the extras written whole at each growth of the store (200 MB at 25M, a 2 s step): both fixed.
     - **The poses within the device's buffer limit**: the poses' room was the store's instances twice (the two phases):
       31M of the yeast at 16³ tiles 2 GB, past WebGPU's 2 GiB (`Device.createBuffer` refused, nothing drawn). Now
       within `ComputeContext.maxStorageBufferSize` (three quarters the first phase's where both do not fit, the
       buffers' growth capped too); a pose past its phase's room not written (`uPoseRoom`), a record's count cut at it
       (`paged-cull-args`): a view of more instances than fit draws those that do. Tests: `webgpu-paged-cull.test.ts`
       (a limit of 4000 poses: the buffer within it, the room's poses all written). The tiles' defaults for testing: 8
       cells an axis, a view of the tiles seen by the camera following it, a slot list.
     - **Next: tiles as views** ([particle-tiles.md](particle-tiles.md) Design 10, Steps 8): the store's pages from a
       tile's view of the list's own blocks (tile-ordered, a directory of ranges, boxes and runs), colors by member,
       steps by time: the first fill's per-particle CPU work (a tile's gather, the page boxes, the color layers) gone.
     - Open: the regions' memory (32 B a visible instance, the count pass bounding it); the raster's threads over all
       the poses' room (a slot's base read on the GPU would bound them); WebGL (no compute: 6e's
       fallback).
7. **Editable lists**: allocate, free, compact (its remap); a procedural example (a state transform adding particles
   as it runs). Tests: the slots kept across edits, a compaction's remap.

## Open questions

- **Editable lists' pages**: the size they allocate (decided: any size; an editable list's choice), a particle freed
  inside a page (a gap in it: its live ranges, or the page split).
- **Instance gaps or compaction**: a removed range left as a gap the culling skips, or the instances after it moved
  (their indices in the merged object's ranges); when to compact.
- **The hierarchy's colors across loads** (decided, done in 6d5): by the catalog's compartments and their entities
  (`ParticleCompartmentInfo.entities`), stable across views; its structure still collected again a version where the
  catalog has none (a slot list not of a tile set).
- **WebGL**: its CPU culling over gaps (the instance grid's cells by live range).
