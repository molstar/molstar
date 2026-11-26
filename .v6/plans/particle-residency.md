# Plan: Particle residency (pages made resident by what the cull draws)

## Status

- **Steps 1 to 3 done** (a template's leaf pages by its copies' levels, drawn coarse till written, its copies' poses
  from its values, the values by page, the leaf pages by the cull's wants; templates tile-ordered), 4 open, step 3's
  load cost open. One residency for every list a page store draws: a nested list's templates (resident whole before step
  1), a flat list (whole) and a tiled list (by the CPU's camera pick, [particle-tiles.md](particle-tiles.md)). Builds on
  [particle-slots.md](particle-slots.md) (the page store and its cull), [particle-nesting.md](particle-nesting.md)
  (templates, their coarse levels, the copies' level by distance) and the tiles' loader.

## Question

What a frame draws of a store against what it holds. A template is written whole at the first sync: every page of its
particles and of its coarse levels, their positions and rotations uploaded (the yeast nested 2x2: 4,391 pages, 576 MB,
about 0.8 s of its load) and given pose slots (1.16 GB). A view needs far less: the far copies draw a coarse level
alone, the near ones what the frustum and the Hi-Z pyramid leave. Residency by what the cull draws takes the load and
the memory from the template's size to the view's.

## Measured

The yeast nested 2x2 (a template of 31.5 M particles, levels of 3.94 M, 494k and 65k; 4 copies), 720p, 45° field of
view; by the CPU (the cull's level rule, `copyLevel`: a copy's template spacing at its center's distance in pixels, a
level each halving below 4; a page's bounding sphere in each copy's pose against the frustum: no occlusion), and on the
GPU (headless WebGPU, 1280x720, Hi-Z on, both phases: the pages of the poses the regions hold). The pages a view needs
of 4,391; spatial: the same drawn instances' pages were the leaves stored in Morton order across targets:

| View                              | Copies' level | Drawn  | Frustum, by target | After Hi-Z, by target | After Hi-Z, spatial |
| --------------------------------- | ------------- | ------ | ------------------ | --------------------- | ------------------- |
| whole                             | all 2         | 0.22 M | 61 (1.4%)          | 56 (1.3%, 7 MB)       | (no leaves)         |
| half the distance                 | all 1         | 0.78 M | 481 (11%)          | 422 (9.6%)            | (no leaves)         |
| a quarter (four cells, all whole) | all 0         | 1.52 M | 3,828 (87%)        | 3,118 (71%)           | 2,120 leaves (55%)  |
| a copy filling the view           | 3 at 0, 1 at 1 | 0.78 M | 3,564 (81%)        | 2,679 (61%)           | 862 leaves (22%)    |
| at a copy's surface               | 0             | 0.39 M | 3,745 (85%)        | 2,673 (61%)           | 1,095 leaves (28%)  |
| inside a copy, at its center      | 0             | 1.40 M | 2,570 (58%)        | 2,030 (46%)           | 506 leaves (13%)    |

Far, the levels carry it: a load at the whole view needs 7 MB of 576. Near, the Hi-Z pyramid leaves about 1% of the
144 M instances, but by target a page is a target's particles spread over the cell: the few drawn touch most pages; by
space 2-4 times fewer (13-28% a copy near; four cells whole in view 55%: their surfaces are). The isg grids: nothing
to gain (their template and levels 67 pages, 9 MB; 85-98% needed once a copy is near).

## Design

- **Coarse resident, leaves on demand**: a template's coarse levels always resident (small: the yeast's 13% of its
  particles), its leaf pages where wanted. A page not resident is not drawn; what draws in its place is the fallback
  (below). A tile set's (a large one's levels a seventh of it) its top levels resident, the rest on demand too (step
  4d).
- **Wanted by the cull**: the cull marks a page wanted where a copy at that page's level has it in view (a bit a page,
  a buffer read back as the cursors are: one read at a time), before occlusion of a page not resident (its particles'
  depth not known: tested by its bounds against the pyramid). The CPU makes the wanted pages resident (written from
  the pack's bytes: views of a pack in memory; a remote pack's by the loader's range reads) and lets go of the pages
  unwanted for a while past a budget (bytes, least recently wanted first).
- **The fallback**: a copy's leaf pages wanted and not resident drawn by its next coarser level. By copy (the copy at
  level 1 until its wanted pages are resident: simple, a near copy coarse for the frames a load takes) or by region
  (a tile's coarse node in place of its leaf pages: needs levels aligned with the pages, a tile order's).
- **The page order: tile-ordered** (decided): spatial chunks of about a store page's particles, by target within each
  (`templateTilesOf`), not by target (a target's particles one run: an entity's markers on few pages, its runs fewer).
  Why: a near view's pages 2-4 times fewer, and the cull's group tests (a group of 256 particles a part of a tile, not a
  target across the cell) cull whole what a pyramid hides: the isg 16x16x16 near and turning (1280x1130), its frames
  39.9 ms by target, 29.2 tile-ordered (16 tiles an axis 28.8: the default kept;
  [particle-nesting.md](particle-nesting.md) 8a). A pack written before tile order culls as by target: written again. A
  template tile-ordered is a tile set: its tiles' coarse nodes the fallback by region.
- **The pack** (not released: its format free to change): what residency reads written as it is read. A template's
  coarse levels blocks of their own (no sample at load), its particles in the page order chosen, a page's positions
  and rotations the store's page (8,192 particles: uploaded as views of the file, a remote pack's page a range read),
  the pages' bounds for the cull's wants before they are resident.
- **What stays**: keys and loci (a particle's index its own, resident or not); the markers and layers by store index
  (a page's written with it); a copy's blocks of sparse overrides by its template's store indices.
- **WebGL**: no compute cull: the CPU's level rule by copy (its center's distance) and frustum by page, or the tiles'
  camera pick as today.

## Steps

1. **Leaves by the copies' levels** (no readback): a template's leaf pages resident where a copy in view is at level
   0 (the cull's level rule on the CPU, a copy's sphere against the frustum), else let go after a while; the coarse
   levels always. The far views' gain alone (the whole view 7 MB, not 576). Measured: the load, a zoom's first frames,
   the frames a page's upload takes. Done (`paged-cull.ts`): a leaf block (a template with coarse levels) not written
   nor given pose slots, its entries none, till a cull wants it (`getLeafWanted`: a copy's level by `copyLevel`'s rule
   from 0.75 of its threshold on, its padded sphere against the frustum), written in that cull (no copy drawn without
   its pages), let go 10 s after the last cull that wanted it (its pose slots: none since step 2, the let-go with them);
   the frame's
   slots by its instances resident or not (a copy drawn coarse meanwhile). The yeast nested 2x2 (WebGPU, Chrome, 720p,
   whole): its add and commit 0.51 to 0.28-0.31 s, GPU storage 2,208 to 1,245 MB (the leaves' pose slots; their
   positions' and rotations' buffers still by the store's capacity, not written), its frames as before (still 1.13
   and 0.92 ms, rotating 1.68 and 1.55, hover 1.78 and 1.73). The first zoom in pays the upload: headless, 1280x720,
   the frames after a jump from the whole view to a quarter of its distance 2.2 s before (the cull's rooms grown), 3.1
   s now; the load's first frames 1.7 to 1.3 s; later zooms as before. Tests: `webgpu-paged-cull.test.ts` (far: the
   coarse level drawn, the leaf pages unwritten; near: made resident in that cull, drawn; far and near again; the
   frame's slots not counted with the leaves fails it).
2. **The fallback by copy, the upload spread, no template poses** (done; the cull's wants by page in 2c). Profiled
   first: the first zoom's hitch not the upload but the poses' buffer grown by the leaves' slots (`createBuffer` 1.7 s,
   the submit copying the poses 1.4 s: a template particle's composed pose kept, 32 B, the yeast's 1 GB). Done
   (`paged-cull.ts`, `paged-cull.wgsl.ts`): a copy's least level in the copies' table (`copyLeastLevel`: 1 where its
   template's leaf block is not resident, `copyLevel` not finer: drawn coarse, no holes); a wanted leaf block written
   32 MB a cull (`PagedLeafUploadBytes`), the view culled again till all are (`frameVersion`), resident then; a
   template's pages no pose slots: a copy's visible instance composed by the cull from the values as stored
   (`composeCopy`: `t_x = p - s R_x c` from the page's rotation and its member's center, as the compose kernel's), the
   compose's of a template's page its group boxes alone; a page's positions and rotations one block in one buffer
   (`pageBase`, `rotationOf`: the cull's eight bindings kept); no let-go (nothing freed: the values' buffer by the
   store's capacity). The yeast nested 2x2: the frames still, rotating and hover as before (still 0.90-0.92 ms against
   0.93-1.03, rotating 1.60 and 1.68, hover 1.70 and 1.75; the isg 16x16x16 0.90 and 0.90, rotating 3.29 and 3.21),
   the near view's 2 ms as before; headless, the jump to a quarter: 15 frames of 50 ms (the upload, coarse), then
   0.85-1.05 s and 0.45 s (2.3-2.5 s in all, 2.9 before; at most 1.13-1.21 s before). What is left of it the first
   frame of many leaf instances' buffers grown (the sphere raster's scratch 0.34 s, the cull's indices and occlusion
   bits 0.17, the frame's slots 0.04, the submit clearing them 0.64): the same before step 1 (a first zoom's), not a
   template's. Tests: `webgpu-paged-cull.test.ts` (a block of three pages a page a cull, its copies coarse meanwhile,
   another cull asked for; the copies' least level 0 fails it), the copies' poses against the CPU's as before.
   - **2a. A first near view's growth** (done). Measured first (headless, logging the buffers made and the cull's
     counts): the first cull drawing the leaves wanted 27.6-38.9 M frame slots (1.5 M from the third on) and the
     buffers grew to it (the frame's slots 1 GB, the indices 609 MiB, the sphere raster's chunks 1.6-1.7 GB, cleared
     at the submit: a second): its selected entries past the occlusion bits' room (untested), its pyramid of the
     copies drawn coarse (a sample's gaps: little hidden, the raster off the same); the indices' and the raster's
     sizes by `PagedRegionMinRoom` (4,096 a region, the yeast's 5,639 members at four levels and the raster's: 92 M
     slots each, whatever a view drew). Done: the occlusion bits grown ahead by a resident block's entries (its copies
     drawing it times its pages); a cull's counts growing the rooms neither where its selected entries passed the
     bits' room nor in the cull making a block resident (culled again instead); the floor 512. The jump to a quarter:
     15 frames of 50 ms (the upload), then at most 161 ms (1.37 s in all; 2.3-2.5 s and 0.85-1.05 s before), the
     frame's slots 201 MiB, the raster's chunks 494-589 MiB; the frames the same (the yeast tiled 1.76-1.77 ms against
     1.80-1.91, physical 1.43-1.45 and 1.42-1.43, nested 2x2 still 0.89-0.90 and 0.95, rotating 1.55-1.64 and
     1.57-1.61, hover 1.77-1.82 and 1.73-1.80; the isg 4x4x4 rotating 1.71 and 1.90, 16x16x16 3.32 and 3.31). Tests:
     `webgpu-paged-cull.test.ts` (the cull making a block resident past the frame's room not growing it, the next
     growing it: the rule off fails it; the regions' first rooms by the floor).
   - **2b. The values by page** (done; `paged-cull.ts`, `pageBase`): a written page's positions and rotations in a
     value block of its own (the tables' section from `uValueBlocksBase`, a word a page; blocks let go reused, the
     buffer grown once a sync or a leaf block's upload for the blocks it writes), not at its place in the store's
     capacity. The yeast nested 2x2 whole: GPU storage after the runs 1,728 to 773 MB (its leaves' 576 MB and the
     buffer's slack not made), its frames and load as before (still 1.16 and 1.24 ms, rotating 1.70 and 1.70); the
     yeast tiled (its pages all written) as before (4,654.5 MB both, 1.82 and 1.83 ms). Tests:
     `webgpu-paged-cull.test.ts` (the blocks the coarse level's alone till the leaves are wanted, then both).
   - **2c. Leaf pages by the frustum** (done; `paged-cull.ts`, `updateResidency`): per page, not per block: a leaf
     page written where a near copy (`getNearCopies`) reaches it (`addReachedPages`: its record's box placed by the
     copy, its sphere grown by half its radius, `PagedLeafPagePadding`, against the frustum: on the CPU, no readback),
     within the budget a cull (all the wanted pages' blocks made room for at once), let go 10 s after the last cull
     that wanted it (its value block, `PagedLeafReleaseMs`); a near copy with a page it reaches not written held at
     its first coarse level (its least level by copy, `leastLevels`), the others as the cull's rule draws them (a copy
     not near by this view's rule at a coarse level there: held by it, `particles/nested/far` drew a copy coarse the
     reference drew whole). The yeast nested 2x2 (headless): the jump to a quarter at most 93-100 ms a frame (1.3-1.4
     s in all; 161 ms and 1.37 s before), the rotating frames as before (the isg 16x16x16 3.35 and 3.37 ms, 4x4x4
     1.58 and 1.60, the yeast 2x2 1.62 and 1.65). Its value blocks by view: whole and half 551 (72 MB, the levels'), a
     quarter, a copy filling the view, at its surface, inside it 4,389-4,392 (575-576 MB): every leaf page reached, as
     the CPU measurement said (a page a target's particles spread over the cell: its box the cell's). So near the
     saving waits for a page order by space (step 3: by space without padding the frustum's 17% inside a copy, 55% a
     copy filling the view, 72% at its surface). Tests: `webgpu-paged-cull.test.ts` (a template's three pages apart:
     the one in view written alone, drawn; another as the view moves; the first let go past the release time).
   - **2d. Wanted by the cull** (done; `paged-cull-wants`, `readWants`): the wants on the GPU, not the CPU's guess
     (2c's: a copy the GPU drew at level 0 not near by the CPU's margin drew no leaves, its pages not written: the
     yeast's a quarter 187 MB, wrong). The leaf pages' entries active again; a thread an entry of a copy at level 0
     whose page is not written (its value block none): its record's box placed by the copy, padded, in the frustum its
     page wanted and its copy's group (`copyGroup`: its top particle's copies across the template's blocks) lacking it,
     split by the select in the same cull (no readback's wait) where the copy is at level 0: both levels selected, a
     particle of its written leaves drawn where the template's spacing at its own depth is at least the coarse pixels
     (`isNear`), the first coarse level's sample beyond; else held at its first coarse level. Not split while wanted
     pages find no room (`uSplitNear`: they may never come). Grown by half its radius wanted alone (written before it
     is in view); the select draws no page not written. Why split, not held coarse whole (as first): a sample's
     particle is drawn larger about its center, so a loaded copy a few pages made coarse drew its near particles
     enlarged; a molecule focused and turned about (its slab's other pages let go after 10 s, back in view as it
     turned) jumped about 3 times larger for the frames its pages took, as if zoomed into. Near, a missing page's
     particles none till written (a hole, a few frames; a view come near an unloaded copy too); beyond, the sample as
     before. The cost only while a copy lacks pages: its leaves culled as when loaded, the sample's eighth more. The wants read back
     (one read at a time), written by the next culls within the budget, let go 10 s after the last want (a cull to
     run: a still view's skipped). Tried and not kept: the wants by the pyramid too: a page's sphere hid little (the
     yeast 2x2, each view culled a frame, tile-ordered: a quarter 539 MB either way, inside a copy 192 and 193, at its
     surface 420 and 515), and a copy held coarse drew a pyramid of a sample's gaps (more wanted, written, let go:
     a feedback). By the frustum, tile-ordered against by target: a quarter 539 and 542 MB, a copy filling the view
     344 and 476, at its surface 515 and 541, inside it 193 and 421; whole and half 72 MB. Tests:
     `webgpu-paged-cull.test.ts` (the copies split while their pages are wanted, written by the cull after the read;
     a page a cull within the budget; the pages apart: the one in view written, another as the view moves, the first
     let go; a copy lacking a page in view: its written leaves near, the sample beyond its near distance alone).
3. **Templates tile-ordered** (done; the load's cost and the fallback by region open). Done: the writer's tile order
   for a template (`templateTiles`, `templateTilesOf`: a tile about a store page's particles, 16 an axis the yeast
   cell's; a tile's by target, each target's in its sample order; its density grid kept: without it computed where
   read, 5.3 s); the levels' pick by runs across a target's runs (`pickByRuns`), the runs' values gathered by an index
   of a run a 64 particles (`getRunBlocks`: 0.75 to 0.38 s); the CPU's wants padded by a block's largest scale, not
   the store's (its coarse levels' up to tens of times larger: every page reached; a member a particle's by its
   extras'). The yeast nested 2x2 (headless, each view alone): its leaf values written, by target then tile-ordered, a
   quarter 422 and 187 MB, at a copy's surface 572 and 541, inside a copy 448 and 214 (575-576 MB each before the
   padding's fix); whole and half 72 MB both. The benchmark: representations 0.72-0.82 to 1.42-1.46 s, add and commit
   0.21-0.23 to 0.31-0.33, GPU storage whole 773 to 840 MB, frames still 0.89 and 0.98 ms, rotating 1.56 and 1.58,
   hover 1.69 and 1.75; the pack 714 to 732 MB (its runs, 2.3 M: thousands of targets a tile). Then the levels' cost
   of its runs (profiled against the pack by target): a sample's targets by its runs, its entities, compartments and
   radii by run where a target's particles have one each (`getPerTarget`: the runs' ends walked together, a table by
   target id; else gathered), the pick by runs with a table by target id: representations 1.08-1.11 s tile-ordered,
   0.61-0.63 by target (1.42-1.46 and 0.72-0.82 before), frames as before. Open: the rest (the store's add walking
   the runs and a page's members a particle where its runs are many, its extras: the GPU's tables 67 MB more; the
   positions' and rotations' gathers), the entity blocks of sparse overrides measured, the fallback by region. Tests:
   `pack.test.ts` (a template tile-ordered: its instances the written list's, its targets' runs apart, its pick by runs
   a pass's, its picks over its box).
   The load's cost of tile order (the yeast nested 2x2, the benchmark's setup, two rounds): representations 1,117-1,153
   ms tile-ordered, 633-653 by target; add and commit 343-449 and 201-233; first frame 449-492 and 415-420; GPU storage
   211 and 144 MB; frames alike (1 ms). By a CPU profile, the difference: the template's coarse levels sampled at load
   (`getTemplateLevels`, 408 ms: the pick by runs 131, the positions' and rotations' gathers 134 (by index,
   tile-ordered: scattered), the per-target values by run 93, the sample order 34); the store's add walking the 2.3 M
   runs (`forEachRun` and its callback, about 200 ms: `accepts` 21, the pool's `acquire` 15 a run); the merged culler's
   sync (about 180 ms more: the members' spheres by page 94, `getPageMembers` and `growBoxes`; the tables and layers
   larger). Then: the add by the runs' arrays (4a), the template's levels stored (4b); the members' spheres a template's
   placed box once a member, not once a page (its pages share it), a page's members by its runs inline: add and commit
   about 355 to 324 ms. A template an older writer wrote warned of where read (not tile-ordered, its levels not stored:
   `warnOlderTemplate`).
4. **One pack form: a template a tile set** (WebGPU's residency; a flat list as it is: written whole, its pages drawn as
   they are; WebGL the CPU's pick and the flattened list). Nested and tiled are two axes, not two kinds: a cell placed
   by copies, and a cell's particles in tiles with coarse nodes. The form: a top list of copies and its templates, each
   a tile set (its tiles, their directory, coarse nodes by tile); a flat tiled pack the case of one template of one
   identity copy. Today a tile view's tiles go into the same page store and cull as a template's pages but as pages
   drawn as they are (`paged.ts`: a tile a block, written whole, persistent pose slots): none of the residency
   (`isLeafBlock`), the near-first order (bucket 0), the select's tests or the composed poses; a template's levels are
   samples of it whole (a copy one level), not by tile; a template is in memory whole (the yeast's 714 MB), a tile set
   read by range. What it buys: GPU memory (the yeast tiled or flat 4,654 MB, its pages all written and 32 B poses; the
   nested 2x2, four of them, 775-842 MB), a tile's cost by what a view needs, the near jump's fallback by region, one
   path on WebGPU, a single cell streamed as a template. What it costs: the CPU's pick (occlusion by optical depth, its
   `maxParticles` budget, a region switched all or none) made again on the GPU; a remote pack's latency (wants read
   back, then fetched); the poses composed every frame; CPU consumers of a template partly loaded; the tile views'
   trajectory frames, editing by slots and targets on demand moved or given up; packs written again.
   a. **The store's add by the runs' arrays** (done): a target's member, acceptance and scaling once a target (a map),
      the runs' arrays walked, a run's scale by the same runs where its scales are; not a callback a run. The yeast
      nested 2x2's representations 1,117-1,153 to 1,064-1,114 ms tile-ordered (by target as before, 640-719): the
      profile's 200 ms mostly the runs split at the pages' ends (2.3 M), not the callbacks; the add 422 to 327 ms.
   b. **A template's levels in the pack** (done, whole: by tile measured, left for c). Its coarse levels written as the
      load samples them (`ParticlePackDocument.templateLevels`, `level.*`: of each target every `8^L`-th particle as
      written, a run a target, the square root's scale; levels 1 to `MaxTemplateLevels`, kept by the load's rules), read
      as views (`Particle.StoredLevel`, `getTemplateLevels`): the same particles, runs and scales as sampled. The yeast
      nested 2x2: the levels 461 to 118 ms, representations 1,322-1,329 to 816-817 ms, frames and storage as before
      (rotating 1.67-1.71 ms, 843 MB); the pack 731.7 to 822.1 MB (the isg's 17.2 to 18.6). By tile (a tile set's
      `node.*`: a run a target a node, at least one a target a node) measured and not taken: the yeast's levels 11-42%
      larger and 724 k runs, not 5,639 (add and commit 250-342 to 349-464 ms, rotating 1.70 to 1.97 ms, 916 MB), the
      isg's 2 levels by its 4 tiles an axis, not 3 (far rotating 2.40 to 3.93 ms, 594 MB). The levels by region (and a
      template's tiles page-aligned) with the fallback by region (c), its gain against that cost. Tests: `pack.test.ts`
      (a template's levels read the samples its load takes).
   c. **The fallback by region** (its levels by region and a template's tiles page-aligned first: b's cost measured):
      `lacking` by (copy, node): a node lacking a leaf page in view held coarse whole, its tiles all at its level (not
      one tile lifted: its siblings drawn with it at both, doubles; a tile view's `getSwitchGroups` the same rule).
      Test: no instance drawn twice, no hole, while pages come in. Measure: the yeast nested 2x2 near jump (93-100 ms
      frames, 187-541 MB written), the isg 16x16x16 rotating as it is.
   d. **A template's pages on demand** (1 and 2 done): a template partly in memory, its pages read by range when the
      cull wants them. What reads a template's particles (searched): on the WebGPU nested path the store's add (its
      pages' values and boxes: `getPageValues`, `writeBox`; a member a particle's radii), its levels' gathers (`subset`:
      its attributes, its radii where not a target's each) and the colors of the attribute and index themes
      (`fillParticleColors` over it whole); else only runs, its boundary, density grid and stored levels. A focus or a
      measurement reads a loci's own particles (a copy held whole its boundary); flattening (spacefill, orientation,
      fibers, a draw not on the GPU's cull, `toFlatLoci`), WebGL's expanded copies and the writer read it all. A store
      page is a range of the template's particles, its blocks' bytes one range each: read by page, not by tile (no
      page-aligned tiles needed).
      1. **The pack opened by range**: a nested pack's top list and its templates' small blocks read whole (the targets'
         runs, the tile directory, `level.*`, the density grid, the boundary, the sites as now), their particles' blocks
         (position, rotation, attributes, a particle's radii) not (`ParticlePackStream`: a file's, a URL's, bytes'). A
         template a list of deferred particles: its counts and runs, no arrays (a consumer of all of it awaits them
         read: flattening, WebGL's copies, the writer). The store's add from the directory: a page's box the union of
         its tiles' (`tile.ranges`, `tile.boxes`: conservative), its values none till read (the cull skips them:
         `updateResidency`). Measure: the load of the yeast nested 2x2 and its heap (today the pack whole in memory, 822
         MB).
      2. **Pages read when wanted**: a page the cull wants and has no values read (the store's page range of the
         particles' blocks: batched, in flight once, `ParticlePackSource.batched`), its values filled, written within
         the budget; its values let go on the CPU least recently wanted past a byte budget (read again when wanted). The
         budget decided here: bytes resident on the GPU (a share of the device's storage) and on the CPU, least recently
         wanted first. Measure: a zoom into the yeast nested 2x2 (bytes read, culls till full detail, the coarse
         fallback meanwhile), a local file and a URL.
         Done of 1 and 2: a nested pack opened by range (`openParticleListFromPackStream`; the plugin's
         `ParticleListFromPackStream` where the GPU's cull draws it, `isParticlePackListStreamable`: a template's levels
         stored, its tile directory), its blocks read apart (`ParticlePackBlocks`) but its templates' particles
         (`Particle.DeferredParticles`: empty arrays of their kinds, read by range); a deferred template's pages boxed
         by their tiles, their values read when the cull wants them (`ParticlePageStore.requestValues`, the view culled
         again as each arrives; a template of no coarse level wanted whole), what reads all of it refusing
         (`assertRead`: a subset gathered, a sample order, flattening, the writer, WebGL's copies). The yeast nested 2x2
         opened: 304 MB read, not 784 (its sites, levels and runs; the sites the most), its buffers 371 MB, not 855, 173
         ms, not 275. A headless view of a template's corner: 20 of its 25 pages read (their boxes the union of the
         tiles each spans, padded). Then the budget: a deferred page's values given up on the CPU once written
         (`ParticlePageStore.releaseValues`: the GPU's copy), kept in a cache of 256 MB (`PageValueCacheBytes`, the
         least recently given up let go first; read again past it); the value blocks within a byte budget
         (`PagedValueBytes`: the device's storage buffer limit by default, the values' buffer one binding), past it the
         leaf pages the last read of the wants did not want let go for room (the least recently wanted first), a wanted
         page without room left unwritten (its copy coarse). The yeast nested 2x2 (headless, 1280x720, opened by range
         from the file): whole, a quarter and a copy read none (their levels'); inside a copy 77 MB (about 2 s till its
         wants settle), panned across it 106 MB in all; its values on the CPU the cache's (106 MB; all kept without
         the release: the same here, a tour past 256 MB capped); the GPU's storage not shrunk by a release (its blocks
         reused). A deferred template's attributes none (none read for a list with templates: 3).
      3. **The levels' and the themes' particles** (not needed: nothing to do): no consumer of a list with templates
         reads its templates' attributes (not kept for one: `ParticleNesting.flatten`; the attribute theme the top
         list's, none), the index theme colors by the top list's indices (a template's particles past its count the
         default color; no particle read), a particle's radii are read at open (not deferred: the levels gather them),
         the yeast nested 2x2 has no attribute blocks (its radii and entities a target's). Themes of a template's own
         particles (its attributes, its indices) a feature of their own for lists with templates, by page then.
      4. **A loci's particles** (done): a focus or a measurement of a copy held in part by its
         particles' pages read (or, not yet read, its tiles' boxes: conservative, at once), not the template read whole.
         Done: by its tiles' boxes alone (`ParticleNesting.getDeferredBoundingSphere`: the tiles holding the loci's
         particles, each grown by their largest radius; `Particle.getBoundingSphere`, `forEachLociSphere`), a deferred
         list's boundary not stored refused. The yeast nested 2x2, a copy's: its largest target (8.2 M particles) 26,536
         Å not 23,508, one of 398 particles 13,747 not 10,967, a particle 2,360-2,440 not 44-123 (a click's focus a
         tile's). Then the pages read kept in the model (`ParticlePageCache`): a deferred list's
         particles read by page (`DeferredPageSize`, the cull's page: `StorePageBits` checked equal), once while in
         flight, kept least recently used past its bytes (256 MB), a failure not kept; one a pack opened by range
         (`pageCacheBytes`), its templates' and their leaf lists' (no module state: the test files share their worker's
         modules); `DeferredParticles.getPage` (kept, at once) and `readPage`. The page store reads through it (its own
         cache gone; a page's values a view of the read where its format is the pack's; its retry kept): two stores of a
         template read a page once. A loci's sphere by page: a kept page's particles by their positions and radii, the
         rest by their tiles' boxes; no read begun. The yeast nested 2x2, a particle of a copy: 44-66 Å once its page is
         kept (the whole read's), 358-2,366 by its tiles, under 1 ms either way, nothing read. Then a particle of a
         stored level by its position there where its page is not kept (a coarse copy's click: the pick's index the
         template's; a page's loci of at most 64 looked for, `LevelLookupMax`): the yeast's levels 1 to 3 68, 76 and 44
         Å (their pages read alike), the levels' index by page 61 ms at a template's first, 0.3 ms a focus then; one in
         no level still its tile's (2,366). Not measured: the heap with two representations. Not here: reads begun by a focus or a measurement (its sphere refined then), one
         budget across packs. Tests: `particle-page-cache.test.ts`, `pack-stream.test.ts` (a page read once, a loci by
         its particles where kept, by its tiles past the cache's bytes), `page-store.test.ts` (a second store's page
         from the cache).
      5. **The comparison** (measured; e's decision open): a single yeast cell as a template of one copy
         (`yeast_cell_atoms_nested-1x1.mspk`, the grid script at 1x1, 822 MB) against its tile view
         (`yeast_cell_atoms.tiles16.levels.mspk`), in the Prism viewer (Chrome, WebGPU, 1920x1080:
         `scripts/benchmark-particle-tiles.mjs`, its path outside, an orbit, a dolly in, an orbit inside, still inside;
         two rounds alike). Not alike in what is drawn: outside the template's copy is at full detail by the level rule
         (its pages all read, 480 MiB), the tile view's CPU pick about 10 M particles (`maxParticles`). Opened (its
         bounding sphere): tiles 1.2 s and 157-161 MiB read, the template 2.0-2.2 s and 307 MiB (its sites, levels,
         runs). Outside to full detail: tiles 9.2-9.4 s after the view rests (about 470 MiB read in all), the template
         2.0-2.2 s (787 MiB in all); inside nothing more read by either. GPU buffers at the end 1.78-1.82 GB both (peaks
         2.08-2.12 and 2.33). Frames moving (p50 and the display's frames missed in 4 s): orbit outside tiles 2.3-4.4
         ms and 336-698 missed (long tasks 285-711 ms), the template 6.9-7.0 ms and 12-24 (none; its GPU 5.3 ms a frame
         against 2.9-3.0, three times the particles); the dolly in 2.0-2.5 ms and 1,152-1,416 missed (long tasks
         1.4-1.6 s) against 6.3-6.4 ms and 71-78 (0.36-0.39 s); orbit inside 2.1 ms and 1,393-1,454 (1.4-1.5 s) against
         5.0 ms and 15-18 (none). The JS heap at the end (ArrayBuffers not counted) 175 MiB against 148. Not
         measured: the CPU's ArrayBuffers (the tiles' loader 430 MiB, the template's page cache at most 256 MiB), a
         remote pack's latency, the template under a particle budget. The features only the tile view has: trajectory
         frames, editing by slots and targets, targets on demand.
         Then a budget (`maxParticles`, the representation's: `PagedCullView.maxInstances`): the copies' coarse pixels
         doubled a step of a bias, begun at the coarsest (no leaf page wanted before a read back), as many steps finer
         a read back as their growths keep within (eight a step, more where one was seen: the visible not scaling with
         the level, 11,618 to 118,529 once), coarser where over; a read back about a second where the GPU is busy.
         The yeast as a template: 10 M does not bind outside (about 0.5 M visible a step coarse, Hi-Z), full detail in
         four read backs (3.7 s) and its 480 MiB read then (the wants the frustum's, not the budget's); 100 k coarse
         throughout, nothing read, its GPU 8.8 ms a frame orbiting inside against 5.0 (the coarse level's particles
         larger: more fill). Open: a budget of bytes read (at full detail the frustum's wants read the occluded pages
         too), a default, a budget across representations.
   e. **The tile views' migration**, decided by d's comparison: a tiled pack written as one template of one copy, the
      tile views and their slot lists retired, their features (trajectory frames, editing, targets on demand) moved or
      given up; the CPU's pick replaced by the wants: their occlusion by the occluder grid (what was drawn not read: no
      feedback; near, the CPU's pick keeps 58 of 426 tiles inside, 1,262 of 2,747 whole), a budget by the wants read
      back (the pick's `maxParticles`), a node's region switched whole (c). Measure: bytes fetched and tiles loaded a
      view against the CPU's pick (542 MB, 4,632 responses filling the 16³), the time to the first and the full detail,
      a zoom tour.

## Decisions open

- **Step 4's**: the pack's form (a template a tile set: its levels by tile, its tiles page-aligned); the budget (d); the
  tile views' migration (e, by d's comparison); a remote pack's latency (the wants read back and the reads take culls:
  the region coarse meanwhile; the browser's limits on reads in flight); a page drawable once its targets' members exist
  (new members rebuild the object); the wants' occlusion (the occluder grid's threshold leaves gaps inside); WebGL on
  the same pack (its CPU pick and expanded copies).
- **The fallback**: by copy, split near and far (step 2d), or by region (step 4c, a tile order's: then coarse only
  where a node lacks pages, not past a distance).
- **Predictive loading and release**: the wants are of the view a cull sees; what a view is about to see is guessed
  only by a page's padding (half its radius). Let go 10 s after the last want, the pages a rotation brings back are
  read and written again (a focused molecule's slab: about 1,000 pages let go, 28 back in view as it turned), their
  near particles holes meanwhile (2d). Candidates: pages within the camera's orbit about its target wanted alone
  (`min(radius, distance)` about the target: what turning about it reaches, never let go while there; a wheel zoom
  without a focus keeps the whole scene's radius, so the distance bounds it); wants ahead of the camera's motion (its
  velocity); release by reachability, not age alone. Measure: the pages written again a turn, the near holes' frames,
  the value blocks held.
- **The budget**: the value blocks within the device's storage buffer limit (a lower share: the device's memory not
  known to WebGPU), an unwanted page let go 10 s after the last cull that wanted it or for room; a page's values on the
  CPU a 256 MB cache (a pack's, in the model); the values' buffer never shrunk.
