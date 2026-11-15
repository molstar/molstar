# Plan: Tiles of a particle list (a mesoscale model subdivided, loaded on demand)

## Status

- **Step 1 done**: `particle-tiles.ts` (`ParticleTileCatalog`, `ParticleTileGrid`, `ParticleTile`,
  `ParticleTileProvider`, `ParticleTileSet`, `ParticleTiles.fromList` and `fromTrajectory`: in memory, an optional delay
  as a fetch's), `Particle.subset` (a list of some particles: positions exact on their grid, values as runs where
  fewer, the targets the same objects; no fibers: a fiber's particles may be in other tiles). Tests:
  `particle-tiles.test.ts`.
- **Step 2 done**: `particle-tile-loader.ts` (`createParticleTileLoader`): `setWanted` by priority (a load no longer
  wanted cancelled: its task aborted at its next update, its result dropped), at most `concurrency` at once (4), a
  tile's targets with it and once each for all tiles and frames; ready kept until `release`, `releaseTargets` only of
  the unused (`getUnusedTargets`), a failed tile told (`onError`) and loaded again where wanted again
  (`ParticleTileLoadAttempts`), then failed until released; `getLoaded` (tiles, particles, targets,
  bytes), `subscribe`, `whenIdle`, `dispose`. Tests: `particle-tile-loader.test.ts` (timeouts faked, the tasks'
  scheduler immediates real).
- **Step 3 done**: `particle-tile-visibility.ts` (`ParticleTileVisibility.get`, `ParticleTileCamera`): the frustum, the
  clip objects of what draws the tiles (`ParticleTileDrawing`, a representation each: any shapes; a tile all clipped
  away left out, a density block cut exactly by planes and by the other shapes split in eight a few times, a part a
  surface still cuts hiding nothing, none split finer than a ray apart on screen, a tile's cell tested before its
  blocks; what hides what every opaque one keeps; measured on the yeast at 16³ from outside: none 15 ms, a plane
  13, a sphere 19, a cube 27; inside 17, 4, 31, 27; at 2 × 2 cells, 32³: 114, 94, 131, 121), the near plane of the
  canvas camera's clipping (its focus' radius: `Camera.getClipDistances`, nothing nearer hides), the size on screen
  (the catalog's `radius`, `ParticleTileCatalog.getTargetRadius`; the tiles' largest cached per set: the yeast's tiles
  use thousands of targets each), the occlusion by the tiles' optical depth (6). Tests: `particle-tile-visibility.test.ts`.
- **Step 4 done**: `particle-tile-opaque.ts` (`ParticleTileOpaque`), `ParticleTile.opacity`: by the content's density
  (decided). Molecules of volume fraction φ and mean chord ℓ (4V/S) let a line of length L pass with probability
  exp(-φL/ℓ), the mean over lines (a line's own count of molecules is noise: at the yeast's density a 42 Å voxel holds
  0.1 molecule, and the least line of a block was always empty); a molecule a solid sphere of its volume (fewer, larger
  objects than its atoms: errs towards passing), its volume its mass's (`ParticleTileCatalog.getTargetVolume`: element
  weights, a hydrogen a heavy atom where none, 1.21 Å³ a dalton; without elements its sites' spheres times 0.7; the
  catalog's `volume`). `forEachMolecule` (a sphere a particle), `computeOpacity` (per block of the cell, 4³ by
  default, its molecules' cross-sections S/4 over its volume: the optical depth a length), `computeOpacities` (a set's,
  a trajectory's each block's least over its frames). Tests: `particle-tile-opaque.test.ts`.
- **Measured, the yeast cell** (31.5M particles, 5,639 targets): its molecules 1.9·10¹² Å³, about 5% of the cell; a
  target's volume median 39,000 Å³ (a sphere of 21 Å); a block's density median 1.2·10⁻³ an Å (p90 1.4·10⁻³): an
  optical depth of 7 after about 6,000 Å of cytoplasm, deeper than a tile (8³: ~5,400 Å; 16³: ~2,700 Å), so tiles
  hide together. Opacities 3.9 s (8³), 4.7 s (16³). Kept, of what is in view, and the call's time:

  | tiles | whole cell           | zoomed               | inside               |
  | ----- | -------------------- | -------------------- | -------------------- |
  | 8³    | 257 of 418, 18.1M of 31.5M particles, 3.4 ms | 26 of 80, 2.7M of 7.3M, 1.4 ms | 29 of 89, 3.4M of 7.5M, 5.2 ms |
  | 16³   | 1,250 of 2,747, 12.9M of 31.5M, 15 ms | 77 of 437, 1.1M of 5.5M, 3.4 ms | 58 of 426, 0.85M of 4.7M, 11 ms |

  Blocks of 8³ a tile: about the same kept, 3-6 times the time. Superseded on the way: opaque boxes (a box every
  axis-aligned line through which crosses the content; seams between boxes merged where their union is a box): from
  filled voxels of stand-in spheres 155 of 418 tiles with a box and 10-30% of what is in view hidden; by optical depth
  no yeast tile opaque alone.
- **Step 5 done**: `particle-tile-view.ts` (`ParticleTileView.create`, `of`, `isSamePart`): the tiles of a set at a
  frame from a source of ready tiles and targets (a loader's), by target (a target one run), within a target by tile,
  within a tile in its order; keys the tiles' (their base and index); the composition per target (its tiles, their
  ranges); a frame step with the previous view of the same tiles `Particle.withPoses` of its list (`isSamePoses`, the
  composition the same object); positions on one grid kept exact, rotations as stored where of one kind (else
  float32); entities, compartments, radii, attributes concatenated where every tile has them. `ParticleValues.compact`
  (an array as runs where fewer than half; `Particle.subset`'s too). The view on the list (`ParticleTileView.of`), not
  a chain of views (`tilesFrom` not kept: each view would hold the one before). Tests: `particle-tile-view.test.ts`.
  Measured on the yeast (8³, the tiles' lists made): 56 ns a particle (100 tiles, 5.3M particles: 336 ms; 30 tiles:
  40 ms; all: 1.8 s), a frame step 25 ns; a tile's particles by target sorted once a tile (1.9 s for all, first).
  Values as runs straight from runs (no array, no `compact`) would save about a third.
- **Step 6 done**: `SO.Particle.TileSet` (a set and its loader, one a set; the loader disposed with it), the plugin's
  `state/transforms/particles/tiles.ts` (with the particles: `ParticlesRegistry`): `ParticleTileSetFromList` (a list
  or trajectory, `cells` an axis, a `delay` as a fetch would have, the loader's `concurrency`, `opacity`: the
  opacities computed) and `ParticleTilesView` (all, by id, in a box, or seen by the camera at a frame: the tiles wanted,
  their loads awaited, the failed left out; an update a frame step of the same tiles where only the frame changed);
  `state/helpers/particle-tiles.ts` (`getParticleTileCamera`: a canvas camera's snapshot and viewport as a tile
  camera, its far plane deep: the canvas' own is the loaded scene's, which tiles not loaded are not of). Tests:
  `state/transforms/particles/_test/tiles.test.ts` (the objects, a view's tiles, its representation's visuals the same
  for a tile and for all, a box, a frame step, the camera moved, the occlusion with the opacities). The opacities
  spread a molecule over the blocks whose centers it covers (a molecule larger than a block, as in small tiles). The
  camera mode following the camera: its `follow` (on unless a snapshot is given), the plugin's
  `behavior/dynamic/particle-tiles.ts` (`ParticleTilesFollowCamera`, in the Viewer's spec): as the canvas' camera
  moves (`Camera.changed`, throttled by `intervalMs`) and as tiles load (the loader's `subscribe`), a following view
  behind what the camera picks (`isParticleTilesViewBehind`) updated (the same params as a new object: applied again)
  a step at a time, a frame drawn between steps. A step (`budget`, particles): the tiles shown still picked, then the
  picked ones ready by priority up to the budget (at least one), not waiting for the rest; the first apply waits for
  one. Why: adding a view's worth of tiles at once blocked for a third of a second (the store's runs, the uploads, the
  colors). The pick kept per tile set by its camera (`pickTiles`: a camera's visibility once, the behavior's check
  and the view's step alike). The camera following stays an option (decided): a view fixed to the camera it was made
  by (a snapshot) is what tests and comparisons need. Tests: `tiles.test.ts` (a following view's tiles again as the
  camera moves, a view of a camera given and one not following as made, no update at rest; a budget a tile a step).
  Not yet: the visibility's own cost while moving (~55 ms a camera at the yeast's 16³: the occlusion's rays), the
  clip objects in it, letting go of tiles no view wants (the app's).
- **Step 7 done**: the target visual (graphics `target/common.ts`) keeps a target's instances where the list before and
  the list now are views of one set at one frame and its part is the same (`ParticleTileView.isSamePart`: the same
  tiles, the same counts): its transforms, bounds, colors, sizes and markers as they are, only its particles' indices
  in the list the new ones (picking, marking, a later color update). Its colors kept where the theme colors a particle
  by what it is: `particle-index` in a view by its key of the set's count (a particle the same color in views of other
  tiles; a list not a view by its index), the attribute theme's range the tiles' union (a list's tiles each the
  whole's range). Tests: the plugin's `tiles.test.ts` (a target's transforms the same version where a tile of another
  target is added, its picked particle the new index; built again where a tile of its own is; the index colors).
  Measured on the yeast (8³, 30 tiles, 3,345 targets, headless WebGPU): a tile added or let go 450 -> 290-390 ms. The
  rest: the targets the tile has particles of built again whole (about 95 ms: their transforms and colors; adding a
  tile's particles to a target's instances is the editable lists' part), the merged object made again (about 105 ms:
  every target one merged object, made again where a target comes or goes; an update of its members in place would
  keep it), the view (about 45 ms). Next: a tile added or let go under 10 ms (the goal), by
  [particle-slots.md](particle-slots.md): slots kept by a particle in pages, a tile its own pages, the updates by what
  changed.
- **Step 8 designed, tiles as views of the pack** (Design 10, Steps 8): a tile a range of the list's own blocks, as
  the whole list is (its grid words and snorm16 rotations uploaded as the file holds them, its targets runs); adding it
  a few `writeBuffer`s and records, no pass over its particles. Why, measured in Chrome (the Viewer, WebGPU, the
  yeast in 16³ tiles, a slot view by the camera): the whole cell filled in 16 s of CPU (2,557 tiles, 31.4M
  particles, 222 steps of 70-75 ms), where the whole list without tiles loads in under 3 s end to end. The fill's
  time: a tile gathered from the list (`Particle.subset`: positions and rotations by index, targets, entities and
  radii expanded per particle then compacted) 24%, the layers' uploads 17%, the colors on the CPU 9%, the tile set
  made (`fromFrames`: a pass binning every particle, a set of targets a tile) 13%, the store's per-page box (a pass
  over its grid words) and the culler's sync the rest; the uploads themselves (`writeBuffer`) 18%.
  - **8a done, the tile order in memory**: `ParticleTiles` (`fromFrames`) orders the list by tile once (`Particle.subset`
    of the counting sort's order: stable, so by target within a tile where the list is; values gathered as runs, not
    an array of all), a tile's list a view of its range (`Particle.range`: the words, rotations, attributes' columns
    windows, the values' runs sliced, `ParticleValues.slice`; its keys a view of the shared serial array), its bounds
    from the tiling's pass (its centers' box padded by its largest radius: `Particle.setBoundary`, no pass of its own),
    its targets by its runs. Measured (the yeast at 16³ in Chrome, as above): the fill 16.4 -> 10.8 s, a step's CPU
    74 -> 61 ms; the tile set made 2.1 -> 3.6 s (the one ordering: its values' gather 1.3 s, a run lookup a particle).
    Left of a step: the colors (`fillColors` 1.4 s, the color layer's uploads: written whole as the store grows, the
    layer by capacity) 8c's, the store's add (0.5 s), the members' counts (0.3 s). Tests: `particle-tiles.test.ts` (a
    tile's values the list's, its bounds holding its particles, the tiles views of one ordering; `ParticleValues.slice`).
  - **8c done, colors by member**: `ColorTheme.targetColor` (a target's color where every particle of it has it, by
    what it is: the entity theme's by `ParticleTargetInfo.entity`, which all its particles share; flat lists and slot
    lists the same palette indices); a slot list's paged object, where the theme gives every target's, colors its pool
    members `object` (a texel a member, read by its segment: `createObjectColor`'s values in the member's own cells),
    no color layer (none filled, none grown, none uploaded); the theme changed a member's again, the other way (a theme
    without, `particle-index`) the object made again with the layer. Measured (the yeast at 16³ in Chrome, the entity
    theme): the fill 10.8 -> about 8 s (a step's CPU 61 -> 46 ms; without the tile set's ordering in the profile about
    30 ms), frames while filling p50 18 ms. Left: the tile set's ordering (3.5 s: 8e), the store's add, the members'
    counts and creation, the markers' layer. Tests: `webgpu-particle-paged.test.ts` (the entity theme's object
    `object` colored, the flat view's frames; `particle-index` a layer again, its frames too).
  - **8b done, the store's boxes from the tiling**: `ParticlePageStore.add`'s `box` (a block's pages each that box: a
    tile's, larger than a page's own, what predicts its members' draws, never fewer; the cull kernel reads none), the
    paged object passing a tile list's bounds where known (8a's), no pass over a page's words (`writeBox` where none).
    The runs were by run already. Measured (16³, Chrome): the store's add 0.44 -> 0.35 s over the fill (1.6 ms a step:
    a tile's target runs). Left of a step (about 27 ms of CPU, the tile set's ordering apart): the object's update
    (7 ms), the merged members' segmented layers written again (3 ms, `SegmentedArray.read` copying segments), the
    frames between (3.5 ms). Tests: `page-store.test.ts` (a box given every page's, else the particles').
  - **8d done, steps by time (off by default)**: a following view's `stepMs` (0 by default: the budget alone), the
    behavior timing each step's update (the view's and its representation's) by the particles it added, a running
    mean of particles a ms a tile set, its next step's `rate × stepMs` particles at most the budget
    (`setParticleTilesStep`). Measured (16³, Chrome, `stepMs` 10): the fill's frames over 33 ms 68 -> 35 (p90 54 ->
    39 ms), but its time about 8 -> 13 s (539 steps for 222): a step has a fixed cost of about 20 ms whatever its
    particles (the culler's sync 2.1 ms, the merged members' segmented layers written again 1.9 ms, the source's
    instances laid out again 1.3 ms, the view made 1 ms and the catalog's `getTargetInfo` 0.9 ms a step, the
    members' counts 0.8 ms, the frames between), so smaller steps are only more of them, each still a long task.
    Next (8d2): a step's fixed cost, then `stepMs` on. Tests: the plugin's `tiles.test.ts` (a step's particles by a
    rate; a view's step its set's with `stepMs`, its budget without).
  - **8d2 done, a step's layers not written again**: counted (16³, the render item's data layers over a fill), the
    members' sites (`tSite`, 86M elements at the end) were written whole 15 times (1.3G elements), their changed
    ranges 19M: the object made again as the pool grew, the layers' page outgrown (moved: all written), the merged
    source laid out whole as its room ran out. Now: the pool's room for every target of the catalog
    (`targetInfo`, the object made once); a render item's layer a slot of its own (grown by half again, a layer
    grown not moving those after it) and a page outgrown copied on the GPU (`BufferPages.copy`, submitted before
    the writes after; through a staging buffer where one buffer is read and written), not uploaded again; a merged
    source's arrays extended longer where the room ran out (`SegmentedArray.extend`'s `length`: the parts where
    they were), laid out whole only where the unused draw more than the used; an extension line cut behind its
    last extensions, not at the newest (an uploader one behind wrote all every 16th). `tSite` written: 67M at
    the first upload, then 36M of changed ranges. With the object made once, no scene commit resets the camera as
    tiles come: a tile set's `boundingSphere` (its tiles' spheres', from the tiling's pass), a slot list's
    `ParticlePages.sphere`, the paged object's bounds from the first step (the whole cell framed). Measured (16³,
    Chrome, the same framing): the fill 10.3-11.2 -> 8.7-9.0 s, its CPU 10.7 -> 9.2 s (a step about 32 -> 24 ms,
    the tiling apart), long tasks 16 -> 4; its frames GPU-bound (p90 about 90 ms, the whole cell drawn as it
    fills). `stepMs` 10: frames p90 23 ms, 2 long tasks, but the fill 14.7 s (a step a frame): left off. Tests:
    `segmented-array.test.ts` (longer, the line cut), `merged-source.test.ts` (grown in place),
    `webgpu-buffer-pages.test.ts` (copies), `particle-tiles.test.ts` and `particle-slots.test.ts` (the
    sphere).
  - **8e done, the pack tile-ordered** (`ParticlePackVersion` 6, opt-in: `writeParticlePack`'s `tiles`, the CLI's
    `--tiles`, a pack written again pack to pack): the particles by tile of the written grid, the tiles in Morton order,
    a tile's by target; the targets' runs a block (`target.runs`: a tile's targets its runs' slice, no runs table of its
    own), the directory two (`tile.ranges`, `tile.boxes`); `ParticleTiles.fromList` from the directory for that grid
    or a power of two coarser (a range a tile, the list itself the ordered one, its keys its indices), else ordered
    in memory (8a). Not a fine grid: the yeast's runs (a target in a cell) at 16³ 2.29 M, 32³ 5.9 M, 64³ 11.6 M (a
    run a store record: a 16³ tile from a 32³ or 64³ grid 2.6-5× the runs). The catalog's targets' reach and volume
    stored in the pack (`reach`, `volume`: of the sites as written), its pass over 20 M sites gone (0.8 s);
    `forEachRun2` by the runs' values, no search a run. Measured (Chrome, 16³): the tile set 3.8 -> 0.4 s (Node:
    3.9 -> 0.39 s, the catalog's compartments over 2.29 M runs the most left), the fill's CPU 9.1 -> 6.1 s, the
    pack 713 -> 731 MB (its runs, 18 MB), written in 7 s. Its list drawn whole: 5.3 s against 2.1 s (a target an index
    set: its poses gathered, 0.74 s; its runs searched, 1.1 s; its bounds its runs' united,
    `Particle.getIndicesBounds`; its color constant where each of its runs' is, `ParticleValues.forEachRunOf`: else
    the members' colors of two kinds, not merged). Tests: `pack.test.ts` (tile order, the directory, the tile sets
    from it at 4, 2 and 1 as by a pass, 3 by a pass; reach and volume; version 5 read), `particle-list.test.ts`
    (`forEachRunOf`, `getIndicesBounds`).
  - **8h and 8i done, opacities in the pack, an occlusion threshold**: a tile-ordered pack's tiles' opacities written
    with it (`tile.opacity`, 4³ blocks a tile; `opacityBlocks` 0: none), read onto the directory's tiles, a coarser
    set's tile its cells' blocks side by side; `ParticleTileSetFromList`'s `opacity` computing only where the tiles
    have none. The camera view's `opticalDepth` (7 by default; its pick's key). Measured (16³, Chrome, the
    tile-ordered pack, occlusion on by default now that it has opacities): the fill 9.1 -> 4.0 s (2,636 -> 1,262
    tiles, 31.5 M -> 13.4 M particles); `opticalDepth` 3: 3.6 s (1,082 tiles, 10.7 M), the cell from outside alike.
    The pack written in 9.9 s (the opacities 2.8 s), 731.5 MB. Tests: `pack.test.ts` (the read opacities the
    computed at the written grid, a coarser one's cross-sections; none written), the plugin's `tiles.test.ts` (four
    layers, each 2 deep: all kept at 7, the two farthest left out at 1).
  - **8g measured, what a run costs on the GPU** (Chrome, timestamp queries in timing mode, the yeast's 31.5 M
    particles all shown by a slot view, occlusion of tiles off, the GPU's two-phase occlusion on; an orbit): a frame's
    GPU 78.6-81.7 ms from the 16³ pack (2.29 M runs), 60.8 ms from the 8³ (0.72 M); zooming 64.4-65.5 against
    39.5-41.0 ms. The cull about 80% of a frame, the raster 18 ms alike. Per pass the count passes read 9.8 ms against
    2.7 (3.6×, as the runs) and the second phase's write 32 against 29 (the poses, 32 B a visible particle: by the
    bandwidth). A pass's own time is not reliable here: tried, none changing the 16³ frame (78.6-79.4 ms): a run found
    from a window of the workgroup's runs in workgroup memory (one search a workgroup; slower, the 8³'s count 2.7 ->
    9.8 ms: thread 0's two searches stall the rest), 1,024 workgroup keys (no key past them, 4 levels: none), the
    count's global atomics left out (its count read 9.9 ms at 8³ too: the timestamps shift with the passes around).
    So the runs cost about 18-25 ms of GPU a frame of all particles at 16³, not yet attributed; next a pass alone (the
    count dispatched N times in a frame) or a GPU profiler, before a fix. The poses' bandwidth the larger lever
    (8j: an index a visible instance, 4 B, the pose read where drawn, not 32 B written a frame).
  - **8j started, the fill's CPU by theme, the budget**: with the pack's opacities the fill's frames are no longer the
    GPU's (the 16³ fill's at the frame cadence but for its steps' long tasks). The hierarchy theme a target's color
    where its info holds its compartment (`ParticleTargetInfo.compartment`, a pack's: stored a target) and so its
    entity's or function's (`targetColor`, as `indexHierarchy` colors its particles): the paged object colored by
    member, no per-particle layer filled. Measured (16³, Chrome, the tile-ordered pack, hierarchy): the fill's CPU
    6.2 -> 5.7 s (`fillByIndex` 272 ms gone), its long tasks 18 -> 8; moving a step 37.9 -> 30.0 ms, `writeBuffer`
    1.36 -> 0.53 s. The following view's `budget` 150k -> 1M by default: the fill 4.2 -> 2.5 s (95 -> 15 steps), the
    frames moving at the cadence (orbit p50 40.5 -> 32.5 ms, p90 49.9 -> 39.7); 5M: 1.8 s in 4 steps of 0.3-0.6 s
    each, moving worse (a fast orbit's steps large; a slow spin's smooth, no tile seen loading). Tests:
    `particle-slots.test.ts` (a target's hierarchy color its particles', by entity and function; none without its
    compartment), `pack.test.ts` (a target's compartment in its info).
  - **The least occluded first**: the visibility's tiles by their exposure (`ParticleTileVisible.exposure`: the mean
    over their pixels' rays of exp(-the optical depth the tiles in front gathered), by tenths), then larger on screen,
    then nearer: a following view's steps fill what is seen first, what shows through gaps after (the yeast's 16³ fill
    at 0.9 s: the cell's front whole, where by size a scatter of tiles). With it a following view's loader given the
    camera's pick whenever it changes (`wantTiles`, in the check whether the view is behind), not only when the view
    steps: a pick that kept the tiles shown and brought others (the first view's, made before the camera framed the
    set) had them never loaded, the view stuck. Tests: `particle-tile-visibility.test.ts` (a tile seen through a thin
    one after a smaller one in clear view; without occlusion by size), the plugin's `tiles.test.ts` (a pick keeping
    the shown tiles and bringing others: those loaded and shown).
  - **A tile-ordered pack drawn whole by its tiles**: such a list (`ParticleTiles.getDirectory`), where culled on the
    GPU, drawn as a slot list of all its tiles (`ParticleTileView.ofWhole`: a tile a page, a view of its range, a slot
    its index; kept per list) by the page store's object, its themes made again for it (`getSlotsTheme`), a picked
    slot's loci the list's (its index), its signatures none (compiled when drawn, as a slot list's): not a target an
    index set (its particles a run a tile). Measured (Chrome, the yeast whole): 5.3 -> 3.7 s, the target-ordered pack
    2.2 s; left the paged object made for all 2,753 tiles at once (the pages' members and boxes 255 ms, the store's adds
    236 ms, the pool of 5,639 members 74 ms) and the list's own theme the plugin makes (149 ms). One file for both: a
    pack written tile-ordered drawn whole within 1.7× of one by target. Tests: `pack.test.ts` (the slot list's pages
    views of the list's ranges, a slot its index; none without a directory).
  - **8k step 1 done, the pyramid primed from the tiles' densities**: a paged object of a tile view (whole or tiled)
    registers its shown tiles' opacities side by side as one grid (`ParticleTileOpaque.getDensityGrid`,
    `setOccluderGrid` by its render object's state: a merged object's values are copied), the frame hands a visible
    object's to the depth pyramid's pass (`HiZPass.setOccluder`), whose first level on WebGPU (`hi-z-occluder`:
    WebGL builds it plain) is the drawn depth's farthest or, nearer, where the texel's center ray through the grid
    reaches the threshold (7, as the tiles' picks), at the far end of its step; the camera's matrices tagged
    (`depthRange`) so the backend maps them to its reversed depth (untagged the ray ran from the far plane: all past a
    few hundred Å hidden). Measured (the yeast's 31.5 M particles all shown, 16³, Chrome): a frame's GPU orbiting
    79.6 -> 46.0 ms, zooming 58.5 -> 39.5 ms (the second phase's write 33 -> 11 ms, its raster 14 -> 2.4 ms; the
    pyramid's pass a fraction of a ms). From outside and near the surface as drawn before; inside the cell a few gaps
    where far content showed through the molecules' gaps (the threshold's trade-off, as the tiles' picks). Tests:
    `particle-tile-opaque.test.ts` (the grid of the shown tiles' blocks), the native suite. Left: the grid for
    packs by target (step 2-3), for lists without one (step 4), a threshold of its own (now 7).
  - **8k steps 2-4 done, the grid of a list drawn by its targets**: a pack by target writes its particles' density
    grid (`density`, 64 blocks along its longest side, near-cubic: the yeast's 62×58×64, 0.9 MB, 1.3 s more to write;
    `densityBlocks` 0 none; a tile-ordered pack none, its tiles' opacities are), read onto the list
    (`ParticleDensityGrid.set`, `get`); `ParticleDensityGrid.compute` the same of a list without one (its particles'
    centers' box, the catalog's volumes; `ParticleTileOpaque.computeDensity` the core `computeOpacity` shares). The
    target representation sets it on its objects (`occluder`: given, the default; computed where none, once a list;
    off) where unmoved (the grid in the list's frame) and the backend uses it (`usesOccluderGrid`: WebGPU); a tile
    view's object `off` too. Computing at load opt-in: 1.3-1.9 s on the main thread for the yeast, again for every
    derived list (a trajectory's frame). Measured (the yeast by target drawn whole, Chrome): a frame's GPU orbiting
    99.3 -> 41.0 ms (the second phase's cull 55.7 -> 12.4 ms, its raster 21.9 -> 6.0 ms), zooming 55.6 -> 29.5 ms;
    the views from outside, close and inside as before (52, 61 and 0 pixels of 926k differ). Tests: `pack.test.ts`
    (the written grid the computed one, none tile-ordered or where asked), `particle-tile-opaque.test.ts` (the blocks
    an axis, a list's grid).
  - **8k measured again** (a restart: the GPU's memory short before, the figures of 8k 2-4× high): the yeast by
    target drawn whole, a frame's GPU orbiting 13.9 -> 10.2 ms with the grid, zooming 15.0 -> 10.6 ms.
  - **8k done, the priming's threshold and the grid computed in a task**: the depth pyramid's pass's
    `occluderOpticalDepth` (`HiZProps`, 1-20, 7 by default: the canvas's, not the tiles' `opticalDepth`). Measured (the
    yeast by target drawn whole, Chrome): a frame's GPU orbiting at 3 / 7 / 12 34.7 / 40.2 / 48.4 ms, zooming 26.3 /
    30.0 / 34.7 ms; at 3 a few holes inside the cell (733 pixels of 926k), at 7 and 12 none. A list's grid computed in
    the representation's task (`ParticleDensityGrid.computeTask`: its box, then its molecules, 2^18 particles a chunk
    between the task's updates; `ParticleTileOpaque.createDensity` the core in parts): the yeast's 2.6 s more, its long
    tasks the task's update cadence (about 110 ms), not one of 1.3-1.9 s. Not a worker: none in the repository yet (an
    entry bundled per app). Tests: `particle-tile-opaque.test.ts` (the task's grid the one at once, past a chunk).
  - **8j, a tile-ordered pack drawn whole: its per-run costs**: what the representation spent a run (the yeast's
    2.29 M at 16³) without a pass's own need: the paged culler's pages' members a stamp a member, not a set a page, a
    page's box read once for its members (`growBoxes`), the members touched flags, not a set; the page store's add its
    runs in one scratch by page (in particle order: a page's consecutive), its members' instances an array by member,
    not a map; the member pool's target lookup a table where ids are small; the tile set's targets a tile its runs'
    where ascending (a pack's), no set or sort; the catalog's compartments' entities a pair a flag of a table, its sets
    a pair, not a run; the shared serial keys made when a tile's list is first asked; `ParticleValues.map` of runs a
    loop (not `Int32Array.from` with a mapping), a run's value as the last one's reused; the entity theme's index by a
    table, not a map a run. Measured (Chrome, the yeast whole): the representation 1.58 -> 1.18 s, drawn 3.7 -> 3.0 s
    (by target 2.2 s); the views from outside, close and inside the same pixels. Left: the pool's empty members made
    and then filled (5,639: about 170 ms), the store's add (about 230 ms, a run's callback), the upload. Tests:
    `page-store.test.ts` (a target past the members' table, freed and back), `particle-tiles.test.ts` (a compartment's
    entities, ids small and past the table), `particle-list.test.ts` (equal runs side by side mapped once).
  - **8j, the paged cull an index a visible instance**: as the merged culler, its regions of instance indices laid out
    on the CPU (a block a level, a member's room its instances; a level's raster region the same room from its end: a
    particle at one level, drawn by one), so one pass a phase, no count and prefix sum before nor a write pass doing the
    count's work again; a particle's pose composed on the GPU when its page is written (`paged-cull-compose`, a slot a
    particle of its page's range, at the next cull), the draws and the raster reading it by the index (the
    raster's paged branch gone). The cull reads a particle's run and position as before, not its pose (tried: a pass
    3.0-3.3 ms against 2.6-3.0, 32 B a particle against 8). Measured (the yeast's tile-ordered pack drawn whole, Chrome,
    RTX 4070 Ti): a frame's GPU orbiting 15.8 -> 13.6 ms, zooming 15.5 -> 13.5 ms; the cull 11.2 -> 5.9 ms a frame
    (count and write a phase before), the raster 4.4 -> 6.8 ms (its threads a slot of every level's block, 126 M, as
    the merged culler's); by target 10.3 and 10.5 ms (its cull 3.3 ms: a cell's instances skipped whole). The images
    the same pixels. Memory: the poses 32 B a particle (1 GB the yeast), the indices 4 B a particle a level (0.5 GB),
    the positions and rotations kept (16 B a store slot, 0.8 GB); the poses' room twice the instances before (2 GB).
    The streamed views (16³, the following view): the fill as before (`writeBuffer` 1.1 s, frames p90 18-34 ms), zoom
    frames p50 22 -> 21 ms, orbit 62 -> 59 ms but its p90 80 -> 107 ms (10 long tasks against 7: not attributed). The
    first measurements of this step and of 8k ran with the GPU's memory short (other processes' 6.8 GB): 2-4× slower
    and misleading (a pose read 31 ms, the fill worse); measured again after a restart. Left: the raster's threads a
    level's raster instances alone (2.4 ms a frame), the orbit's tail. Tests: `webgpu-paged-cull.test.ts` (the
    regions' slots' poses as the CPU composes them, the ghost's, within the device's limit).
  - **8j, the sphere raster a workgroup a block of its instances**: a merged object's raster regions (the merged
    culler's and the paged one's) dispatched by their records' blocks of 64 instances (`sphere-raster-blocks`, their
    firsts a prefix sum, the dispatch indirect from the total), not a thread a slot of every level's regions (126 M for
    the yeast, 2.7 ms a phase with nothing to draw). Measured (the yeast drawn whole, Chrome, RTX 4070 Ti): the raster
    pass 3.5 -> 0.85 ms a phase; a frame's GPU tile-ordered 13.6 -> 8.0 ms, by target 10.3 -> 5.1 ms. The images the
    same pixels tile-ordered; by target within its own run-to-run noise (2.7k pixels of 926k between two runs before,
    3.1k after). The streamed views' orbit with it as before the index (frames p90 75 ms, 4 long tasks against 7-8;
    zoom p50 17.7 against 22 ms): the tail the raster's threads. Tests: the native suite (the raster's merged objects).
  - **8j, the following view's layers by range**: an update's colors (its tiles' pages, scattered over the store's
    layer) uploaded whole each step: the changes' ranges past 8 handed out as none (`TextureChanges`), the yeast's 268
    MB of colors written 90 times an orbit. Past 64 now merged across their smallest gaps (`mergeSmallestGaps`; half the
    layer still all). Measured (16³, the following view, Chrome): orbiting frames p50 60 -> 21 ms, p90 75 -> 24 ms, no
    long task (4 before), `writeBuffer` 1.75 -> 0.4-0.7 s of the moving phase; zoom and the fill as before. Tests:
    `texture-changes.test.ts` (scattered ranges each, many merged across their smallest gaps, all covered).
  - **8j, the paged cull a group's box, a tile's order kept**: a cull workgroup's particles (256 of a page) tested
    whole by their centers' box (padded by the largest reach), against the frustum and the pyramid, not the page's (a
    tile's); the boxes reduced by the compose pass (a workgroup a group, no CPU), the regions' firsts and the first
    phase's counts one buffer (a binding of the stage's eight). Measured (the yeast drawn whole, tile-ordered): a
    frame's GPU 8.0 -> 7.4 ms; the images the same pixels. Tried, a tile's particles by sub-cell then target (in
    memory, 16³ tiles, no priming): 2³ sub-cells a frame 10.0 -> 9.5 ms, 4³ 9.2 ms, the build 4.3 -> 6.0 and 7.5 s (a
    target a run a sub-cell: 2.5-5× the runs); not taken: the gain under a ms, the cost every tile's runs. The
    tiles stay 16³ (the streaming's unit: coarser tiles fewer runs, but fewer and larger loads). Left: the gap to
    the merged culler (5.1 ms by target: its cells skipped whole). Tests: the native suite (`webgpu-paged-cull`: the
    regions as the CPU classifies, a pyramid, the clip objects).
  - **The grid by default its data's**: `ParticleTileSetFromList`'s `cells` "auto" (the default) or a power of two
    (1-64): a tile-ordered list's written grid (`ParticleTiles.getAutoCells`: its tiles views of its ranges; another
    grid, or a non-power-of-two coarser one, orders every particle again: the yeast's tile set 3.8 s and a copy, not 0.4
    s), else the power of two of about `AutoTileParticles` (8192, a store page) a cell: the yeast 16. Not coarser for
    fewer runs: the tiles are what streams. Tests: `particle-tiles.test.ts` (the grid by count, a directory's), the
    plugin's `tiles.test.ts` (the default the automatic grid, named in the description).
  - **8f done, a tile-ordered pack read by range**: the model's `pack-stream.ts` (`ParticlePackSource`: bytes, a Blob by
    slice; `batched`: six reads in flight, ranges within 256 kB one read; `openParticlePackStream`,
    `createParticleTileSetFromPack`), the pack reader split to share it (`readParticlePackHeader`,
    `getParticlePackInfo`, `createParticlePackSitesTarget`); the plugin's `ParticleTileSetFromPack` (from Root: a URL by
    HTTP ranges, a file by slices) and a format's `stream` hook tried by `OpenFiles` and `DownloadFile` before a read or
    download (the particle pack's: a tile-ordered pack opened by range with a view following the camera and its targets'
    representation; else, or a server ignoring ranges, read whole); the dev server answering single ranges. Measured
    (the yeast's tile-ordered pack from a local server, Chrome): opened in 0.16 s reading 22 MB of 731 (its directory,
    opacities and runs; Node), its first tiles drawn after 1.3 s, the following view filled (1,262 tiles, 6 M particles)
    after 12 s, 542 MB fetched in 4,632 responses, the heap 1.2 GB at most (whole: 1.1 GB, the file held). Tried
    unbatched: a view's targets (thousands of site ranges at once) refused by the browser (ERR_INSUFFICIENT_RESOURCES),
    nothing drawn. Left: fewer requests (a tile's blocks apart, a read each), the heap, a tile let go its bytes too.
    Tests: `pack-stream.test.ts` (a streamed set as the whole read's, its tiles' and targets' data, read when asked,
    once a target; a blob; a pack by target refused; fibers not read; batched reads merged within the gap, at most so
    many in flight), the plugin's `pack-stream.test.ts` (a tile-ordered pack opened by range with its view and
    representation, a pack by target whole).
  - **A camera view's particles capped**: the camera pick's `maxParticles` (0, the default: none), its tiles by priority
    (the least occluded, the largest on screen, the nearest) while within it, at least one. Measured on an iPhone 16e
    (Safari, the streamed views from a local server): the Mycoplasma cell's tile-ordered pack (2.1 M particles, 848
    targets, 2.9 M sites; 194 MB of heap and 163 MB of GPU buffers on desktop) drawn; the yeast cell's crashing the page
    at any cap tried (3 M particles: 578 MB of heap, 421 MB of GPU buffers): its targets the floor (5,639, 20 M sites,
    nearly all of them in any few hundred tiles), not its tiles. Left: targets drawn coarser for a phone (done:
    coarse sites, below). Tests: the plugin's `tiles.test.ts` (the capped pick the uncapped's first tiles within the
    cap, one below a tile).
  - **Coarse sites**: a pack's tile set's `siteStride` (1, 10, 100): a target's sites every stride-th, read and copied
    (`ParticleSites.coarsen`), drawn as that stride's level of detail near the camera too (`SphereSites.stride`, the
    levels' `base`). Measured (the yeast, desktop): a tile of 6.8 k particles 2,058 targets and 9.5 M sites; at 10 its
    heap 321 -> 253 MB, at a 300k cap 419 -> 314 MB, its textures 323 -> 113 MB. Its representation's members the rest
    (a pool member a target of the catalog, about 8 KB; a target's own about 16 KB). Tests: `particle-list.test.ts`
    (`coarsen`), `pack-stream.test.ts` (a target every stride-th), `lod-levels.test.ts` (coarse levels).
  - **Coarse levels**: a tile-ordered pack's octree (formats.md: `tiles.levels`, its nodes copies of a sample of their
    tiles' particles), read as tiles of the set (`ParticleTileLevels`, a node's particles drawn `2^L` larger: the page
    store's runs by `ParticleList.scale`, the sphere raster's radii by the pose's scale). A camera view picks a cut of
    the octree (`pickNodes`): a visible tile its coarsest node all of whose visible tiles want it (`detailPixels`, 4:
    a tile's largest molecule's radius on screen from which it is all drawn, the sample `(size / detailPixels)³`
    below); under `maxParticles` the detail lowered by bisection, at most 4 times (`MaxDetailRaise`: molecules drawn
    larger near the camera filled a phone's screen many layers deep), then the least important left out. Tried first:
    levels a tile each, not a region's: a page store page (8,192 slots) a level of a tile however few its particles,
    the yeast at 300k 1,271 tiles, 309 MB of storage. Measured (the yeast, 1600×900, still, the steady state): uncapped
    the same tiles and frames as without (13.2 M particles at 4 px, 1.6 ms of GPU whole, 2.6 ms near); at a 1M cap
    996 k in 150 nodes, 0.9-1.8 ms, heap 0.55 GB and GPU 0.67 GB (from 1.05 and 2.0); at 300k 293 k in 77 nodes but
    6.2 ms (the sphere raster 3.2 ms, the pyramid built 11 times a frame); the pick's cut 0.1-1.5 ms (the visibility
    15-25 ms). On an iPhone 16e the yeast whole and zoomed, borderline. Left: step 9 (the coarse nodes' GPU time, a
    phone's defaults, the members, the look).
    Tests: `pack-stream.test.ts` (a node a sample of its tiles' particles, copies, its slots and scale; none without),
    the plugin's `tiles.test.ts` (the cut: a visible tile shown once; the root where small; within a cap),
    `page-store.test.ts` (a list's scale).
  - **A following view's steps**: a step whose pick is all new and none of it ready keeps the view as it is (it was
    empty: a blank frame, and the target representation built again twice, its color mode by member lost with an empty
    view's theme and back: the yeast's store, pool and object twice at each zoom step, the old alive until committed);
    an empty slot list keeps the representation's color mode. A capped view's loader lets go of the ready tiles neither
    wanted nor shown past twice the cap (its targets kept: one read again a new object, the members all made again).
    Measured: a zoom tour's resets 7 -> 1, its heap 329 -> 387 MB (was to about 500). Tests: the plugin's
    `tiles.test.ts` (never empty; the tiles passed let go).
- Written without the tile spec it serves (not in this repository): the model below is
  generic, its open questions are where the spec decides. The plugin state drives it (which tiles, when); this plan
  is the model below it, and example transforms to test it.
- Decided: a regular grid (no octree), its particles ordered by a fine grid's cells in Morton order (Design 10: a
  power-of-two coarser grid's cells ranges too, the tile size chosen when the set is made); the memory measured and let go by the app, the model giving what it
  needs to know (a tile's particle count and bytes); a particle in the tile of its center, in one tile only; a tiled
  trajectory's tiles fixed, its frames a tile's poses; a light catalog up front, the targets' data on demand.

## Why

- **A model larger than memory or than one download**: a spec elsewhere subdivides a mesoscale model into tiles, each a
  part of the space with the particles in it. A viewer loads the tiles it needs (in view, near the camera), when it
  needs them, and lets others go.
- **The targets are the model's, not a tile's**: every tile places particles of the same targets (proteins, lipids:
  their structures or sites). A tile must not bring its own copy: the targets' geometry, its GPU data (the yeast cell's
  sites: 163 MB) and its programs are built once and drawn for the particles of every loaded tile.
- **What exists**: a `ParticleList` holds particles and their targets (`targetMapping`, `targetInfo`, entities,
  compartments); a target representation draws each target once, instanced at its particles (`target/
  representation.ts`: a visual a target id, its geometry kept while the target is the same object). A
  `ParticleTrajectory` serves frames on demand (`getFrameAtIndex`: a list or a task); a frame of the same particles with
  other poses (`Particle.withPoses`) updates the targets' transforms alone (`instance-index.md` step 4). Nothing
  subdivides a list or loads part of it.

## Design

1. **The catalog, light and up front** (`ParticleTileCatalog`): what tiles share and what is known of every target before
   any is loaded: its id, kind (structure, sites, shape, volume), name, entity, compartment, its radius (how far its data
   reaches from where a particle places it: what a tile's visibility needs), and what decides its
   programs (sites quantized, their elements known: the precompile's signatures, `getSitesSignatureDefines`); the
   entities and compartments. Kilobytes for thousands of targets: themes, legends, the precompile and the keys from it,
   their palettes over all entities from the start (no color shifting as tiles load). Why not all targets' data up front:
   the first view would wait for every target (the yeast cell's sites, ~160 MB, most out of view) and memory would hold
   all of them; why not each tile bringing its targets: the same target from two tiles would have to be merged by id (one
   object, or its geometry and GPU data twice), the catalog would change as tiles load (palettes, the precompile late)
   and a tile let go would need counts of who uses its targets.
2. **The targets' data on demand**: a target's data (its structure, sites, shape or volume) loaded once, when the first
   tile using it is wanted (`getTarget(id): ParticleTarget | Task<ParticleTarget>` of the provider), and the same object
   from then on: a tile's list carries the loaded targets themselves, so whatever is keyed by a target (its visual, its
   geometry, its GPU data) is one for all tiles. A tile is ready where its particles and the targets they use are
   loaded: a view holds ready tiles only, no target without data drawn (no placeholder). Which targets no loaded tile
   uses is the model's to tell (`unusedTargets`), letting them go the app's.
3. **The grid** (`ParticleTileGrid`): an origin, a cell size and the cells a dimension; a tile a cell, its id the
   cell's linear index (`x + y nx + z nx ny`), its box the cell's (a trajectory's: where its particles go in any frame,
   so boxes of neighbours may overlap). A particle in the tile of its center (in the first frame where frames move it):
   in one tile only, never repeated in a neighbour.
4. **The tile set** (`ParticleTileSet`, model `particles/particle-tiles.ts`): the catalog, the grid, per tile its
   particle count, the targets it uses, its opacity where known (per block of its cell the optical depth a length of
   its content: the occlusion's, 6) and an
   estimate of its bytes where known (what the app measures memory by; a
   target's bytes alike), the frame count (1 for no trajectory), and a provider: `getTile(id, frame)` (a
   `ParticleList` or a task of one, as `ParticleTrajectory.getFrameAtIndex`) and `getTarget(id)`: sync where in
   memory, a task where fetched or decoded; a tile's particles' targets by the catalog's ids. A tile's frames are the same particles in the same order (their poses only), so a frame step is a
   pose change, not another list.
5. **The loader** (`ParticleTileLoader`): the tiles' states (none, loading, ready, failed) per frame and the targets'
   (none, loading, loaded, failed), loads at most a few at a time, a tile's targets with it (a target once for all the
   tiles wanting it), a load no longer wanted cancelled (the task's abort); `want(ids, frame)`, `release(ids)`,
   `releaseTargets(ids)`, a `changed` event, and what is loaded with its counts and bytes. No budget of its own: what to let go is the app's (it measures
   memory). Playing a trajectory: the next frames of the loaded tiles wanted ahead (`want` of the frame after).
6. **View-dependent visibility** (`ParticleTileVisibility.get`): which tiles would show something to a camera, so that
   only those are loaded. A view given as numbers, not a graphics camera (the model is below graphics): the frustum's
   planes (`Frustum3D`), the camera's position and plane (the view direction), the pixels an Å makes at a distance of 1
   (the projection's focal scale times half the viewport's height, as the levels of detail measure it; orthographic: a
   scale alone) and the smallest size on screen worth drawing (pixels). Per tile:
   - **its reach**: its box grown by the largest radius of its targets (the catalog's): a particle at the box's border
     draws its target's atoms beyond it;
   - **in view**: its reach intersects the frustum;
   - **its size on screen**: its largest target's radius in pixels at the reach's nearest point to the camera; below
     the smallest worth drawing nothing of it shows (sub-pixel), even in view;
   - **occluded**: hidden behind the tiles in front of it, their content's optical depth gathered along the rays. The
     tiles in view, front to back (by their reach's nearest distance): each tile's reach tested against rays through
     the pixels' corners of a small grid of the camera's (64², on the CPU), then, where it is wanted, its blocks added
     to the rays through them (their density times the ray's path in them, from the near plane; a tile's blocks all
     within 1.5 times of each other one box of the least; nothing added to a ray opaque already). A ray opaque where it
     reaches an optical depth of 7 (under 0.1% passes), from the deepest back of the blocks it crossed (what is drawn in
     them may be anywhere in them); a pixel where all its corners' rays are, from the deepest of them; a reach hidden
     where all its pixels are opaque nearer than its nearest depth. So only wanted tiles occlude, tiles that each let
     much through hide together, and a tile's reach within the blocks that made a ray opaque is not hidden by them (the
     layer behind a dense front kept). A heuristic: the densities a mean (holes finer than a block in it only), a gap
     narrower than the rays' spacing not seen;
   - **a priority**: its size (larger first: what is near and large shows most), and its distance.
   What makes a tile not occlude, from the caller: the clip objects (its blocks cut by a plane to their largest
   axis-aligned part on the kept side, exact for an axis-aligned plane, shrunk along the normal otherwise; cut by
   another shape: no occlusion), the tiles drawn not opaque (an entity of theirs hidden, transparency). A camera inside
   the model, any direction, perspective at any distance: the rays are the camera's own, no set of directions. A
   trajectory's opacity one for all frames (each block's least) or a frame's.
   Out: the tiles that would show something, by priority, with their sizes and distances. What to want and what to let
   go is the caller's: margins against flicker (want with the frustum grown and the size lowered, let go only beyond
   wider ones: hysteresis), prefetch along the camera's motion, the memory. Not here: levels (a coarse stand-in of a
   far tile: the spec's, a grid of one level for now); the depth of what is drawn (the Hi-Z pyramid: a GPU query, finer
   than densities, later where they leave much loaded that is hidden).
7. **The view** (`ParticleTileView.create`): one `ParticleList` of the ready tiles at a frame, the loaded targets, the tiles'
   particles concatenated (their runs by target kept: a pack's tile a range a target), with
   - **keys**: a particle's identity, its index in the whole (its tile's base, the counts of the tiles before it, plus
     its index in the tile; `ParticleList.keys`, int32: the counts known with the set): the view's indices change as
     tiles load and go (what is concatenated before a particle); what is kept by particle (a selection, marks, a focused
     particle, a label) is kept by key, so it survives loads;
   - **its composition**: per target the tiles and their ranges it is made of;
   - **its lineage**: for a frame step the same composition with other poses (`Particle.withPoses`: `posesFrom`); the
     view a list is on the list (`ParticleTileView.of`), the previous one the caller's (each holding the one before
     would hold every view).
   The target representation draws a view as it draws any list, a visual a target: the targets' geometry shared by all
   tiles.
8. **Updates by composition**: a tile loaded or let go changes the targets it has particles of; the others' instances
   are the same particles (moved in the view's indices only). The representation, given a view of the same lineage,
   keeps the visuals of targets whose composition is the same (their transforms, colors, the culler's layout of them)
   and updates the others (a check by the composition, not the indices). A frame step keeps every composition: the
   pose-only path as it is. Adds and removes within a target's instances: `instance-index.md` step 4's open part, and
   the editable particle lists' plan: until then a changed target's instances built again.
9. **Plugin state** (examples, to test): `ParticleTileSetFromList` (a list subdivided by a grid of `n³` cells, by its
   particles' centers, its catalog from the list's targets; its provider slices the list and gives its targets, with an
   optional delay as a fetch would have; from a trajectory,
   its frames sliced alike), and `ParticleTilesView` (from a tile set: the tiles wanted, all or those in a box or by id,
   and the frame; their loads awaited, a view out as `SO.Particle.List`: the representations as they are); a mode by the
   camera (the tiles `ParticleTileVisibility.get` gives for the canvas's camera or a snapshot of one; updated as it
   moves: a behavior, later).
   The memory is the app's.

10. **Tiles as views** (Steps 8; the pack's principle, [particle-pack.md](particle-pack.md): the list's storage,
    blocks aligned for `writeBuffer` from the file's bytes, read as views). The whole list's path visits no particle
    on the CPU: its poses composed on the GPU from the pack's grid words and snorm16 rotations, a target its particles'
    range (`ParticleRuns`), its bounds the pack's. A tile is to be the same: a range of every block.
    - **Tile order**: the particles ordered by tile, and within a tile by target (a tile's targets a few runs, as a
      pack's are); a tile then `[first, first + count)` of the positions, the rotations and every value block, and its
      targets' runs a slice of a runs table (target, first, count). A list's particle index is its order: a pack
      written tile-ordered has those indices (its keys; no keys stored, as now); a list ordered in memory keeps its
      source's index as its key (4 B a particle: picking, labels and `particle-index` colors as the untiled list's).
    - **The tile directory** (with the list, the pack's where it holds one): per tile its range, its box (of its
      particles' centers: the pack's bounds kind), its runs' range, its byte count; per run its target, first and
      count. The tile set made from it: no pass over the particles (`fromFrames` today 2 s for the yeast).
    - **A tile's list a view** (`Particle.range`): subarrays of the grid words and rotations, its targets, entities,
      compartments and radii the runs' slice (a pack's per-target values are runs already), its bounds the directory's;
      no gather, no keys array, no copy kept per tile (`firsts` today: the tiles' copies held for the set's life).
    - **The page store from the view**: a store page's values the view's subarrays (as now, the copies gone), its box
      the directory's tiles' (no decode of its words: `writeBox`), its runs the directory's.
    - **Colors by member**: a theme whose color is a target's (or its entity's, compartment's: constant over the
      target's runs) a color a member, read in the shader by the instance's member (the merged objects' `object` color
      type: a table read by the segment), no per-particle layer; only a per-particle theme (`particle-index`, an
      attribute) fills a layer, by range. Markers' aggregates kept as marks change (no scan of the capacity an add).
    - **Steps by time**: a follow step adds tiles until its CPU time (not its particles) reaches a budget (about 4 ms a
      frame), the uploads by bytes: the fill bounded by the upload's bandwidth, frames drawn between.
    - **Levels without an octree** (decided): a pack ordered by its tiles in Morton order holds every power-of-two
      coarser grid's cells as ranges too (8³ from 16³), a coarse tile the union of its tiles' runs; the tile size
      chosen when the pack is written, not a finer grid (a fine cell's targets runs of their own: 2.6-5× the runs at
      16³, 8e).

## Steps

1. **The model**: `ParticleTileCatalog` (from a list's targets: `ParticleTileCatalog.fromList`), `ParticleTileGrid`,
   `ParticleTileSet`, a provider of tiles and targets in memory (`ParticleTiles.fromList`, `fromTrajectory`: subdivided
   by a grid). Tests: every particle in one tile, by its center; a tile's targets; a tile's list valid; a tile's frames
   the same particles with other poses.
2. **The loader** (done): states per frame, the targets loaded once for the tiles using them, the concurrency limit,
   cancellation, `release`, `releaseTargets`, `getUnusedTargets`, `subscribe`, the loaded's counts and bytes. Tests with a
   provider of delayed tasks (fake timers): order, a target loaded once for two tiles, cancellation, release, a frame
   wanted ahead.
3. **The visibility** (done): the catalog's target radii (from a list's targets: their data's bounds around their
   centers), `ParticleTileVisibility.get` with the occlusion by optical depth (front to back, the clip objects and the
   tiles not occluding). Tests: a tile out of the frustum, one in it but sub-pixel, the order by size on screen, a
   tile's reach past its box, orthographic; a tile behind a dense one left out, behind one thin kept, behind two thin
   left out, behind a hole kept, behind neighbours meeting at a seam left out, behind one a plane cuts away kept,
   behind a tile not drawn opaque kept; the camera inside.
4. **The opacities** (done): computed where the spec gives none (`ParticleTileOpaque.computeOpacities`, offline: the
   pack writer), by the content's density. Measured: the yeast cell in tiles seen from outside whole, zoomed and from
   inside, the tiles kept with and without the occlusion (Status).
5. **The view** (done): the concatenation with its keys, composition and lineage. Tests: a view of tiles equal to the list's
   particles of them; its targets the same objects in views of other tiles; the keys of a particle the same in views of other tiles; the composition per target; a frame
   step's view the same composition (`isSamePoses`).
6. **Plugin transforms** (done): `ParticleTileSetFromList`, `ParticleTilesView` (by id, box or camera, the clip objects). Tests: the state's
   objects, a view's representation drawing the targets once (headless), a frame step, the camera moved.
7. **The representation by composition** (done; graphics): the visuals of targets of the same composition kept. Measured:
   the yeast cell in tiles (e.g. 4³), a tile loaded or let go: the update's time against a full one.

8. **Tiles as views** (designed; Design 10). In order, each measured in Chrome (the yeast, 8³ and 16³: a known tile's
   add, the whole cell's first fill, frames while filling):
   - 8a. **Tile order in memory**: `ParticleTiles.fromList` orders a list by tile and target once (one counting sort:
     the grid words and rotations gathered once, the values as runs per tile and target, the source index the key),
     builds the directory in the same pass, and gives a tile's list as a view. No format change; `subset` and the
     tiles' copies gone. Tests: a tile's view the same particles (positions, rotations, targets, radii, keys) as
     `subset` gives; the slot view's frames as before.
   - 8b. **The store from the directory**: a page's box and runs from it (`writeBox` gone), a page's values its view's
     subarrays. Tests: the store's records as computed from the particles.
   - 8c. **Colors by member** for themes constant over a target's runs (entity, target, compartment), the layer only for
     per-particle themes; markers' aggregates kept, no capacity scan. Tests: the frames as the per-particle layer's.
   - 8d. **Steps by time and bytes** in the follow behavior. Tests: a step's time within its budget (fake clock).
   - 8e. **The pack tile-ordered** (done; `ParticlePackVersion` 6): the writer (`writeParticlePack`'s `tiles`, the
     `particle-pack` CLI's `--tiles <n>`) writes the particles by tile in Morton order and the directory as blocks
     (`tile.ranges`: cell id, first, count; `tile.boxes`; the targets' runs `target.runs`: first, target), the reader
     makes the tile set from them (no pass, no copy: every tile a view of the file's buffer). A pack without them: 8a's
     order in memory. Tests: a written pack's tiles the views 8a makes; versions 2-5 read as before.
   - 8f. **Ranges read on demand** (done; particle-pack.md's `Blob.slice` and HTTP ranges): a tile's bytes read when
     wanted, a pack larger than an `ArrayBuffer`.
   - 8g. **What a run costs, measured first**. A run is a target in a tile; ordered by target within a tile, a tile's
     runs are its distinct targets, the fewest any layout where a page knows its members has (the yeast: 2.29 M at 16³,
     14 particles a run; 0.72 M at 8³, 43). Where a run costs: a store record (CPU, the add: 1.6 ms a step at 16³),
     the cull kernel's binary search of a page's runs a particle (about 10 steps at 800 runs a page), its workgroups'
     keys by run (short runs: more keys a workgroup), the pack's runs block (8 B a run), the catalog's compartments
     over the runs (0.18 s: per target instead, a pack's values are a target's), a whole draw's runs walked a target
     (8e). Measure the paged cull's GPU time by pass (timestamps) and the store's add at 16³ and at 8³ (an 8³ pack),
     then, by what costs:
     - fewer runs: the pack written at 8³ (16³ tiles not then a view: an 8³ tile's runs, served as one tile or split
       by the store's pages), the streaming coarser;
     - cheaper runs: a page's run of a particle without a search (a byte a particle, a page's runs at most 256; the
       store's per-particle members, 4 B, where runs are denser, `RunsShare`), the workgroups keyed by member;
     - not joined: a run is a target, other targets' geometry other; joining needs targets drawn alike (none here).
     Tests: as the step's (the records the particles', the frames the same).
   - 8h. **Opacities in the pack**: a tile-ordered pack's tiles' opacities written with it (`tile.opacity`: per tile
     its blocks' optical density a length, `ParticleTileOpaque.computeOpacity` at write, 4³ blocks a tile: 64 values,
     the yeast's 2,753 tiles 0.7 MB as float32, less as float16 or a log byte), read onto the directory's tiles; a
     coarser set's tile its cells' blocks (8 of b³: one of (2b)³), so the set made from the directory has them, no
     pass (`ParticleTileSetFromList`'s `opacity` computed only where the pack has none). Tests: the read opacities the
     computed ones; a coarser set's its cells'.
   - 8i. **An occlusion threshold**: the view's `opticalDepth` (`ParticleTileVisibilityOptions`, 7 by default: under
     0.1% of light passes) a parameter of the camera tiles (e.g. 1-10, labeled by what passes: 37% at 1, 5% at 3,
     0.1% at 7): lower leaves out more tiles behind thin content (faster, sometimes a gap), higher only behind dense.
     Its key in the pick's cache. Tests: a tile behind a thin one kept at 7, left out at 1.
   - 8j. **The fill's frames** (about 90 ms at p90 while the whole cell fills, 17.5 ms after): the GPU's time by pass
     over the fill (the store's uploads, the culler's sync, the occlusion pyramid, the draws), then what dominates;
     with it `stepMs` again. And a tile-ordered pack drawn whole by the paged path (a tile a range: 5.3 s -> about
     2 s, one path for both).
   - 8k. **Occlusion without tiles**: the opacities as precomputed occlusion for any list drawn culled on the GPU,
     tiled or whole, by the paged path or the whole list's (its instance cells). Measured first: the CPU's occlusion
     pass over the yeast's 2,753 tiles (`ParticleTileVisibility`) 16-18 ms from outside, 3-6 ms near (its blocks' adds,
     not its resolution: 16-18 ms at 64², 32² and 16²): not a frame's, and throttled its stale picks holes while moving
     (no loading to hide them), so not the CPU. A page tested whole in the cull kernel against the pyramid (its box's
     sphere padded by the largest reach times the largest scale, `uMaxReach`; a page bit after the particles' in
     `occludedBits`: the first phase's occluded marked, the second's retested, all its particles where no longer)
     changes little alone (the yeast whole at 16³: orbit 79.6 ms of GPU as 79-82, zoom 58.5 as 64-66): a page's box a
     tile's (2,700 Å, its sphere reaching in front of the surface) and the pyramid's farthest depth seen through the
     gaps between molecules (a gap's background a region's depth: nothing behind it occluded). The opacities model
     those gaps (how little comes through, not a closed surface), so:
     1. **The pyramid primed from a density grid**: a compute pass a frame, a ray a pixel of the pyramid's base
        (a few hundred a side) through a grid of optical densities (64³ for the yeast: the tiles' 4³ blocks side by
        side), the depth where its optical depth reaches the view's threshold (`opticalDepth`, as the tiles' picks)
        merged into the pyramid (the nearer): from the first frame and after a camera jump, the gaps closed, so the
        page test above and the instances' own (the whole path's cells too) leave out what is hidden. Cheap (about 40k
        rays of up to about 128 steps), no CPU, no lag; conservative as the tiles' occlusion (a density a lower
        bound).
     2. **The grid**: a tile-ordered pack's tiles' opacities (`tile.opacity`, side by side); else a pack's own
        density grid over its box, written with it whether tile-ordered or not (`ParticleTileOpaque.computeOpacity`
        over the grid's cells, about 1 MB at 64³, a few seconds more to write).
     3. **The pack writer** writing that grid for packs by target too (version 6, an optional block).
     4. **Lists without one** (an mmCIF assembly, a list made in memory): computed once at load (a pass over the
        particles, about a second for the yeast), opt-in or off the main thread.
     Measured: the yeast whole and tiled, from outside and zoomed in, the frames' GPU against the two-phase occlusion
     alone, what the pages and instances leave out. Tests: the pyramid's depth where the grid is dense (a ray reaching
     the threshold), none where clear; a page behind a dense grid left out, one seen through a thin one kept.
   Expected: a known tile's add its bytes' upload and a few records (under 1 ms of CPU for 100k particles), the whole
   cell's fill bounded by the upload of 16 B a particle (31.4M: 500 MB), the tile set made in milliseconds from a
   tile-ordered pack.
9. **A whole cell on a phone** (coarse sites and coarse levels done: Status). What is left, each measured as Status's
   are (the yeast, desktop and an iPhone 16e):
   - 9a. **Coarse nodes' GPU time**: the yeast whole at a 300k cap 6.2 ms of GPU a frame against 1.6 ms for 13 M
     particles uncapped (still, the steady state): the sphere raster 3.2 ms, the occlusion pyramid built 11 times a
     frame (about 2.3 ms; not in the uncapped frames). By pass and moving: a scaled molecule's spheres past the raster's
     size (impostors), their overlap, the occlusion a coarse view gets. What a phone pays most.
   - 9b. **A phone's defaults**: `siteStride` and `maxParticles` chosen for the device, not by hand (Safari gives no
     `deviceMemory`: by a touch screen's size, a setting, or the memory a view takes measured as it grows); the render
     targets at a phone's pixel ratio (3: about 2.8 times a 1200×900 canvas's pixels) in its budget.
   - 9c. **The representation's members**: a pool member a target of the catalog made up front (the yeast's 5,640: 45-80
     MB of heap, `createEmptyMergeMember`'s cells) and a target's own member (about 16 KB): lighter empty members (cells
     shared until filled), or a pool of the targets in view with headroom (made again past it).
   - 9d. **Targets let go**: a capped view's loader keeps its targets (a target read again a new object, the
     representation's members all made again, `isSameTargets`): a target's identity kept across reads (its id), so its
     sites can go with the tiles that used them.
   - 9e. **The coarse levels' look**: a coarse region lighter on screen than its neighbours (the volume kept, `2^L`, not
     the coverage); a level's scale by eye (an exponent, or its coverage measured), screenshots at each level.
   - 9f. **Levels on the desktop**: the pick's size a tile's largest molecule (rarely under 4 px at desktop sizes: the
     yeast whole uncapped all its tiles), so levels change nothing there; a typical molecule's size (a tile's median
     reach), so far regions draw coarse and faster where the eye cannot tell.
   - 9g. **A region replaced while it loads**: a step keeps the view only where none of its pick is ready; a node
     replacing its tiles (or tiles their node) partly ready leaves the rest of the region blank until it is. The shown
     kept per region until what replaces it is ready.
   - 9h. **The other paths**: levels only from a pack written with them (a list in memory, `ParticleTileSetFromList`,
     has none: its nodes sampled in memory, its far view reading every tile); `ParticleList.scale` drawn by the page
     store alone (the paged WebGPU path: a tile set's views elsewhere draw nodes unscaled); coarse sites drawn as a
     level only from quantized sites (`Spheres.fromSites`: WebGL1's builder path draws them sparse).
   - 9i. **What levels read up front**: the nodes' directory and runs read at once with the tiles' (not measured; the
     yeast's node runs about 1.3 M); a node's runs read with it instead, its targets from its directory entry.
   - 9j. **The occluder grid and the clip objects** (done): the pyramid primed from a density grid (8k,
     `occluder-grid.ts`) marches each texel's ray through it, a sample the grid's object's clip objects remove adding
     nothing (`hi_z_occluder_wgsl`, `clip_storage_wgsl` on the object's packed clip objects, `HiZPass.setOccluder`'s
     `clip`): any shapes, a cut moved (a drag) its uniforms alone, the grid not uploaded again; on with clip objects
     (`isTargetOccluderOff`: as asked only). An object moved (a preview transform) primes nothing (its grid at rest).
     Tests: `webgpu-hi-z.test.ts` (a grid cut by a plane, clipped whole by a cube, the cut moved).
   - 9k. **A cut moved at the yeast's scale** (a clipping gizmo's drag, a few clip objects): the store's object takes a
     clip-only update (its uniforms and its cull's clip objects, no page written again: `TargetSlotsObject.updateClip`);
     what draws a set registered by the representation itself on every update, a drag's preview too
     (`ParticleTileDrawings`, the model's), the follow behavior told per set (at most every `intervalMs`, once more as it
     settles); a view's step keeps the representation's props (a data-only update: the preview not undone); while the
     drawings change the pick's occlusion at 32² (`ParticleTileDrawings.SettleMs`): inside the yeast at 16³ a sphere's
     28 → 11 ms, from outside about the same (19 ms: the tiles' own work, a block about a ray there). Left: the pick
     off the main thread or the tiles' work cheaper where many are in view (a step a frame then), measured frames
     of a drag in the app.

## Open questions (the spec's)

- **The occlusion's time where many tiles are in view** (the yeast whole at 16³: 15 ms a call): a tile's blocks
  marched by each ray (one slab test a tile, not a block), fewer rays (32²), or the call throttled to the camera's
  rest. Also a tile's reach by its own particles (their spheres' box at tiling time) instead of its largest target's
  (650 Å in nearly every yeast tile).
- **The opacities in the spec**: a tile's blocks' densities given, or computed by its reader (from the catalog's
  volumes: the particles' positions alone).
- **The targets' radii in the catalog**: given by the spec (proposed: each target's bounding radius around where a
  particle places it), else known only as their data loads (a tile's reach its box alone until then: its border's
  particles' atoms may show before their tile is wanted).
- **Particles appearing or going between frames** (a frame of a tile with other particles): not a pose change; a
  tile's frames then other lists, its targets' instances built again (or the editable lists' adds and removes).
- **The catalog in the spec**: a file of its own, or the tiles listing their targets' ids (the light catalog gathered
  from them); and whether it gives the targets' bounds and kinds (else known as their data loads, the precompile late).
- **Keys from the spec**: where its particles have ids of their own, the keys those (proposed: the index in the whole
  where it has none).
