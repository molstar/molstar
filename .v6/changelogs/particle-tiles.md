## Particles: particle tiles

- Particle tiles (`model/particles/particle-tiles.ts`): a particle list or trajectory subdivided on a regular grid, a
  light catalog of its targets up front, a tile's particles and the targets' data served on demand (sync or a task); a
  loader (`particle-tile-loader.ts`: the tiles wanted by priority, a few at once, a load no longer wanted cancelled, a
  target loaded once for all tiles); views of some tiles as one particle list (`particle-tile-view.ts`: keys kept
  across views, per target the tiles it is made of, a frame step the same particles posed otherwise).
- The tiles a camera would see something of (`particle-tile-visibility.ts`): in view, not sub-pixel, not hidden behind
  the tiles in front by their content's optical depth gathered along the camera's rays (`particle-tile-opaque.ts`: a
  tile's opacity per block of its cell from its molecules' volumes). The yeast cell in 16³ tiles, zoomed: 77 of 437
  tiles in view kept (1.1 M of 5.5 M particles).
- State transforms with the particles (examples to test them with): `particle-tile-set-from-list` (a list or trajectory
  in tiles, an optional delay as a fetch would have, the opacities) and `particle-tiles-view` (all, by id, in a box or
  seen by the camera, at a frame: a particle list the representations draw as any).
- A view of other tiles keeps the instances of the targets whose tiles are the same (their transforms, colors, sizes):
  the yeast cell's 30 tiles in view, a tile added or let go 450 -> 290-390 ms (headless WebGPU). The `particle-index`
  color theme colors a view's particle by its index in the whole, so it keeps its color in views of other tiles.
- A tiles view as a slot list (`particle-tiles-view`'s `slots`): a tile's particles at slots of its own, so the target
  representation on WebGPU adds and lets go of a tile's instances in one object kept across views
  (`target/paged.ts`: a tile's particles uploaded as stored, culled and composed on the GPU), not the targets'
  instances built again.
- A tiles view by the camera follows it (`particle-tiles-view`'s `follow`, on unless a camera is given): the
  `particle-tiles-follow-camera` behavior (in the Viewer) picks its tiles again as the camera moves and adds them as
  they load, in steps of the view's `budget` of particles with frames drawn between (not all at once).
- A tile set's and its views' defaults for working with tiles: 8 cells an axis, a view of the tiles seen by the camera
  (following it), as a slot list. A following view's `stepMs`: its updates by time (off by default).
- A following particle tiles view adds up to 1M particles a step by default (was 150k): the yeast cell filled in
  2.5 s instead of 4.2 s, moving smoother. The hierarchy color theme colors a pack's targets by member (no per-particle
  colors).
- A particle tiles view by the camera adds the least occluded tiles first (what is seen fills before what shows through
  gaps), and loads the tiles the camera picks whenever the pick changes.
- WebGPU: tiled particle views and tile-ordered packs drawn whole prime the occlusion depth pyramid from their tiles'
  densities, so what lies behind dense content is culled through the gaps between molecules: the yeast cell's frames
  orbiting 80 -> 46 ms of GPU.
- Particle packs ordered by target store a density grid (`densityBlocks`, 64 by default; about 1 MB), and the target
  representation primes WebGPU's occlusion from it when drawing the list whole: the yeast cell's frames orbiting
  99 -> 41 ms of GPU. Lists without one compute it with the new `occluder` param's "computed where none" (opt-in).
- WebGPU: the occlusion's priming threshold is a canvas setting (`hiZ.occluderOpticalDepth`, 7 by default; lower is
  faster with occasional gaps), and a computed density grid is made in chunks without freezing the page.
- A tile-ordered particle pack opens read by range by default (dropped, opened or loaded from a URL whose server
  answers range requests): its tiles and targets read as a view following the camera wants them, never whole. The
  yeast cell's first tiles drawn after 1.3 s, 22 MB read up front of its 731 MB (a local server); other packs, or a
  server ignoring ranges, read whole as before. New: the `particle-tile-set-from-pack` transform.
- Particle Tiles: the grid's `cells` is "Auto" by default: a tile-ordered pack's written grid (its tiles read as they
  are, where another grid orders every particle again), else about 8k particles a tile; or a power of two (1-64). The
  tile set's description names the grid.
- Particle Tiles View: `maxParticles` caps the particles a camera view shows (0, the default: no cap), the tiles past it
  by priority left out: what bounds a view's memory on a phone. A capped view lets go of the tiles it passed.
- Particle packs: coarse levels (`writeParticlePack`'s `levels`, the particle-pack CLI's `--levels`): a tile-ordered
  pack written with an octree of samples of its tiles' particles (about a seventh larger), so a view by the camera
  draws far regions as a few particles drawn larger, read as one range each: the yeast cell whole as 7.7 k particles at
  the coarsest. A view's `detailPixels` (4) chooses where; within `maxParticles` it draws coarser before leaving
  tiles out. The pack's tile set's `levels` (on) uses them where written.
- Particle Pack Tiles: `siteStride` keeps every 10th or 100th of a target's sites (drawn larger, as a coarser level of
  detail): a tenth of the memory a view's targets take.
