# Plan: shared streaming

## Status

- **Current**: proposed, nothing built. Decisions open (see Open).
- One way to stream model data by view for every kind that needs it: particle tiles and pages today, streamed SOG
  next ([gl/gaussian-splats-lod.md](gl/gaussian-splats-lod.md)), then meso\*'s volumes, segmentations and meshes.
  Shared where the kinds do the same thing (reading, scheduling, caching, budgets, upload pacing), each kind's own
  where they differ (decoding, GPU layout, identity, what decides a want).
- Related: [particle-tiles.md](particle-tiles.md) (the loader, the visibility, the pack read by ranges),
  [particle-residency.md](particle-residency.md) (pages by the cull's wants), [gl/mesoscale-lod.md](gl/mesoscale-lod.md)
  (the error metric, fades, settling).

## Question

Each streamer does the same five things its own way, and on one page they share memory and bandwidth without
knowing of each other:

| What                | Particle tiles                                          | Particle pages (residency)                        | Streamed SOG (planned)                    | Volumes (density server)        |
| ------------------- | ------------------------------------------------------- | ------------------------------------------------- | ----------------------------------------- | ------------------------------- |
| reads               | `ParticlePackSource` (a file's slices, HTTP ranges, memory) | the same pack's ranges                         | files by URL relative to a manifest, a zip | a box by query                  |
| what is wanted      | `ParticleTileVisibility` on the CPU: frustum, size on screen, occlusion by opacity blocks | the cull's wants on the GPU, read back | a level a leaf by distance within a splat budget | a box about the camera or a selection |
| scheduling          | `ParticleTileLoader`: states, concurrency, aborts, retries, targets once for all tiles | its own reads, failures backed off | its own: two at a time, coarse and near first | one query at a time            |
| cache and release   | none of its own: "what to let go is the app's"          | a CPU cache a pack (`ParticlePageCache`: 256 MB, least recently used); its GPU pages let go 10 s after the last want | a cooldown of 100 updates            | the last box                    |
| GPU                 | the page store's pages                                  | the page store's pages, the fallback by copy      | texture array layers, a free list         | a volume's texture             |

meso\* adds Zarr pyramids (tomograms from the CryoET Data Portal), segmentations, meshes, and its package format
(a head and packs, read lazily by range requests from a zip, a directory or an object store); its runtime has a
streaming sketch behind a tile-source port. Each written as its own streamer is a sixth and seventh copy.

## Design

What is shared, by layer, and what stays a kind's:

1. **Sources** (`@molstar/core`, Mol\*-free so meso-io implements them): a `ByteSource` (`read(offset, length,
   signal)`, its byte length: `ParticlePackSource` generalized) and a `KeyedSource` (a store of keys to bytes and their
   ranges: a zip's entries, a directory, an object store, URLs relative to a manifest). Shared behavior: adjacent
   ranges coalesced, requests in flight limited by origin, retries (`retry-if`), aborts that drop a read nothing else
   waits for. The plugin's asset manager gives the browser's (a `File`'s slices, fetch with ranges); meso-io's package
   stores theirs. Why: network behavior written once, and the meso\* package format the one port every kind reads by.
   A zip's entry is read by range only where it is stored, not deflated (a pack's blocks and SOG's images are
   compressed already: storing them costs nothing); a deflated entry is read whole (Asked of meso\* 1).
2. **The unit graph** (a kind's adapter, model side): what can be streamed, declared up front and light. A unit has
   bounds (in the data's frame, with its reach), a level and its error in Å ([gl/mesoscale-lod.md](gl/mesoscale-lod.md)
   1), its bytes, its frames, its reads (keys and ranges of a source) and its dependencies (units shared by many: a
   particle target, a SOG file's codebooks and palette, a volume's header). Refinement is either `replace` (a finer
   unit replaces its parent's region: SOG levels, a Zarr pyramid, mesh levels, a template's coarse samples) or `add`
   (it adds to the coarser: particle tiles, nested strides). A graph, not a tree: a template's units are wanted by any
   of its copies, a target by every tile. Adapters: a `ParticleTileSet` (its directory), a SOG `lod-meta.json`, Zarr
   multiscales, a meso\* manifest's resources. A unit is keyed by its resource's content hash and its range (the data
   IR's key with a part: a unit cached survives a recompile, an undo and a restore, the adapter keeping no record of
   what it drew). The graph's types (units, wants, refinement, dependencies) in `@molstar/core` with the sources:
   meso-runtime imports nothing else of Mol\*, and its tile-source port is this graph (step 8); the adapters
   implementing it are model's and graphics'.
3. **Wants** (a kind's, through one interface): each consumer says which units it wants and at what priority, at
   each update. A shared traversal by camera covers the common case: frustum, error under `D` pixels, a priority by
   size on screen, hysteresis margins, prefetch along the camera's motion. A kind adds hooks (the tiles' occlusion by
   opacity blocks) or bypasses it (the residency's wants read back from the GPU cull: the GPU knows better). Why not
   one selector: what decides a want differs by kind; how it is served does not.
4. **The scheduler and cache** (one a page): unit states (none, loading, ready, failed) a key and frame; one priority
   queue across consumers, a request no longer wanted cancelled, dependencies loaded once and counted by their users;
   one CPU byte budget, released least recently wanted first after a cooldown; a trajectory's next frames prefetched.
   A read is one for every consumer wanting the same bytes (keyed by source and range: two representations of a
   template, a target wanted by tiles and by pages), not one a unit. The CPU budget holds what the GPU let go and what
   is not uploaded yet, not a second count of what the GPU holds (the residency's two tiers). The compile's loads
   (meso\*'s: whole resources, in the spec's order) go through the same sources and queue, ranked before the view's
   wants they block (Asked of meso\* 3).
   Generalized from `ParticleTileLoader`, which has all of it but the budget. Why one a page: a budget each lets the
   kinds take each other's memory; one queue puts a near particle tile before a far splat chunk.
5. **GPU residency** (graphics): one GPU byte budget across consumers; upload pacing by bytes a frame (the tiles'
   time-bounded steps, about 4 ms) shared, so loads of several kinds don't land in one frame; the fallback (a `replace`
   unit's nearest resident ancestor, an `add` unit's coarser levels) and the fade-in on arrival
   ([gl/mesoscale-lod.md](gl/mesoscale-lod.md) 4) by unit; outstanding loads in the settled state that a frame at rest
   waits for. A paged allocator (fixed-size slots in a large buffer or texture array, a free list) shared where the
   layouts allow: page store pages, SOG layers, volume bricks.
6. **What stays a kind's**: decoding (the pack's blocks, SOG's images, Zarr's codecs), the GPU layout (pages in
   storage buffers, images in texture arrays), identity (a particle's key, a splat's index), and its want producer.

## Why

- **One memory budget and one bandwidth a page**, and priority across kinds by what shows.
- **A new kind is an adapter** (graph, decode, GPU placement), not a streamer: meso\*'s tomograms, segmentations and
  meshes come one at a time.
- **meso\***: its package format is a `KeyedSource`, its runtime's tile-source port this graph. Data reads by
  reference stay the app's sources, as meso\* has it ("the runtime never fetches by itself").

## Risks

- **Over-generalizing what is tuned**: the tiles' occlusion and the residency's GPU wants took measurement; they stay
  want producers of their kinds, not folded into a generic selector.
- **Priority across kinds**: error in pixels times area on screen is comparable in principle, a heuristic in
  practice; a weight a consumer if needed.
- **Regressing what works**: the tile loader and residency are migrated behind their own tests and measured equal (the
  yeast tiled load), not rewritten.
- **Designed top-down**: a shared layer designed before its second user fits the first. So it is extracted from the
  tile loader and checked by SOG, the second consumer, before it is called general.

## Steps

1. **Inventory and baseline**: the table above checked against the code; the yeast cell tiled and the yeast nested
   2x2 measured as the baseline (time to the first frame and to settled, bytes read, requests, CPU and GPU peaks).
2. **Sources in core**: `ByteSource` and `KeyedSource` from `ParticlePackSource` and the plugin's URL range source
   (`getUrlRangeSource`); coalescing, limits by origin, retries, aborts; tests by a fake source (latency, failures,
   aborts, coalesced reads). The pack stream reads through them.
3. **The scheduler and cache**: generic keys, frames and dependencies, the CPU budget and release; `ParticleTileLoader`
   on it (its API and tests kept), the residency's page reads on it (its 256 MB cache and 10 s release the shared
   ones); measured against the baseline.
4. **SOG as the second consumer** ([gl/gaussian-splats-lod.md](gl/gaussian-splats-lod.md) step 5 built on it): its tree
   a unit graph (`replace`), its selection the shared traversal with its splat budget, its files a `KeyedSource`; what
   the shared layer got wrong for it fixed here.
5. **The shared traversal and the metric**: error in pixels ([gl/mesoscale-lod.md](gl/mesoscale-lod.md) step 2), the
   tiles' visibility as the traversal with its occlusion hook.
6. **GPU**: upload pacing across consumers, the GPU budget, fallback and fade by unit; the paged allocator where page
   store pages and SOG layers fit one (else each its own, the budget shared).
7. **A third kind**: a Zarr multiscale volume (a CryoET Data Portal tomogram), an adapter only: what proves a new kind
   is one.
8. **meso\***: the package store as a `KeyedSource` (meso-io), the tile-source port as the unit graph, described for
   the adapter.

## Asked of meso\*

What this plan needs meso\* to change or decide:

1. **Streamable resources stored, not deflated, in a `.mesox`**: the package writer stores a pack (`.mspk`), a SOG
   file's images and any resource read by range as zip entries without compression (they are compressed already); a
   deflated entry is read whole. A head and separate packs (a file system, S3) need nothing.
2. **The tile-source port is the core graph**: meso-runtime's streaming port takes the unit graph's types from
   `@molstar/core` (Design 2), so a streamed population's units, wants and budget are the same as Mol\*'s and the
   adapter implements them once.
3. **The load through the app's scheduler**: the sources the app installs are the shared ones (Design 1), so the
   compile's reads and the view's wants share one queue and one budget; each resource's content hash (the manifest's,
   else its bytes') given to the units that key by it (Design 2).

## Open

- **Where it lives**: the sources, the graph's types and the scheduler in `@molstar/core` (meso-runtime imports only
  core: its port the graph's types, Design 2), the adapters in model, GPU residency in graphics; open whether the
  scheduler's implementation, not only its interface, is core's.
- **A cache across sessions**: units keyed by content hash (Design 2, the manifest's hashes) kept in the browser's
  storage (Cache API, IndexedDB) for the next session, within a budget of its own.
- **The budget's owner**: the plugin (one a page) or the canvas (several canvases on a page share memory).
- **Priority across kinds**: one metric (error in pixels times area) or weights by consumer.
- **The paged allocator**: shared between storage-buffer pages and texture-array layers, or the budget alone shared.
- **The density server's volume streaming** (a box by query, the plugin's `volume-streaming` behavior): migrated as a
  kind of its own, or left as it is.
- **Zarr in scope here** or in a volumes plan of its own, with this as its streamer.
