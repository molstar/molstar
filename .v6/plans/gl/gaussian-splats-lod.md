# Plan: Streamed SOG (levels of detail)

## Status

- **Current**: planned, nothing built. Decisions open (see Open).
- Related: [gaussian-splats.md](gaussian-splats.md) (step 10: SOG read in place, the `sog` encoding this extends; the
  GPU sort, its centers kernel), [storage-buffers.md](storage-buffers.md) (data pages),
  [particle-residency.md](../particle-residency.md) (pages of a store resident by the view, a budget),
  [streaming.md](../streaming.md) (proposed: step 5's loader, residency and cooldown on a streamer shared with the
  particle tiles; this plan its second consumer).

## Question

Draw splat scenes too large to hold (millions to tens of millions of splats) as PlayCanvas's streamed SOG: a
`lod-meta.json` over chunk files of SOG, each region of space at several levels of detail, the level by the distance
to the camera within a splat budget, the files loaded and released as the view moves.

## The format

What `splat-transform` writes (`src/lib/writers/write-lod.ts`), as read from a file of its v3.10.3 (AlterEyes' carder
bee: 7.0 M splats in 3 levels of 4.0, 2.0 and 1.0 M, 12 files, 138 MB):

- `lod-meta.json`: `version` 1; `asset` (`chunkGaussians`, `chunkExtent`, `chunkMinGaussians`: the writer's split
  limits); `count`, `counts` a level, `lodLevels`; `filenames` (`<level>_<i>/meta.json`); `tree`: nodes of a `bound`
  (min, max: the union of the children's), inner ones `children`, leaves `lods`: per level `{ file, offset, count }`, a
  range of a file. The carder bee: 16 leaves at depth 4, every leaf with every level (233 K to 270 K splats at level 0).
  Optional `environment` (a SOG of the background, `env/meta.json`), `extras` and `comments` (version 3 units only).
- The files: SOG as a folder (`meta.json` and its WebP images, not zipped), each with its own codebooks and palette
  (the carder bee's: version 2, palettes of 65 536 entries of degree 3); a file holds the ranges of many leaves of one
  level, packed until the writer's unit size (the carder bee's 0.58 to 0.75 M splats, images up to 868 × 864).
- The levels are what the source has (a streamed SOG, LCC: trained or simplified each), not made by the writer: a
  level a full representation of a leaf's region, not nested in the finer ones. The frame the 3DGS PLY's (y-down), the
  tree's bounds in it.

## PlayCanvas

The engine (`gsplat-octree.resource.js`, `gsplat-octree.js`, `gsplat-lod-tree.js`, `gsplat-budget-balancer.js`,
`gsplat-world.js`):

- A leaf's level by its distance `d` (to the nearest point of its bounds, the camera's FOV compensated, behind it
  penalized ×1.5): `1 + log_m(d / (s b))`, floored and clamped; `b` 5 and `m` 3 by default, `s` one scale for the scene,
  the largest whose total fits the splat budget (1 M by default). Re-evaluated when the camera moved a unit, turned
  90° or changed its FOV by 2%, or a load completed. A region keeps its loaded coarser level until the finer loads (and
  the finer while the coarser does).
- Files loaded whole (2 at a time, coarse levels first, then the near), referenced by the leaves selected, released a
  cooldown (100 updates) after their last reference.
- All selected ranges in one world: an incremental work buffer (the changed blocks rendered again, a copy into its own
  layout), one sort. No fade between levels.

## Proposal

Keep reading SOG in place (no copy into another layout): the loaded files' images as layers of texture arrays, the
selected ranges' splats drawn and sorted as one object.

- **Sources**: a `lod-meta.json` by URL, its files fetched relative to it on demand (the plugin's asset manager); a
  zip of the folder (all of it in memory: the levels and the budget still apply to the GPU's memory). A folder of
  loose files on disk later. The single SOG as a folder (a `meta.json` and its images) the same reader: the archive's
  entries by name, from a zip or by fetching (`parseSogArchive` given an accessor).
- **Residency**: the images of each loaded file in a layer of a `texture_2d_array` per image kind (means low and high,
  scales, rotations, base color, labels; the palettes in arrays of their own), layers of the largest file's size (1024 ×
  1024 holds a unit of 1 M), allocated in blocks, released to a free list. Per layer its file's codebooks (a row of a
  `256 × layers` float texture), means ranges, image width, palette layout (a small layer table). WebGPU uploads a
  bitmap into a layer (`copyExternalImageToTexture`, its origin's z), WebGL2 by `texSubImage3D`. WebGL1 has no arrays:
  the coarsest level alone, as one `sog` object a file (or not supported: open).
- **Selection** (as PlayCanvas, the same parameters): per leaf a level by distance within the budget; the files the
  selection needs loaded, coarse first, near first, a few at a time; a region drawn at its loaded level nearest the
  wanted. Evaluated in `Renderer.prepareGaussianSplats` (outside the passes, where the sort is), on the view's change by
  the thresholds above or a load's completion. Loads asynchronous through a loader the geometry holds (the plugin's
  fetch and decode), their completion requesting a frame.
- **Drawing**: one render object, the `sog` encoding with layers (`dSplatEncoding` `sogLayers`): the draw's splats a
  list of `(layer, index)` (the selected ranges', `tSplatOrder` and `sSplatOrder` over them), `uVertexCount` the selected
  count; the shader's decode as step 10's, its texel in the splat's layer, the layer's codebook and ranges. The sort over
  the selected (the GPU's centers kernel per layer; WebGL2's CPU sort from each file's centers, read back as step 10 at
  its load).
- **Groups**: a group a splat of the whole scene (picking, loci, overpaint by group as today) or a leaf (cheaper: 16 for
  the carder bee, a region picked and themed as one): open.
- **Bounds**: the tree's (exact, known before any file), a leaf's loci bounds its bound; upright as the shape's
  transform (step 10).
- **Memory**: the budget's splats resident (1 M at about 30 bytes, plus each file's palette: up to 4 MB), the files'
  layers released a cooldown after their last selection.

## Steps

1. **The folder SOG**: `parseSogArchive` over an accessor of entries by name (a zip's, a fetch's), the SOG provider
   for a `meta.json` with its images; tests: an archive by accessor, the zip unchanged.
2. **The manifest**: `io/src/reader/sog/lod.ts` (`lod-meta.json` parsed: tree, leaves, files), `model` the leaves'
   bounds and ranges by level, the selection as a function (camera position, FOV, budget → a level a leaf) with its
   tests (bands, the budget's scale, a child never nearer than its parent).
3. **Layers**: texture arrays on both backends (`TextureSpec` of a `2d-array` dimension: the WebGL2 declarations, WGSL
   `texture_2d_array`, uploads into a layer from a bitmap), with tests on WebGPU (Dawn: a layer written, read back).
4. **The `sogLayers` encoding**: values, the layer table, the shaders' decode by layer, the draw's list; the GPU sort's
   centers per layer; a scene of two files of synthetic SOG drawn as their single `sog` objects (headless WebGPU).
5. **Streaming**: the loader (URL relative to the manifest, a zip), residency (layers, free list, cooldown), the
   selection in `prepareGaussianSplats`, loads requested and completed into frames; the provider `lod-meta.json`.
6. **Checks**: the carder bee by URL on WebGL2 and WebGPU (moving through it: the levels, the budget kept, no holes),
   time to first frame and GPU bytes; a visual case of a synthetic manifest (levels by distance).

## Open

- Groups: a splat each, or a leaf each.
- WebGL1: the coarsest level as plain `sog` objects, or not supported.
- The environment SOG (a background of splats): drawn as an object of its own, or later.
- A fade between levels (PlayCanvas has none): later if popping shows.
- Sources: a folder of loose files on disk (the browser's directory picker) after the URL and the zip.
