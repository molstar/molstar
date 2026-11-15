# Formats

`src/formats/`: readers into models, coordinates, topologies, volumes and particle lists; writers.

## 1. Structures

- `ModelFormat<T> = { kind, name, data }`; each format has a namespace with `is()` (by `kind`). PDB models are
  `MmcifFormat` with the PDB kept in `data.source`.
- Every structure reader builds `BasicData` (`basic/schema.ts`, a subset of mmCIF tables; `createBasic` fills missing
  tables and, normalizing, substitutes `label_*`/`auth_*` columns for each other) and calls
  `createModels(basic, format, ctx)` (`basic/parser.ts`), which returns an `ArrayTrajectory`:
  - standard: `atom_site` split by `pdbx_PDB_model_num`, each run sorted (`sortAtomSite`) and made a model; a model
    whose tables equal the previous one's shares its topology (model.md §3);
  - integrative (`ihm_model_list` present): atomic plus coarse spheres and gaussians.
- Readers: mmCIF/BinaryCIF (`mmcif.ts`; `MmcifFormat` sits apart in `mmcif-format.ts` so code can recognize mmCIF models
  without loading the parser), CCD (ideal and model coordinates), PDB/PQR/PDBQT (via `pdbToMmCif`), GRO, MOL, SDF, MOL2,
  XYZ, CIF core, cube, LAMMPS data and trajectory.
- Coordinates only: DCD, XTC, TRR, NetCDF, LAMMPS trajectory. Topologies only: PSF, PRMTOP, TOP (bonds as 0-based index
  pairs).

## 2. Format properties

`FormatPropertyProvider.create(descriptor, { asDynamic?, cachePerFormat? })` (`structure/common/property.ts`): a
registry from `ModelFormat.kind` to an `obtain(model)`.

- `get` stores the result, `undefined` included, in the static or dynamic bag under the descriptor's name and checks
  that bag, not `customProperties`: same-topology models share `customProperties` but not dynamic data, so the
  descriptor's presence would hide the value of every model after the first one asked.
- `cachePerFormat` memoizes per `model.sourceData` (a `WeakMap`): one value for every model of a file.
- Registration is an import side effect of the reader module (`mmcif.ts`, `cif-core.ts`; `pdb.ts` imports `mmcif.ts` for
  it).

| Provider                  | Bag                 | From                                                                         |
| ------------------------- | ------------------- | ---------------------------------------------------------------------------- |
| `ModelSymmetry`           | dynamic, per format | mmCIF (assemblies, spacegroup, NCS); `fromCell` by other readers, frames     |
| `Model.PeriodicBox`       | dynamic             | frames with a cell or `Frame.box`; GRO; LAMMPS; NAMD `.xsc` (`namd-xsc.ts`)  |
| `ModelSecondaryStructure` | static              | mmCIF `struct_conf`, `struct_sheet_range`                                    |
| `StructConn`              | static              | mmCIF `struct_conn`                                                          |
| `ComponentBond`           | static, per format  | mmCIF `chem_comp_bond`; set by CCD                                           |
| `IndexPairBonds`          | dynamic             | mmCIF `molstar_bond_site`; set by MOL/SDF/MOL2, CIF core, LAMMPS, topologies |
| `AtomSiteAnisotrop`       | dynamic             | mmCIF (by each model's atom ids), CIF core (by atom label)                   |
| `AtomPartialCharge`       | static              | set by PDBQT, MOL2, LAMMPS                                                   |

`IndexPairBonds` index atoms by source index (`atomSourceIndex`); `maxDistance` (default -1: an element-based threshold;
`Infinity`: trust the table; a per-bond `distance` takes precedence), `cacheable`, and operators per bond (assemblies
written with operator ids). The frames of `trajectoryFromModelAndCoordinates` share the model's value: no format
provider gives a value a reader set.

## 3. Volumes and particles

- Volumes: CCP4/MRC, DSN6, DX, cube, density server, segmentation, structure factors and MTZ (through
  `shared/reflections.ts`). Each gives a `Volume` with one identity instance. CCP4 header stats are trusted unless they
  are all zero or NaN.
- Particles: Artiatomi, Dynamo, RELION STAR, CryoET Data Portal NDJSON, Simularium (a lazy `ParticleTrajectory`), mmCIF
  assemblies (`createParticleListFromMmcifAssembly`), particle packs (`pack.ts`).
- A particle pack holds what a `ParticleList` and its `sites` targets use, nothing else: positions, and rotations and
  radii where the list has them; fibers, numeric attributes, the cell; a row a target (its id, name, entity,
  compartment, particle range, radius, its sites and their chains as runs; the sites quantized: int16 positions by the
  target's origin and scale, 8 bytes a site with a type indexing a table of radii and, where every target's sites have
  them, elements (atomic numbers, version 5; a type a distinct radius and element: S and P share a radius); the
  rotations snorm16, 8 bytes a particle within about 1e-4 rad, read as they are: `ParticleRotations`; the positions on a
  21-bit grid of the list's box, 8 bytes a particle within half a step, the yeast cell's 0.011 Å: `GridPositions`). Why:
  it is the list's storage, read as views, not a document for other tools (mmCIF stays the exchange format).
  - The particles sorted by target, a target's particles a range: its id, entity and compartment stored once, not per
    particle; no keys (a particle's key is its index in the pack). The writer refuses a target whose particles differ in
    what it stores once. A target's particles' radius is one by definition: stored once (the largest where a list's
    differ: an operator's scale read from its quantized matrix, 0.15% in MG_1184).
  - A 64-byte header (`MSPK`, version, byte order mark, the JSON's offset and length), the blocks at multiples of 256
    bytes (WebGPU's copy alignment), then the JSON of the small tables and the blocks' directory: raw little-endian
    arrays, no encoding, so they are views of the file (a compressed block cannot be viewed; transport compression is
    the server's). The JSON after the blocks: their offsets known when it is written.
  - The reader's list holds the targets, entities, compartments and radii as runs (`ParticleValues`), a value a target,
    and no keys.
  - A target's bounds (`bounds`, optional in its row: its positions' box and largest distance from its center, rounded
    outward to 0.001 Å), what its placement needs of its positions without a pass over them (`Particle.getRangeBounds`:
    the grid's cell size); a pack without them (older, the same version) has them computed at the first use, one pass
    for all targets. Not a format version: an optional field, read where present.
  - Tile-ordered (version 6, the writer's `tiles`, opt-in): the particles by tile of a grid of n³ cells over their box
    (n a power of two, as `ParticleTiles` makes it), the tiles in Morton order of their cells, a tile's particles by
    target; a target's particles runs, one a tile it is in (`target.runs`: first particle, target id), its row's
    `particleOffset` -1; the tiles a directory (`tile.ranges`: cell id, first particle, count; `tile.boxes`: the
    centers' box, float32 rounded outward) given to the list (`ParticleTiles.setDirectory`). A tile set of that grid or
    a power of two coarser (a coarse cell's cells consecutive in Morton order) is made from it, no pass over the
    particles, a tile a range of every block; another grid orders the list again. Why the written tile size, not a finer
    grid: a fine cell's targets are runs of their own, and the yeast's runs at 16³ are 2.29 M, at 32³ 5.9 M, at 64³ 11.6
    M (a run a store record). No bounds a run (2.29 M rows): computed where asked. The keys are the pack's indices, as
    any pack's. Its list drawn whole: where culled on the GPU by its tiles (`ParticleTileView.ofWhole`: a tile a range,
    the page store's object), the yeast's 3.0 s against 2.2 s by target; else a target an index set, not a range (its
    poses gathered, its bounds its runs' united, `Particle.getIndicesBounds`), 5.3 s.
  - A tile-ordered pack's tiles' opacities (`tile.opacity`, `tiles.opacityBlocks`: 4 by default, 0 none): a tile's
    blocks' optical densities (`ParticleTileOpaque.computeOpacity` of its particles as written, its targets' volumes;
    256 B a tile at 4³), what the camera's occlusion needs without a compute pass over the particles (the yeast's
    written in 2.8 s more); a coarser set's tile its cells' blocks side by side ((2b)³: a molecule spread over its own
    cell's blocks, so not the coarse tile's own computation block by block, the same cross-sections).
  - A pack by target's density grid (`density`, `densityBlocks`: 64 along the box's longest side by default, 0 none; a
    tile-ordered pack none, its tiles' opacities are): its particles' optical densities on near-cubic blocks over their
    centers' box (`ParticleTileOpaque.computeDensity`, as `computeOpacity`; about 1 MB, the yeast's written in 1.3 s
    more), read onto the list (`ParticleDensityGrid`): what the GPU's occlusion is primed from as the list is drawn
    whole, without a pass over the particles at load (computed there, opt-in: `ParticleDensityGrid.computeTask`, in
    chunks, the yeast's 2.6 s).
  - A tile-ordered pack read by range (`pack-stream.ts`): a source of byte ranges (`ParticlePackSource`: bytes, a `Blob`
    by `slice`, a URL by HTTP ranges in the plugin; `batched`: a few reads in flight, ranges near each other one read),
    its header and document, then the tile set (`createParticleTileSetFromPack`, on its written grid): the directory,
    the opacities and the targets' runs read at once (the yeast's 22 MB of 731 in 0.16 s), a tile's range of the
    position, rotation and attribute blocks when its tile is asked, a target's sites when it is (once). A tile's list as
    a whole list's tile is (its keys the pack's indices, runs, its bounds by its box); the catalog from the targets'
    rows (their reach as written: required). Its fibers not read (a tile's list has none, as a whole list's tiles). Why:
    a pack larger than an `ArrayBuffer`, and the first tiles drawn before the rest is read. A target's sites every
    `siteStride`-th where asked (`ParticleSites.coarsen`: read, the kept ones copied; drawn as that stride's level of
    detail, `ParticleSites.stride`): a view's targets a tenth of the memory at 10, its floor on a phone.
  - A URL read by range needs its server to answer `Range` requests (206) and to tell the file's size: the answer's
    `Content-Range`, which another origin reads only where the server sends
    `Access-Control-Expose-Headers: Content-Range` (and allows the `Range` request header), else a `HEAD`'s
    `Content-Length` (CORS-safelisted). A server answering neither: the pack read whole, the plugin's log saying why
    (`getUrlRangeSource`).
  - A block's sizes agree (`parseParticlePackDocument`): its `byteLength` its `count` × `components` values of its type,
    a block of a value a particle (`position`, `rotation`, `attribute.*`) of the pack's `count`, else it is not read,
    named. Why: a view of it reads `count` × `components`, into the next block where they disagree.

  - Coarse levels of a tile-ordered pack (`tiles.levels`, `WriteParticlePackOptions.levels`, the CLI's `--levels`): an
    octree over the tiles' grid, level `L` (1 to `log2(dimensions)`) a node a cell of the grid of `dimensions >> L`
    cells an axis, every `8^L`-th of its tiles' particles (in Morton order of a 512³ grid over the box: spread over it),
    copies in blocks of their own (`node.position`, `node.rotation`, `node.attribute.<name>`; `node.keys` a copy's index
    in the pack, `node.runs` its targets' runs, `node.ranges` a node's level, cell, first and count, `node.boxes` its
    centers' box). Read by range a node is a tile of the set (`ParticleTileLevels`: its slots after the particles', its
    particles drawn `2^L` larger, `ParticleList.scale`); a reader of none, and the whole read, ignore them. Why: a far
    region drawn as a few thousand particles read as one range (the yeast's root 7.7 k of 31.5 M), a whole cell within a
    phone's memory; copies, not the tiles' prefixes (a node's sample of 512 tiles one read, not 512). About a seventh
    more particles written (the yeast's 731 -> 825 MB, 17 s).
  - A target's reach and volume (`reach`, `volume`, version 6: `ParticleTileCatalog.getTargetRadius` and
    `getTargetVolume` of its sites as written, what a reader computes): a tile set's catalog without a pass over the
    sites (the yeast's 0.8 s); computed where absent.
  - A target's sites' extremes (`extremes`, optional in its `sites`: the quantized positions' min and max of each chain,
    of all sites without chains; `ParticleSitesQuantized.extremes`): a chain's box, its boundary's, decoded as the sites
    are with no pass over them (the yeast cell's 20 M sites: 113 ms at load); computed where absent. An optional field
    as the bounds are.
- A list's per-particle values (`ParticleList.targets`, `entities`, `compartments`, `radii`) are `ParticleValues`: an
  array, or runs (a value a run of consecutive particles). Why: a list sorted by target holds a target's values once,
  not a particle (the yeast cell: 0.63 GB). Readers take either: `ParticleValues.get` for a particle (a binary search,
  the run found last tried first), `forEachRun` (`forEachRun2` for two) to visit all, `map` to derive per-particle
  values in the same form (a theme's color index: a value a run), `toArray` only where a whole array is needed. `keys`
  is optional (none: a particle's key is its index).

## 4. Writers and export

- `writer/`: MOL, SDF and MOL2 encoders over `atom_site` (component atoms and bonds supplied by the caller), and
  `CCP4Writer.writeMrc`.
- `model/structure/export/mmcif.ts`: `to_mmCIF`, `encode_mmCIF_categories` for structures of one model. Categories are
  copied from the source or generated (`_entity`, `_atom_site`, secondary structure, `_struct_asym`, `_entity_poly*`,
  `_chem_comp*`, `_pdbx_nonpoly_scheme`), plus custom properties with `cifExport`. Mol\* extensions:
  `molstar_atom_site_operator_mapping` (read back by the parser) and, if asked, `molstar_bond_site`.
