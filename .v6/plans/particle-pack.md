# Plan: a particle pack file (mmCIF models, transforms and sites as views)

## Status

- **In progress**: steps 8 and 9 done (the rotations snorm16, the positions on a 21-bit grid: the yeast cell's pack
  1,089 -> 838 -> 712 MB); step 6 done (the impostors drawing the pack's quantized sites as they are, sized by their types: the
  yeast cell's scene 1.80 -> 1.50 GB of JS heap, its GPU data 531 -> 376 MB); step 7's lists by use, the GPU culler's scene cells and the targets' coarse grids done (the yeast
  cell's buffers 10.68 -> 7.56 GB by heap snapshot); step 5 on the CPU and WebGPU (`instance-poses.md`: the buffers
  7.56 -> 4.75 GB, the heap 6.31 GB); steps 1 to 4 (the format's reader and writer, `formats/particles/pack.ts`; the converter, `cli/particle-pack`;
  the plugin's `particle_pack`; the list by runs: the yeast cell's heap 7.29 -> 3.17 GB); the yeast cell's memory
  inventoried (below). Builds on
  [mesoscale-scale.md](mesoscale-scale.md) (the yeast milestone: its table of what stays and its step 4, the file's
  bytes) and [mesoscale-load.md](mesoscale-load.md).

## Question

`yeast_cell_atoms.bcif` (`.benchmark/`, 1.17 GB) is mostly transforms, stored as mmCIF assembly operators:

| category                    | rows   | bytes  | what                                                                      |
| --------------------------- | ------ | ------ | ------------------------------------------------------------------------- |
| `_pdbx_struct_oper_list`    | 31.5 M | 944 MB | 12 columns: the rotation FixedPoint 1e-4 as int16, the vector 1e-2 as int32 |
| `_atom_site`                | 20.3 M | 206 MB | the 5,639 targets' reference atoms                                        |
| `_pdbx_struct_assembly_gen` | 857 k  | 16 MB  | `oper_expression` as id ranges (`(1-976)(10018299-10018332)…`)            |
| `_entity` and the rest      | 5.6 k  | 0.3 MB |                                                                           |

One operator a particle, none shared: the assembly's indirection (gen rows, expressions, operator ids) costs a parse and
buys nothing here. Reading it copies every transform several times before the GPU sees it:

1. the downloaded bytes, kept by the asset cache;
2. `decodeMsgPack`'s copy of every binary column (on purpose, `msgpack/decode.ts`);
3. the FixedPoint columns decoded to floats (`readParticleOperators`, 12 columns);
4. the list's loop: the expressions parsed, each operator's matrix turned into a quaternion, the list's
   positions and rotations (28 B) and its per-particle indices (keys, targets, entities, compartments: 16 B) written
   (9.2 s for the yeast cell);
5. the targets' transforms (`writeParticleTransform`, 64 B an instance), then the GPU's rows.

`mesoscale-scale.md` counts 1–3 at 2.2 GB and the list at 1.75 GB.

## The format

What a `ParticleList` and its `sites` targets use, nothing else (`model/_spec/formats.md` §3, `formats/particles/pack.ts`).
One file, little-endian:

```
offset 0   header, 64 bytes: 'MSPK', version u32 5, byte order mark u32 0x01020304, the JSON's offset u64 and length u64
256 ...    the blocks, each at a multiple of 256: raw arrays, no encoding
after      the JSON: the list's label, entry id and cell, its entities, compartments and attributes, a row a target,
           the blocks' directory (name, type, components, count, offset, length)
```

- **A row a target**: `id` (the list's target id), `name`, `entity`, `compartment` (indices, -1 for none),
  `particleOffset`, `particleCount`, `radius` (where its particles' are the same) and `sites` (their range of the
  site blocks, label, center: the point a particle's position places, `ParticleTarget.getCenter`; the chains as runs:
  labels and counts; each chain's quantized extremes, `extremes`, its box with no pass). The particles are sorted by target, so a target's are a range: targets, entities, compartments
  keys and radii are not stored per particle (20 B a particle, 630 MB for the yeast cell; a key is a particle's index,
  a radius one a target by definition, the largest where a list's differ).
- **The blocks**: `position` on a 21-bit grid, two u32s (step 9, the grid in the document's `positionGrid`; float32
  in versions 2 and 3, still read) and `rotation` (a unit quaternion) snorm16 × 4 (step 8; float32 in version 2, still
  read), exactly `ParticleList.coordinates` and `rotations` (16 B a particle, 504 MB for the yeast cell, against 944 MB
  of operators);
  `site` int16 × 4 (8 B a site: its position quantized by its target's origin and scale, its type an index of the
  document's `siteRadii` and, version 5, `siteElements` (atomic numbers, where every target's sites have them: a type
  a distinct radius and element); step 6); `fiber.offsets`, `fiber.indices` (int32, the
  indices of the particles at their packed positions); `attribute.<name>` float32 (numeric attributes).
  - Precision: the yeast cell spans ±23,700 Å, where float32 keeps about 0.002 Å, the grid 0.011 Å (half its step);
    the operators' FixedPoint 0.01 Å.
  - Not 12-float affine rows (the GPU's `aTransform`): 48 B, 1.51 GB for the yeast cell. The GPU builds the rows from
    the poses (`gl/instance-poses.md`).
- **Not mmCIF**: no `_atom_site`, entities or assemblies; the pack is the list's storage, not an exchange document.
  The converter writes standard mmCIF back from a pack (operators from the poses) for tools that need it.
- **Aligned to 256 bytes**: typed-array views need 4 (float32), and 256 is WebGPU's copy-offset alignment
  (`writeBuffer` from the file's bytes at a block's offset).
- **No compression of the blocks**: a compressed block cannot be viewed. Transport compression (gzip, brotli) is the
  server's.

## The sites block (step 6)

The pack's sites were `ParticleSites` as they are (float32 positions and radii, 16 B a site), copied into a `Spheres`
geometry (center and group, 16 B a sphere) and a size texture: about 36 B a site, 730 MB for the yeast cell's 20.3 M
sites. A block in their own layout, read as views and drawn from (the format and the CPU's reads done, the drawing
open):

- **`site`**, a site 8 B: `int16 × 3` its position relative to its target's origin, times the target's scale, and a
  `uint16` type, an index into the pack's radius table. 162 MB for the yeast cell, against 325 MB of float32 sites. It
  is a vertex format as it is (WebGPU `sint16x4`, WebGL a `SHORT` × 4 attribute).
- **Sorted** by target, then chain: a target's sites are a range, its chains runs in it.
- **A target's row** gains the sites' origin and scale (from the
  target's extent: 0.01 Å keeps ±327 Å, a 3,000 Å fiber 0.09 Å, far below a van der Waals radius), the chains' runs and
  the target's bounding sphere (written by the converter, so bounds and levels of detail read no positions).
- **The site types** (the document): a radius each.
- **Drawn** from the block: the sphere impostors' quantized centers (GLSL and WGSL: `vec3(a.xyz) * uScale + uOrigin`,
  the group the sphere's index, a site's groups being 0…n-1), a target's range uploaded from the file's bytes, so no
  `Spheres` copy and 8 instead of 16 B a sphere on the GPU. The size texture stays (4 B a site, from the radius
  table): the size themes keep working.
- **Decoded where read**: the blob surface's density, the sphere-mesh fallback, a site's label or picked position (a
  transient copy, or one site). `ParticleSites` gains the quantized form (the data, origin, scale, types and radius
  table); `positions` and `radii` become decoded on demand.

## Reading

- **Browser**: one `fetch`, one `ArrayBuffer`. The header, then the document (JSON), then views:
  `new Float32Array(buffer, block.byteOffset, block.count * block.components)`. The views keep the bytes alive:
  the file is the particles' storage, not a copy of them. For a file past the `ArrayBuffer` limits (or one not wanted
  whole), a `Blob` or `File` read by `slice` per block, a target's range only when it is drawn; HTTP range requests
  likewise.
- **Node**: one `Buffer` from the file (its own `ArrayBuffer` at offset 0 for a large read, so the views stay aligned;
  checked, else copied once). A true memory map (Bun, Deno, a native module) changes nothing in the format.
- **Big-endian hosts**: refused (none in use; the header's byte order mark checked).
- **The list**: `ParticleList.coordinates`, `rotations` and the targets' sites as views; the targets from their rows
  (ranges); no loop over the particles but the per-particle `keys`, `targets`, `entities` and `compartments` the
  list requires today, filled from the ranges until they are derivable (step 4).

## Writing

- From a `ParticleList` and its `sites` targets: the pack is a `ParticleList` serialized, so any particle format
  (cellPACK, PetWorld, STAR, Simularium) can be packed, not only mmCIF assemblies; other target kinds (structures,
  shapes, volumes) are not packed.
- A converter (`cli/particle-pack`): BinaryCIF (cellPACK, PetWorld) to a pack through
  `createParticleListFromMmcifAssembly`, and a pack back to standard mmCIF (the operators and gen rows written from the
  ranges).

## What the renderer holds, the yeast cell

Measured after step 5's GPU half and step 7 (the instances placed from poses, `instance-poses.md`; the lists by use,
the scene cells): the Viewer, Chrome, WebGPU; the first frames at 33.8 s, the renderer 6.67 GB, all processes
11.06 GB after the first frames, the GPU culler's instance data 1.39 GB (44 B an instance). A heap snapshot after the
first frames (847 MB, an unminified build for the owners' names), each array buffer followed up its retainers to an
owner. The heap is 6.31 GB: array buffers of at least 1 MB 4.75 GB (100), smaller ones about 1.2 GB (mostly per
target), JS objects 0.36 GB. By owner (MB as the snapshot gives them, MiB), and which step takes it:

| owner                                                                                     | MB    | taken by       |
| ----------------------------------------------------------------------------------------- | ----- | -------------- |
| the downloaded bytes (the asset manager's `HeldAsset.bytes`)                              | 1,112 | the pack is it |
| the parsed BinaryCIF's columns: `pdbx_struct_oper_list` 900, `atom_site` 195              | 1,095 | 1-3, 6         |
| the list: `rotations` 480, `coordinates` 360 (the poses read from them)                   | 840   | 1-3 (views)    |
| the list: `keys`, `targets`, `entities`, `compartments`, `radii` (120 each)            | 600   | 4              |
| the particles' label indices (`getParticleLabel`: entries 120, combinations 140)          | 260   | 4 (by range)   |
| the color index                                                                           | 120   | 4 (by range)   |
| the merged targets' poses: `indices` 120, `order` (the grids') 49 at least                | 169   | 4 (ranges)     |
| `aInstance`                                                                               | 125   | 5 (implicit)   |
| the merged values' parts (`SegmentedArray`), the culler's indices, the fibers' indices    | 178   | stays, 7       |
| the sites (the current target's), the spheres' `tPositionGroup`, `tSize`, `tColor`        | ~75   | 6              |
| the font atlas, a value presented (78) and the rest at least 1 MB                         | ~180  | stays          |
| JS objects (5,639 targets' representations, values, renderables)                           | 359   | stays          |

Gone since the first inventory (10.76 GB): the instances' transforms (1.92 GB: `aTransform`, `transform`,
`cellTransform`), the multi-draw lists (2.15 GB), the per-target grids' cell spheres and batches (0.55 GB). The
largest owners left are the format's: the download and the parsed columns (2.2 GB), the list's arrays and indices
(1.8 GB): steps 1-4.

## Expected, the yeast cell

From the inventory (the renderer's heap; the process adds Chrome's own, about 0.45 GB here), to be measured at each
step (Node and Chrome, the yeast cell and the isg grids, what stays and the peak):

| what                 | today    | steps 1-3 | + 4      | + 5 (GPU) | + 6 (sites) | + 7 (scene grid) |
| -------------------- | -------- | --------- | -------- | --------- | ----------- | ---------------- |
| download             | 1.17 GB  | ~1.09 GB  | same     | 0.84-1.09 GB | 1.0-1.25 GB | same          |
| first frame          | ~43 s    | ~34-35 s  | same     | ~27-31 s  | ~22-27 s    | ~22-27 s         |
| renderer heap        | 10.76 GB | ~9.0 GB   | ~7.8 GB  | ~5.7 GB   | ~4.9 GB     | ~2.3 GB          |
| GPU instance data    | 1.76 GB  | same      | same     | ~1.0-1.1 GB | same, spheres -160 MB | same   |

- Steps 1-3: the operators' MessagePack copy (0.90 GB) and the list's positions and rotations (0.84 GB) gone, the
  pack (1.09 GB) in place of the download (1.11 GB); the list's build 9.2 s to milliseconds.
- Step 4: the list's per-particle arrays, the label indices, the target sets and the color index by range (1.16 GB).
- Step 5: the instances' transforms and `aInstance` (2.12 GB); the culler's 56 B an instance to ~32-36 B on the GPU, the
  uploads 3.3 to 1.4 GB, the transforms' writes (est. 3-8 s); `snorm16` rotations shrink the download (0.84 GB).
- Step 6: `atom_site`'s copy, the sites and the spheres' values against a 162 MB block (0.57 GB here, with the block
  in the file), the decode, grouping and spheres' build (est. 3-6 s).
- Step 7: the scene-wide grid: the multi-draw lists (2.15 GB) mostly, the per-target grids and the culler's cell
  table (about 0.9 GB with what of the unattributed it is); fewer CPU culling tests a frame.
- Not changed: the values that stay, the JS objects, Chrome's GPU process (about 5 GB here
  after the first frames: all processes 16.4 GB, the renderer 11.2).

Measured with step 5 (CPU and WebGPU) and 7 done, before steps 1-4 and 6: the first frame 33.8 s, the renderer heap
6.31 GB (the renderer 6.67 GB, all processes 11.06 GB), the GPU culler's instance data 1.39 GB (44 B an instance: a
pose a particle, 32 B, not shared across targets).

## Steps

1. **The spec and its reader and writer** in `@molstar/model` (`formats/particles/pack.ts`): the header, the
   document, the blocks. Done: a round trip in the quick lane (generated lists, an mmCIF assembly); PetWorld and
   MG_1184 converted, their lists equal to the BinaryCIF's, every field exact. MG_1184 (2.12 M particles, 848 targets,
   2.94 M sites): 53.4 MB of BinaryCIF, a 112.3 MB pack (the sites 47 MB, positions and rotations 60 MB; a radius a target,
   the largest of its particles', which differ by up to 0.15% by the scale of their operators' matrices, quantized to
   1e-4), the list built in 1,225 ms from the
   BinaryCIF, 9 ms from the pack. PetWorld (0.52 M particles, 3.68 M sites): 88.2 MB, a 73.4 MB pack, 554 ms and
   11 ms.
2. **The converter** (`cli/particle-pack`, both ways). Done: mmCIF (`.cif`, `.bcif`) to a pack through
   `createParticleListFromMmcifAssembly`; a pack to mmCIF in cellPACK's layout where every target is one chain (an
   entity a target, its description its compartment and name, a fiber a gen row), PetWorld's otherwise (a model a
   target), a site an atom of an element of its radius, a particle an operator (`R`, `p - R c`), read back as the same
   list (MG_1184, PetWorld: positions within 0.002 Å, rotations, sites, fibers exact). The packs in `.benchmark/`: the
   yeast cell 1.25 GB (1.17 GB of BinaryCIF), the isg grids (4×4×4 842 MB, 894 MB), MG_1184, PetWorld. In Node, the
   list from the file: the yeast cell 10.3 s -> 0.1 s, 6.99 -> 1.96 GB (the pack 1.25 GB and the per-particle arrays
   filled from the ranges, 0.63 GB, until step 4); the 64 isg cells 5.6 s -> 0.1 s, 6.52 -> 1.50 GB.
3. **The plugin**: the format registered (sniffed by its magic, extension `.mspk`), loaded in the Viewer and Prism as
   `mmcif_particles` is. Done: `particle_pack` (`plugin/state/formats/particles/pack.ts`, the particles registry, Prism's
   particles feature), a transform from the data to the list, its arrays views of the downloaded bytes. The yeast cell
   in the Viewer (WebGPU), the same run's BinaryCIF beside it: the list 9.5 s -> 0.17 s, the first frames 42.6 -> 34.9 s,
   the heap 7.29 -> 3.77 GB at its peak, the renderer after the first frames 6.09 -> 4.18 GB, all processes 10.42 ->
   7.44 GB; the image the same.
4. **The list by ranges**: `ParticleList.targets`, `keys`, `entities`, `compartments` and `radii` derivable from
   target ranges (an interface change in `particle-list.ts`, its readers checked), and what is built from them per
   particle by range too: the label indices (`getParticleLabel`), the targets' particle sets, the hierarchy color
   theme's index. No format pays for them per particle (1.16 GB here). Done: `ParticleValues` (an array or runs, read
   by `get`, `forEachRun`, `forEachRun2`, `map`), the pack's list of runs and no keys; the target groups an
   `Interval` a target of one run, the color themes' indices runs, a target's radius its poses' uniform `scale`
   (`InstancePoses.scale`); the other formats' lists arrays as before. The yeast cell's pack in Node: array buffers 1.88
   -> 1.25 GB (the file alone); in the Viewer: the heap 3.77 -> 3.17 GB at its peak, the list in 12 ms (172), the
   targets in 7.5 s (10.6), the image the same. The label indices are the mmCIF reader's (the pack's label reads the
   runs).
5. **The GPU from position and rotation** ([gl/instance-poses.md](gl/instance-poses.md)): the targets' instances from
   the list's views, 32 B a pose on the GPU against the 64 B `transform` and the 48 B rows. Not the file's bytes: a
   pose has its target's center and the representation's similarity folded in, and each part has poses of its own (a
   list's targets' particles are scattered in it), so the poses are derived (`instance-poses.md` step 3), streamed to
   the GPU in blocks of 2^18 ids through one reused buffer (8 MB; `CulledData.layout`): no CPU copy of the list's
   size retained. The rotations snorm16 in the file and the list (step 8), float32 in the poses.
6. **The sites block** (above): the float32 sites quantized, `ParticleSites` in its quantized form,
   the sphere impostors drawing quantized centers (the same quantized-buffer pattern as step 5), the CPU readers
   decoding where they read. The format done (version 2): the `site` block int16 × 4, a target's sites' origin and
   scale (the box's center, its larger half extent over 32767), the radii by type (`siteRadii`, the distinct radii);
   `ParticleSites.quantized`, read a site at a time (`getPosition`, `getRadius`: the geometry, the sizes, the
   boundary, the center, the converter), `positions` and `radii` decoded on each read. The yeast cell's pack 1.25 ->
   1.09 GB (MG_1184 112 -> 89 MB, PetWorld 73 -> 44 MB), its sites within 0.0095 Å (MG_1184 0.0046 Å), the radii
   exact; in the Viewer the heap 2.56 GB after loading (3.17 GB measured at step 4, another session's: the pack's
   162 MB the part accounted for), the image the same. Done toward it: the items drawn in their order with levels of detail (a level every
   stride-th, `screen-space-lod.md`: no reorder, so no copy needed to draw from the block). Done: the impostors drawing from the block (`Spheres.fromSites`,
   `sphere-sites.ts`): a sites target's spheres the pack's site view as it is (`tSite`, an `rgba16i` data layer: WebGPU two
   words an element, WebGL2 an `RG32I` texel; `dSites`), the center `origin + q * scale` (unmerged by object uniforms,
   merged by the member's in `tSegment` and the culler's segment data, `CulledData`), the group the site; the raster and
   its resolve decode them as well; WebGL1 (no integer data) keeps a float copy (`Spheres.isSitesSupported`). The
   yeast cell (`merging/yeast-pack`, whole at 720p, WebGPU): the scene's CPU geometry 310.5 -> 0.1 MB (the sites the
   pack's view, 155 MB, counted with the download), its GPU data 531 -> 376 MB, the JS heap it retains 1.80 -> 1.50 GB,
   its representations 6.6 -> 5.9 s, the frame 6.34 -> 6.40 ms; WebGL2 the GPU's 2.93 -> 2.78 GB, the heap 6.14 -> 5.83
   GB, 15.7 -> 14.2 s. PetWorld from its pack (`merging/postsyn-pack-lod-pixel3[-zoom]`) looks as from the BinaryCIF.
   The sizes from the radius table, done: a sites target drawn from its sites sized by their types (`siteType` sizes,
   `createSiteTypeSize`): `tSize` the table of the radii by type (the theme's scale applied, a few entries), read at a
   site's type (its uint16 in `tSite`), merged the members' tables packed (`table` granularity, a member's first entry
   in `tSegment` and the culler's segment data). With the default (a uniform size) there are no sizes a site; with the
   radii (`merging/yeast-pack-physical`, WebGPU) the yeast cell's GPU data 453 -> 376 MB, its CPU layers 228 -> 150 MB,
   the JS heap the scene retains 1.57 -> 1.50 GB, the image PetWorld's (`merging/postsyn-pack-lod-pixel3-physical`) the
   same. Found on the way: on WebGL2 the yeast cell drew nothing, its merged object's instance data texture
   (`tInstanceData`, 80 B an instance: 2.5 GB) beyond what a texture upload takes (fixed: 64 B an instance in pages of
   256 MiB, [merged-render-objects.md](gl/merged-render-objects.md)); merged objects' group indices were floats (`aSegment.z` + group), not
   exact above 2^24 groups (the yeast cell's sites sized by group read their neighbours' sizes there; by type exact):
   integers since [integer-ids.md](gl/integer-ids.md) steps 1 (`aSegment` int32) and 4 (the group), WebGL1 floats.
7. **A scene-wide instance grid** (graphics, independent of the format, for every large scene). Today each target
   (a merged object's member) grids its own instances, and for sparse targets spread over the whole scene the grids
   degenerate: the yeast cell's are about 20 M cells for 31.5 M instances (`cellSpheres` 306 MB at 16 B a cell,
   `cellOffsets` 77 MB). Three things grow with them: the grids (551 MB, and the `cellTransform` copies 302 MB), the
   merged culler's cell table rebuilt from them (`merged-culler.ts`: a sphere a cell, `cellOf` 8 B an instance; likely
   most of the unattributed 377 MB) and the multi-draw lists, sized by the members' cells times the regions and levels
   (`merged.ts` `getCapacity`: 2.15 GB). Instead:
   - **One grid over all instances**, by position, not by target: 10⁵-10⁶ cells (a few MB of spheres), the CPU paths'
     frustum and occlusion tests 20-200× fewer.
   - **A cell's contents as runs** of (target, first, count): with the pack's particles sorted by target, and within a
     target by scene cell (Morton order), a run is a sub-range of the target's range, so step 4's ranges hold and no
     cell index an instance is needed.
   - **Size classes**: a cell's sphere covers its contents' extents, so one holding a 300 Å target is padded by it (the
     per-target grids size their cells by the target); a scene grid per size class (or a cell padded by its largest
     content) keeps the culling tight.
   - **Draws**: the GPU culls and compacts each target's visible instances (`gl/culled-instance-indices.md`), one
     indirect draw a target and level (about 5,639 × the levels); the CPU path (WebGL) draws from the visible cells'
     runs, its lists sized by what is visible, not by a capacity. The multi-draw lists' 2.15 GB mostly gone, the
     shifted copy made where read.
   - **Changes by cell**: a moved, added or removed target updates the cells holding it (particle drags,
     `Particle.transform`), a hidden one is a visibility the culler reads, not a grid change; levels of detail stay per
     target (the runs carry it).
   - **The lists by use** (done): `getMdbData` starts a list at 1,024 entries and `emit` grows it by half again,
     so a list holds what a cull emitted, not a draw a cell of every member's cells; the GPU path emits none. The
     yeast cell: the renderer at the first frame 11.21 -> 9.11 GB, the peak 11.72 -> 10.03 GB, the first frame
     unchanged (42.8 s); the visual suite, the guards unchanged. A CPU-culled frame (a fallback) grows them to what
     it drew, and they stay at that size.
   - **The unattributed closures**, named (a build without minification): the culler's `layoutSpheres` (377 MB) and
     `cellSegments` (94 MB), the member cells split 4×4×4 (`splitCell`, about 24.7 M cells), and a label string table
     (78 MB). The scene-wide grid takes the first two.
   - **The GPU culler's cells from one grid** (done): `buildSceneCells` (`merged-culler.ts`) in place of each member's
     grid cells split 4×4×4 (`splitCell`): a uniform grid over all instances, about 32 instances a cell
     (`MergedCullerLimits.sceneCellInstances`), a member's instances in the size class whose cells hold its sphere (4
     classes, cells 2^k times the side), a cell's sphere the box of its instances' centers expanded by their largest
     radius (the members' invariant radius by their largest scale); a move of all members alike moves the cells, any
     other is laid out again. Since built on the GPU (`scene-cells.ts`, [load-time.md](gl/load-time.md) step 4: each
     instance's own radius, any move a build again). The yeast cell (heap snapshots): buffers of at least 1 MB 6,870 -> 6,406 MB (`layoutSpheres`
     and `cellSegments` gone), the first frame 41-45 s (as before), the image unchanged; MG_1184 the same frames (see
     `merged-render-objects.md` step 8).
   - **The targets' grids coarse where the GPU culls** (done): a target spread over the scene had about a cell an
     instance at 200 Å (20 M cells, 551 MB). Where the representation is culled on the GPU (`isTargetCulledOnGpu`: the
     merge's rule applied ahead), a target's cells grow to hold about 32 instances (`getTargetCellSize`): its grid
     serves only the frames culled on the CPU. Not where the CPU culls: coarse cells cull less (MG_1184 culled on the
     CPU, zoomed: the opaque pass 6.3 -> 18.6 ms), though drawn fewer times (whole: 15.7 -> 5.2 ms). The yeast cell:
     buffers 8.07 -> 7.56 GB (the grids' cells 1 MB), the image the same; MG_1184 both paths as before.
   - **Every view culled on the GPU**: the CPU culls where the GPU did not (WebGL; kernels not ready, frames held
     meanwhile; the clip ghost where its region was not dispatched; objects past the raster budget with `sphereRaster`
     `auto`), and for a frame's second view: GPU culling runs once a frame, for its first view (`isGpuView`), so
     stereo's and XR's second view cull on the CPU every frame (the yeast cell: ~20 M member cells, its lists grown to
     all it draws), and a scaled camera (XR's world scale) skips it (`gpu-culling.ts` `updateScale`). A dispatch a view,
     or one for their joint frustum; the scene-wide grid makes the CPU culls that stay 20-200× cheaper.

8. **snorm16 rotations** (done): the `rotation` block snorm16 × 4 (`ParticleRotations`: a unit component times
   32767, rounded, `w >= 0`; pack version 3, version 2's float32 still read), the list's `rotations` that `Int16Array`
   as a view of the bytes. The readers take either kind as it is and rotate by `2 / |q|²` in place of 2 (exact for
   any length: the core kernels, `writeParticleTransform`; `ParticleRotations.get` where a unit quaternion is
   wanted); a move (`Particle.transform`) gives float32. Within 1e-4 rad (the tests). The packs: the yeast cell 1,089 ->
   838 MB, MG_1184 89 -> 72 MB, PetWorld 44 -> 40 MB; PetWorld's image the same by eye (the WebGPU image of the
   BinaryCIF against the commit before: 3% of the pixels, depth ties between overlapping spheres, mean 0.07 of 255; the
   commit before against itself identical). The GPU keeps float32 poses (32 B, normalized where written). Quantized
   there, the reader runs in every raster chunk and every vertex, its decoding in their latency: snorm16 quaternions (24 B, two ids in three vec4s) saved 240 MB of the yeast cell's
   3.6 GB of storage, its sphere raster 5.9 -> 7.9 ms, its opaque pass 0.19 -> 0.55 ms, MG_1184's `whole-cpu` opaque
   pass 3.9 -> 4.9 ms (two runs each); one aligned vec4 (16 B: the translation 21 bits an axis on a grid of the
   buffer, the quaternion smallest-three, the scale a half float where a part's one) saved 480 MB, its raster 7.7 ms,
   its opaque pass 0.53 ms; the reader's code alone (the 16 B branch compiled, not taken) 6.0 and 0.21 ms, so the
   decoding, not the reads. The rows built once a chunk in the raster (lane 0's, shared through
   workgroup memory) do not help: a chunk's lanes run in lockstep and read one address as one transaction, so a lane
   builds them in the time all do, and the barrier waits on it (the yeast cell's raster 5.9 -> 6.2 ms, MG_1184's as
   before). As kept: the yeast cell's frame 6.35 ms as before, MG_1184's
   cases within their runs' noise, the yeast cell's add and commit 3.0-3.2 -> 3.3-3.5 s (the quaternions normalized
   where written). The JS heap the scene retains unchanged (1.50 GB): the rotations are the download's bytes, 252 MB
   fewer.

9. **Positions on a grid** (done): the `position` block 21 bits an axis in two u32s over the list's box
   (`GridPositions`, `@molstar/core/math/geometry/grid-positions`: `origin + step * u`, the grid in the document;
   pack version 4, versions 2's and 3's float32 still read), the list's `coordinates` that grid as a view of the
   bytes (`Positions`: float32 or a grid). The loops over all positions (the core kernels, `CulledData`'s poses)
   branch on the kind once and decode inline; the others read `GridPositions.get`; a move gives float32. Within half
   a step: the yeast cell's 0.011 Å, MG_1184's and PetWorld's below 0.002 Å. The packs: the yeast cell 838 -> 712 MB,
   MG_1184 72 -> 63 MB, PetWorld 40 -> 38 MB; the GPU's poses as before (float32, the positions decoded where
   written). Measured (two runs each, against the commit before): the yeast cell's frame 6.37 ms, its raster 5.9 ms,
   its add and commit 3.4-3.7 s (3.7), its heap the same (the positions the download's bytes, 126 MB fewer);
   MG_1184 (float32 positions, the branch taken) its frames and setup as before; PetWorld's image the same by eye.

## Decisions open

- **The name and extension** (`MSPK`, `.mspk`).
