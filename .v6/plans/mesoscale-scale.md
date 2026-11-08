# Plan: mesoscale models at scale

## Status

- **In progress**: step 4 done (29.4 M isg particles load, 9.1 GB); next the yeast milestone (below). The GPU side
  holds past 8 M instances (`gl/culled-instance-indices.md`); a structure of every placement does not, particles do.
  The mesoscale explorer is out of scope (it is to be rewritten): mesoscale models load as particles.

## Question

The isg model repeated on a grid (`isg_mature_grid-AxBxC.bcif`, 460k placements of 100 entities a cell): the mesoscale
explorer loads 1 cell (460k) and runs out of memory at 8 (3.68 M): Chrome's renderer reaches 4.9 GB and dies. The
explorer builds the cellPACK assembly as a `Structure` (`buildCellpackAssembly`), a `Unit` and a `SymmetryOperator` per
placement, about 1.7 KB of JS heap each (Node, Mol*'s libraries, the explorer's steps):

| step               | 1 cell (460k) | 8 cells (3.68 M) | per placement |
| ------------------ | ------------- | ---------------- | ------------- |
| parse and model    | +1 MB         | +1 MB            | typed arrays  |
| assembly operators | +436 MB       | +3492 MB         | ~950 B        |
| units              | +345 MB       | +2758 MB         | ~750 B        |
| structure          | +28 MB        | +224 MB          | 64 B          |

V8's heap with pointer compression (Chrome) holds at most 4 GB a page: about 2.3 M placements, whatever the GPU does.
The operators: three `Mat4` arrays each (`getMatrices`' map by id, the composed one, an eager inverse), strings
(`ASM_n`, the canonical name, the suffix, `operList`), the expanded expression lists; past ~8.4 M a `Map` of the
matrices fails outright (3Ã—3Ã—3: `Map maximum size exceeded`).

The Viewer's particle format (`mmcif_particles`, `createParticleListFromMmcifAssembly`) reads the same file as a
`ParticleList`: a particle per placement in typed arrays (target, entity, compartment, position, rotation, radius:
about 48 B), the targets a reference `Structure` per entity (100 here), drawn by target representations (`spacefill`,
`blob-surface`) instanced over the particles. In Node its JS heap stays at 13 MB for 1, 8 and 27 cells (12.4 M
particles in 2.6 s); the typed arrays are mostly the parsed CIF (the decoded operator columns, ~96 B an operator).

## Design

Particles decorated with target representations are what mesoscale models are: no structure of every placement. The
particle path carries them; what is left is its targets and its own limits at scale.

- **Targets as a minimal interface.** A target's geometry needs positions, radii (or element symbols for the size
  theme), group ids for picking and color, and its center: not a model and structure of its atoms (the trajectory,
  the sorted columns, the hierarchy, `Structure.boundary`; 0.6â€“0.8 s of PetWorld's and MG_1184's load,
  `mesoscale-load.md`). A `ParticleTarget` kind `atoms`: typed arrays read from `atom_site` (and the particles' own
  data) directly, `PositionData` for the density (`computeGaussianDensityTexture` takes it) and centers and groups for
  sphere impostors (`Spheres.create`). The `structure` kind stays for targets that need one (labels, interactions,
  a structure representation of a target opened on its own); a target becomes a structure on demand, from the same
  rows.
- **A Gaussian surface target** next to `spacefill` and `blob-surface`, from the minimal interface: a target's
  surface computed once, instanced over its particles.
- **The assembly path's own costs** (other large assemblies loaded as a structure): the operators' matrices in a flat
  `Float32Array` by row (`getMatrices`, no `Map`: the 8.4 M limit), the inverse lazily, names when read. Optional.

## Steps

1. **The particle path in the browser at scale**, measured (Chrome, WebGPU, the Viewer's `mmcif_particles`, a fresh
   browser each; the driver `explorer-scale.mjs --app viewer`, scratch):

   | cells | particles | to the first frame | page memory (heap and buffers) |
   | ----- | --------- | ------------------ | ------------------------------ |
   | 1     | 460k      | 4.9 s              | 0.6 GB                         |
   | 8     | 3.68 M    | 6.6 s              | 2.4 GB                         |
   | 18    | 8.28 M    | 10.3 s             | 4.9 GB                         |
   | 27    | 12.4 M    | 13.8 s             | 7.3 GB                         |
   | 64    | 29.4 M    | fails at 18.5 s    | 14.8 GB (`Array buffer allocation failed`) |

   What stays (8 cells, a heap snapshot after the first frames, the array buffers of at least 1 MB by owner; scratch
   `snapshot-buffers.mjs`): 1.98 GB, 540 B a particle, mostly copies:

   | owner                                                           | per particle |
   | --------------------------------------------------------------- | ------------ |
   | the particle list's cached transforms (`particle-transforms`)   | 64 B         |
   | the merged object's `aTransform`                                | 64 B         |
   | the members' (targets') `transform` values                      | 57 B         |
   | the members' instance grids' `cellTransform` (reordered)        | 55 B         |
   | the parsed file's decoded columns                               | 89 B         |
   | the asset cache's file bytes (`inflatedCache`)                  | 31 B         |
   | the parsed operator columns (`matrix`, `vector`)                | 28 B         |
   | the particle list (keys, targets, entities, ..., rotations)     | 46 B         |
   | the merged `aSegment` (per member data per instance), `aInstance` | 20 B       |
   | `getParticleLabel`'s indices                                    | 11 B         |
   | the CPU lists' copies (`shiftedLists`), the color index         | 18 B         |

   The transforms four times (240 B, one copy 64 B, or 28 B as position and rotation); the parsed file kept after
   the list and targets are built from it (150 B). Loading peaks above what stays (2.4â€“2.9 GB at 8 cells): 64 cells
   failed while building.
2. **Minimal targets**: the `sites` kind from `atom_site`, `spacefill` and `blob-surface` from it, the `structure` kind
   on demand; the load measured against `mesoscale-load.md`'s (PetWorld, MG_1184) and the images unchanged.
3. **Gaussian surface target** from the minimal interface.
4. **The particle path's copies** (step 1). Done, what stays at 8 cells 1981 -> 1181 MB (540 -> 320 B a particle), the
   largest process 2.95 -> 1.67 GB:
   - the targets' transforms written from each particle's position and rotation (`writeParticleTransform`), not
     copied from the list's of all particles (`getParticleTransforms`, cached on the list: kept for the particles' own
     visual), -224 MB;
   - the merged object's `aTransform` its members' as parts (`SegmentedArray`), not a copy; its rows packed where read
     (`packAffineRowsOf`: the instanced attribute of objects culled on the CPU, WebGL2's instance data, the culler's
     data, the raster), -225 MB;
   - the operators' columns read once (`CifCategory.readField`, `readCifField`), not kept decoded in the parsed file
     (its categories' field cache, kept by the models' source data), -351 MB.

   At scale (Chrome, WebGPU, to the first frame, the largest process):

   | cells | particles | before           | now              |
   | ----- | --------- | ---------------- | ---------------- |
   | 27    | 12.4 M    | 13.8 s, 7.47 GB  | 13.8 s, 4.92 GB  |
   | 50    | 23.0 M    |                  | 24.6 s, 8.02 GB  |
   | 64    | 29.4 M    | fails at 14.8 GB | 26.2 s, 10.95 GB |

   The loading peak (64 cells, the largest process sampled every 2.5 s): 10.95 GB against 9.05 GB kept, the extra
   while the targets were made (the culler's sync), not while the list was built. Done:
   - the culler's instance data written from its sources (`writeCulledData`: the rows packed a block of instances at a
     time, `aSegment` and `aInstance` as they are), not composed into one array first (68 B an instance): the peak
     11.6 -> 9.07 GB, what stays 9.03 GB;
   - the operators read as their columns (`ParticleOperators.columns`), not composed into 16 values an operator as
     well (96 instead of 224 B an operator while the list is built).

   Left: the members' `transform` and their grids' `cellTransform` (57 B each, the latter the merged object's parts),
   the downloaded bytes (31 B, the state's data), `aSegment` per instance (16 B), the label indices (11 B). Chrome's GPU
   process holds the uploads until they are processed: all processes 16.1 GB while the targets upload at 64 cells,
   13.4 GB after.
5. **The assembly path's operators** (optional).

## Milestone: the yeast cell

`yeast_cell_atoms.bcif` (`.benchmark/`, not checked in; 1.17 GB): 31.5 M operators (particles), 20.3 M target atoms
over 5,639 entities, 857k gen rows. The Viewer builds the list (9.2 s) and the targets (37 s), then fails at the commit
and first frame (`Array buffer allocation failed`, the renderer at 14.7 GB). The 64 isg cells failed the same way at
14.8 GB before step 4: a renderer's practical limit about 15 GB (Chrome, 96 GB of memory). The goal: about 9 GB, a
margin for transients and for models a little larger.

| owner                                                            | size     | how known  |
| ---------------------------------------------------------------- | -------- | ---------- |
| the file and the MessagePack copy of every column                | 2.2 GB   | Node       |
| the particle list (48 B a particle) and its label indices        | 1.75 GB  | Node       |
| the target models (20.3 M atoms, ~143 B an atom) and their heap  | 3.7 GB   | Node       |
| the targets' transforms and their grids' copies (64 + 64 B)      | 4.0 GB   | estimate   |
| `aInstance`, `cellInstance`, the merged `aSegment` (28 B)        | 0.9 GB   | estimate   |
| the targets' geometry (20.3 M spheres)                           | ~1 GB    | estimate   |

The GPU side's own limit: the culler's instance data (`CulledData`, 68 B an instance) is 2.139 GB in one buffer
against WebGPU's 2.147 GB (`maxBufferSize`): it fits by 0.4%.

Steps, each measured on the yeast cell and the isg grids (what stays, the peak, the load time), the images of the
loaded scene unchanged:

1. **One transform copy a target**, done: `createTransform`'s `oneCopy` (the particle targets'), `transform` the
   grid's `cellTransform` (the same array as `aTransform`, row `i` instance `aInstance[i]`, `isOneCopy`); a rebuild
   (`updateTransformData`: a matrix, an extra transform) restores the instances' order first, a move moves the grid's
   copy, an update writes new arrays (the grid's still read until rebuilt). 8 cells 1181 -> 972 MB. The yeast cell
   loads: the first frame at 74 s, 13.0 GB kept, the peak 13.9 GB (all processes 18.5â€“23 GB: the GPU process holds
   about 10 GB), close to the limit: steps 2 and 3 are still needed for a margin.
2. **Minimal targets** (step 2 above), built: the `sites` target kind (`ParticleSites`: positions, van der Waals
   radii, a label, the chains as runs; the centroid as the target's center), each chain's (cellPACK) or model's
   (PetWorld) read from `atom_site` (`getAtomSiteGroups`), not a model and a structure each (the trajectory's sorted
   columns, the hierarchy). A model's sites by chain as a model holds them (by entity, then `label_asym_id`, else
   `auth_asym_id`; reordered only where out of that order, a stable counting sort); its chains as runs (an offset and a
   label a chain, for the bounds and later for coloring), found by comparing the columns' values (`areValuesEqual`), a
   string read a run. Drawn as a structure's elements are (`target/sites.ts`: impostors, a sphere mesh, the blob surface
   over the positions, `computePositionBlobSurface`), the structure params and auto quality by count
   (`getElementCountQuality`), bounded as its target structure was (`computeStructureBoundary` over a chain's box
   sphere each, a unit each: the impostors' level of detail goes by it). Every render value of the targets as before
   (PetWorld's 150, MG_1184's compared). The yeast cell 12.44 -> 11.11 GB kept, the peak 13.46 -> 11.72 GB, the targets
   37 -> 17 s, the first frame 64 -> 43 s. Chrome, the Viewer, three runs each against the build before:

   | model    | parse                    | targets               | kept            | peak heap       | image     |
   | -------- | ------------------------ | --------------------- | --------------- | --------------- | --------- |
   | PetWorld | 1.03–1.15 -> 0.47 s      | 445–452 -> 420–440 ms | 0.87 -> 0.74 GB | 0.93 -> 0.66 GB | unchanged |
   | MG_1184  | 1.36–1.39 -> 0.87–0.90 s | 930–969 -> 833–913 ms | 1.01 -> 0.88 GB | 1.18 -> 0.96 GB | unchanged |

   PetWorld's chains cost about 100 ms of its parse (3.7 M rows: two columns decoded, 25 ms, and compared a row, 70 ms);
   cellPACK's none (a target a chain). Not exact: radii as float32 move a packed size by one step (0.01 Å) at a rounding
   half. Open: a target's structure on demand.
3. **The culler's segment as an index**, done: `CulledData` holds each instance's segment index (4 B) and each
   segment's groups (a vec4, `MergedSegments.groupCounts`, `groupInstanceOffsets`, `groupBases`), `culledSegment`
   composing `aSegment` from them: 56 B an instance (38 M in a buffer; the yeast cell's 1.76 GB, 18% under
   `maxBufferSize`). The merged `aSegment` made where read (a box whose value is computed once read: objects culled on
   the CPU, WebGL2's instance data), so a culled object never makes it. 8 cells 972 -> 916 MB; the yeast cell 13.03 ->
   12.44 GB kept, the peak 13.94 -> 13.46 GB, the first frame 74 -> 64 s.
4. **The file's bytes**: the MessagePack columns as views of the downloaded bytes (`decodeMsgPack` copies on purpose:
   element access, `view(data.buffer)` callers; measure), or the bytes released once parsed. -1.1 GB.
5. **The label indices** made where a label is asked for. -0.25 GB.

Steps 1 to 3 bring it to about 9 GB by the estimates; 4 and 5 add margin.
