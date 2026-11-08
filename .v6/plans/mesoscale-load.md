# Plan: mesoscale load time

## Status

- **Done**. Final (as below, `prepareObjects` on, three runs; the build to the commit, then from it to the first full
  frame): PetWorld parse 1.24–1.31 s, build 0.49–0.51 s, to the full frame 0.24–0.26 s, about 2.0 s in all (3.0–3.7 s
  before); MG_1184 parse 1.54–1.58 s, build 1.02 s, to the full frame 0.48–0.50 s, about 3.1 s (5.1–5.7 s before).
  Closed with each cost left at most about 120 ms (3–6% of its stage, under the runs' noise of about 10%, so that a
  change to one cannot be shown in the browser: `splitCell`'s, measured and dropped). Left, if ever taken up again,
  measured in Node: the particle targets' own copy and centering of the transforms (113 ms of MG_1184: per target from
  the list's transforms, `getParticleTransforms`, where the particles' positions and rotations could make them
  directly), `splitCell` (168–184 ms in the profile), the uploads (115–130 ms), the instance grid (60 ms), the box
  sphere (55 ms), the transforms' rows (48 ms), the particle hierarchy color theme's index (46 ms); in the parse the
  list's own loop (about 200 ms) and the target models' work per atom (the sorted columns' copies, the coordinates'
  decoding, the hierarchy). A structure target is placed by its centroid (`ParticleTarget.getCenter`; the mmCIF
  formats give the centroid their particles were placed by), not its geometry's bounding sphere center, which had
  misplaced the atoms by up to 126 Å (PetWorld) and 92 Å (MG_1184); a volume target by its box center, not its
  isosurface's. Pipelines compile while it builds (`gl/pipeline-readiness.md` step 7d).

## Where things stand

- **Measured before the steps** (Chrome, RTX 4070 Ti, 1280x800, the Viewer, WebGPU, `mmcif_particles` format: parse,
  then visuals; three runs, a fresh browser each): the parse 1.4–1.8 s (PetWorld `postsynstubby_2`) and 2.0–2.3 s (MG_1184
  cellpack), the build to the commit 1.2–1.6 s and 2.1–2.4 s, the commit's own work before the first frame 255 and
  880 ms, the first frame 70–90 ms.
- **Profiled** by stage (the Viewer built with `--sourcemap`, CPU profiles mapped to the sources), inclusive ms,
  PetWorld / MG_1184:

  | stage, what                                                                  | PetWorld | MG_1184 |
  | ---------------------------------------------------------------------------- | -------- | ------- |
  | parse (busy)                                                                 | 1363     | 2093    |
  | the target structures' model (`readStandard`; `sortAtomSite` 184)            | 835      | 602     |
  | the particle list (`createParticleListFromMmcifAssembly`)                    | 324      | 809     |
  | cellpack's particle list (`buildCellpackStandardParticleList`, operators 296) | –        | 969     |
  | binary CIF decoding (`decodedField`)                                         | 315      | 427     |
  | visuals and commit (busy)                                                    | 1701     | 3613    |
  | render objects (`createRO`)                                                  | 548      | 1561    |
  | of it, the instances' bounding sphere (`calculateTransformBoundingSphere`)   | –        | 616     |
  | of it, the particles' transforms (`createTargetParticleTransform`)           | –        | 519     |
  | target geometry (`createStructureElementSphereImpostor`)                     | 784      | 700     |
  | of it, the structures' boundaries (`Structure.boundary`)                     | 435      | ~650    |
  | of it, the largest size by a size theme lookup per atom (`VdwRadius`, ...)   | ~160     | –       |
  | the sizes (`createStructureTargetSizes`)                                     | 187      | –       |
  | the canvas commit (`Canvas3D.commit`)                                        | 258      | 960     |
  | of it, the GPU culler's sync (`MergedCuller.sync`, `syncInstances` 491 MG)   | 215      | 899     |

- **Bounding spheres by extrema** (`BoundaryHelper`, 49 directions, `computeSphereExtrema`): the structures' boundary
  per unit and the transforms' sphere over every instance, together about 435 ms (PetWorld) and 1.3 s (MG_1184). The
  perf audit's 02/P6 found the coarse helper's spheres no longer enclose; what the targets need here is a sphere that
  encloses (culling, levels), not the tightest.

## Steps

1. **The instances' bounding sphere** (`calculateTransformBoundingSphere`), done: from 1000 instances, for every
   renderable, the sphere around the box of the instances' centers, its center and the largest distance plus the
   invariant radius (two passes). As tight or tighter than the helper's there (a 3PQR lattice 8% smaller, the 1M4X
   capsid 1%); not below, where few instances spread wide (a dozen over the cell's surface) made MG_1184's scene sphere
   1.5x too wide. 616 → 55–70 ms of MG_1184 (2.1 M instances over 848 targets, the 38 above 10k most of them).
2. **The target structures' bounding sphere** (done: fast boundaries); the size lookup per atom (done): the physical
   size theme reads the vdW radius by the model's atomic numbers (`VdwRadiusOfNumber`), not by a lookup of the element
   symbol, 157 → about 65 ms on PetWorld (the theme's own call per atom left). An element's serial index the unit's
   first plus its position (done), not a binary search an element (`getSerialIndex`: the sphere and point builders),
   and the targets' sizes written unit by unit, not a location an element through the serial mapping (which then is
   not built at all): the impostors 354 → 181 ms, PetWorld's visuals busy 1031 → 947 ms. The spheres' level of
   detail order (`forEachInLevelOrder`, done): the items of each stride marked (no modulo an item and level) and placed
   by a counting sort in one pass (not a pass over the items per level), the same order: about 120 → 62 ms on
   PetWorld.
3. **The GPU culler's sync at the commit** (`MergedCuller.sync`), 934 → 455 ms of MG_1184: the instance copy reads
   the identity region, not a packed copy of it (`instanceData`, 80 B an instance uploaded again); the cells' spheres
   into a typed array, not pushed (1.4 M cells); no array and no closure per instance. Left: `splitCell` (168–184 ms
   in the browser's profile), the upload (115–130 ms). Measured and dropped: `splitCell`'s work per call kept to the
   buckets a cell uses (half of MG_1184's 97k split cells hold 2 to 9 instances, the 64 buckets cleared and scanned per
   call), in Node 131 → 85 ms on MG_1184's cells, in the browser no gain beyond the runs' noise; the instances sorted by
   bucket first (each box a run in locals), slower (102 → 128 ms in Node).
4. **The particles' transforms** (`createTargetParticleTransform`), 510 → about 400 ms of MG_1184 in the browser
   (`createTransform` 2.3x faster in Node: 654 → 284 ms for 2 M instances): the reflection test, the grid's centers and
   cell boxes inlined (no `Mat4`, `Vec3`, `Box3D` calls and no `Sphere3D.fromBox3D` extrema per instance or cell); the
   grid reads `transform` (no copy into `aTransform` first, its `cellTransform` replaces it); a new array handed over
   (`createTransform`'s `adopt`), not copied. Left: the particle code's own copy and offset (about 120 ms), the grid's
   reorder and lookup.
5. **Parse**: the particle lists, in Node 1584 → 1112 ms (MG_1184) and 660 → 313 ms (PetWorld), the lists the same
   to the bit: the operators read as whole columns a matrix at a time (not a value per accessor call), operator
   expressions to operator indices directly (`parseOperatorCombinations`: a range not a string per id), a gen row's
   chains read once, a chain's target, compartment and entity once a chain (not a particle), the atom rows' keys once
   (a chain's entity from its first row). The targets' model (3.7 M atoms over 150 models in PetWorld, 2.9 M in
   MG_1184), in Node 1940 → 1201 ms and 2341 → 1336 ms, the models the same (element symbols, atomic numbers, trace
   and direction atoms, molecule and polymer types): a BinaryCIF `type_symbol` kept a dictionary (its values mapped, not
   a string per atom), the atomic numbers once a value, a residue's trace and direction atoms by their name's
   dictionary index; in the browser the parse 1.40 → 1.22–1.28 s and 1.64–1.79 → 1.55–1.60 s. The models must be
   sorted: PetWorld's chains (lipids) and MG_1184's entities interleave. Left: the list's own loop (composing, the
   quaternion, about 200 ms in the browser), the operators' decoding; in the model the sorted columns' copies (150 ms),
   the coordinates' decoding (166 ms), the hierarchy's offsets and data (130 and 155 ms).

Each step measured with the stage profile (scratch `load-profile.mjs` against `node scripts/esbuild/app.mjs apps/viewer
--sourcemap`) and the load times, images of the loaded scene unchanged.
