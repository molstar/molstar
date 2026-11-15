# Plan: instances placed from poses (particle targets)

## Status

- **In progress**: steps 1 to 3 done, step 4 on WebGPU (WebGL2 later). Step 5 of [particle-pack.md](../particle-pack.md); the baseline below measured.

## Question

A particle target's instances are particles: each a position and a rotation (the list's `coordinates`, `rotations`),
times a scale (`scaleByRadius`: the particle's radius) and minus the target's center (the point a particle places,
`ParticleTarget.getCenter`). Today each instance becomes a 16-float matrix on the CPU (`createTargetParticleTransform`,
`TransformData.aTransform`, one copy in the grid's order) and three affine rows on the GPU (an instanced attribute, the
WebGL2 `tInstanceData` texture, the GPU culler's `CulledData`: rows 48 of its 56 B an instance). The yeast cell (31.5 M
instances): about 1.5 GB of matrices in the renderer's heap (`cellTransform` 981 MB, `aTransform` 452, `transform` 59;
heap snapshot), 1.5 GB of rows in `CulledData` on the GPU, written from them at load (`writeCulledData`).

## The encoding

- **A pose an instance**: position (vec3) and uniform scale (f32) as one vec4, rotation as a unit quaternion (vec4):
  32 B, where the matrix is 64 B on the CPU and 48 B on the GPU. Per target (per member of a merged object) its center
  (vec3), in the segment's data (`tSegment`, `segmentData`) or an object uniform: model = T(position) · R(q) · S(scale)
  · T(-center). The representation's `matrix` and a move compose into the poses (rigid, uniform scale), as now into the
  matrices; a mirroring `matrix` keeps the matrices (`hasReflection`, below).
- **Read by index where the GPU can**: the poses of the whole list in two storage buffers (positions with scale,
  quaternions), each instance's particle its member's first particle plus its instance id where the target's particles
  are a range (an `Interval`, the list sorted by target: the pack's order, step 4), else a particle index an instance
  (4 B). On WebGPU no data an instance beyond what `CulledData` has besides the rows (`aInstance`, the segment): 56 to
  8 B an instance, plus 32 B a particle once for the list (the pack's blocks, written to the GPU from the file's bytes).
- **Attributes where it cannot** (WebGL2 instanced, uncull): a pose an instance as two vec4 instanced attributes (32 B,
  in the grid's order), or the `dInstanceData` texture's two texels (4 -> 3 texels an instance). WebGL1 keeps matrices
  (8 attribute locations, no texel fetch).
- **Applied as a pose, not as a matrix**: world = position + scale · rotate(q, local - center), the rotation
  `v + 2 · cross(q.xyz, cross(q.xyz, v) + q.w · v)` (about 15 multiply-adds), no matrix built. An impostor needs its
  center (one rotated vector) and its radius scale (`scale`); a mesh rotates its position and normal by `q` (uniform
  scale: no adjoint); the cull kernels and the CPU readers an instance's center and largest scale (`scale` itself, not
  the square root of the rows' lengths). A matrix only where one is consumed: export (a matrix an instance), and the
  shader chunks that take a model matrix (instance clipping, tumble, volume instance grids, direct-volume's rows to its
  fragment stage), adapted or building it in those variants alone.
- **The CPU readers branch on the encoding, not call an accessor an instance**: those that read all instances run at
  build, update, a move or a theme (the bounds, the grids, the culler's segment scales, the position iterators, moves,
  export), none a frame: each takes the arrays and the encoding (matrices, or poses with the targets' centers) and runs
  one of two monomorphic loops, branched once outside them (as `writePoints` loops). A
  pose's loop costs about a matrix's (a rotation and a translation for a matrix-vector product) on half the bytes. The
  marked lists (`markedSegment`: the marked instances alone, on a change of marks) may read an instance at a time.

## Where poses are read (from the map of the readers)

- Vertex stages: spheres and mesh first (the targets' geometry: impostors, sphere meshes, surfaces), both backends;
  cylinders, lines, points (shape targets) after, or matrices kept for them. A define (`dInstancePose`) in place of
  `aTransform`, as `dInstanceData` replaces it; the parity test's names in both backends.
- The GPU culler: `instanceSphere` (cell-cull) from the pose (its center, the scale), `writeCulledData` without rows,
  `setSegmentScales` by the encoding; the scene cells built on the GPU from the culled data
  ([load-time.md](load-time.md) step 4).
- The sphere raster: `slotTransform` (merged: the culled data; unmerged: its own buffer) from the pose; the spheres'
  resolve reads the raster's transforms.
- Normals: rotation and uniform scale, so the mesh's adjoint stays exact; tumble animation seeded by the pose's center.
- `hasReflection`: poses never reflect; a mirroring representation `matrix` or move keeps the target on matrices.

## Steps, each measured against the baseline

1. **The readers by encoding** (no behaviour change): the CPU readers given the arrays and their encoding, matrices
   alone for now, their loops branched once. The same numbers expected (allocation guards, MG_1184, PetWorld, the
   yeast cell's load, where these loops run). Done: core's `InstanceTransforms` (kernels `writePoints`,
   `getMaxDistanceSq`, `getMaxScaleOfAll`; `getPoint`, `getDirection`, `getMatrix`, `getMaxScale` an instance) read by
   the bounds, the instance grid, the culler's scene cells and segment scales, the marked culling, the mesh and
   center-group position iterators (`getInstanceTransforms`, the one place a `TransformData` names its encoding).
   Not here: the moves and the target producer's box write matrices; for poses a move composes into the poses'
   `matrix` (no instance touched) and the box comes from `writePoints`, so both change with step 2.
2. **Poses on the CPU** for particle targets: no matrix array made; the grid orders an index (`aInstance`) and the
   poses' order, not 64 B rows. The GPU still fed rows, packed from the poses where read (`packAffineRowsOf`'s place),
   so this step is the CPU's memory alone: the yeast cell's heap -1.5 GB expected, the load faster by the matrices'
   writes and the grid's row copies. Built: spheres and mesh targets (`isPosesGeometry`) by
   `createPosesTransform`, the particle list's arrays read where they are; `TransformData.poses` (the base: the
   instances' order, moves composed into its `matrix`) and `instancePoses` (as drawn: the grid's order, `matrix`
   applied; `aTransform` empty, its version bumped with them); the grid from poses (`InstanceData.poses`, no
   `cellTransform`); the row producers (`affineRowsCell`, `writeCulledData`, the sphere raster, WebGL2's instance
   data) by `getTransformParts`, a merged object's members as parts (its `instancePoses` `segmented`); the readers
   that take matrices and run rarely (export, color smoothing, the debug helpers) by `getInstanceMatrices`. An extra
   transform or a mirroring `matrix` turns the poses into matrices. The kernels' pose loops inlined (fields read
   once, `local - center` once) and poses of one shape (`InstanceTransforms.ofPoses`): a call an instance made
   MG_1184's setup 3% slower. Measured: the yeast cell's buffers of 1 MB or more 6.24 -> 4.75 GB (`cellTransform`
   981 MB, `aTransform` 359 MB gone), its heap at the first frame 8.18 -> 6.31 GB, the renderer 8.57 -> 6.69 GB, the
   target made in 9.1 s (11.3), the first frames at 34.8 s (40.6); MG_1184's setup 6-11% faster than before the step,
   its frames as before (the step's tree and the commit before it run alike; both drifted from the baseline since,
   `whole-cpu` 13.2 -> 17 ms on both); PetWorld as the baseline.
3. **WebGPU from poses**: the culled programs (spheres, mesh) read the pose by index and apply it (`dInstancePose`), `CulledData`
   without rows, the cell-cull and raster kernels from poses. GPU storage for the yeast cell 1.76 -> ~0.25 GB plus the
   list's 1.0 GB of poses; the load without the rows' writes. Done otherwise: `CulledData` an instance a slot word
   (u32, a bit for a pose), a pose an id of its part (two vec4s: `t = p - s R c` with `s`, the quaternion; the part's
   matrix folded in where a similarity, else the part on rows), its reader (`culledRows`) building the rows, so the
   culled programs, cell-cull, the raster and its resolve unchanged (no `dInstancePose`); the raster's buffer of
   unmerged objects in the same layout. A pose an id of its part, not a list's once: the targets' particles are
   scattered in the list, a range a target repeated it (MG_1184: 5.2 GB, the buffer failed). Per segment a center and
   a matrix applied in the reader cost the opaque pass 10% and the raster 24% (nine reads a vertex); folded, a read
   more than rows (the slot word) and the quaternion's matrix. Measured: MG_1184's storage 170.8 -> 146.5 MB, its
   zoomed cases' opaque pass +2-4% and raster +1% (their GPU passes within 1% between runs; the whole-view and
   rotating cases 20-35% apart between runs of one tree, not judged), medians as step 2; PetWorld as the baseline; the
   yeast cell 56 -> 44 B an instance (an instance a particle: 1.76 -> 1.39 GB), its processes after the first frames
   11.36 -> 11.06 GB, the first frames at 33.8 s (34.8), its heap as step 2. A pose stays 32 B: quantized poses (24 B,
   16 B; `particle-pack.md` step 8) are decoded in every raster chunk and vertex, the yeast cell's sphere raster 31-34% slower.
4. **WebGL2 from poses**: two instanced vec4s (or the texture's texels); the unculled WebGPU path too. The WebGPU path
   done: objects of poses not culled on the GPU (unmerged, merged culled on the CPU; spheres and mesh) read their
   instances from `CulledData` of their own (`sInstanceData`, `instance-storage.ts`) at `instance_index` (the draw's
   first instance included), by the culled programs' reader without indices (`dInstanceStorage`; `usedIf` conditions
   joined by `&`: the instanced attributes unused with either define), synced with the renderable's update; set where
   the object has poses when made (matrices it falls back to read as rows); precompile by the signature's
   `instancePoses` (a particle target's by its geometry). Measured (MG_1184 culled on the CPU): attributes 137.8 -> 0 MB,
   storage 131.0 -> 220.2 MB (48.6 MB less), add and commit 140 -> 116 ms; the opaque pass 6-7% slower (a storage read
   and the rows built where vertex fetch read rows), medians -5 to +1% (CPU-bound); the GPU-culled cases as step 3.
   WebGL2 not done.
   The vertex stages after (spheres, mesh, image; both backends): the position transformed, not the matrices
   multiplied first (`uModelView * transform * position4` multiplied two 4x4 matrices a vertex, three with the model
   position and `uObjectTransform`), the mesh's normal by each transform's cofactors in turn, the sphere resolve's
   model position from the placed point. With the identity rotation in place of the quaternion's (the image wrong)
   the opaque pass was 4-7% faster, but the zeros let the compiler drop the products after; ordering the products did
   it with the image: MG_1184's opaque pass 2-6% faster than steps 3 and 4 (every case at or below step 2's, culled on
   the GPU or the CPU, the raster-only cases' resolve too), meshes' 1% (`globule-300k-spacefill-mesh`), the WebGL
   images unchanged (geometry cases against the tree before). A pose applied as such (`dInstancePose`: no matrix
   built) not done: the rows' cost was the products, not the quaternion's matrix. The other vertex stages the same (cylinders,
   lines, points, text, instance clipping; direct volume's box not, its fragment stage takes the matrix), GLSL
   cylinders' `vTransform` (a varying no fragment read) dropped: `7qpd-ball-and-stick`'s opaque pass on WebGL2
   0.087 -> 0.054 ms, on WebGPU as before; the WebGL images unchanged (geometry, clipping, shader variants).
5. **The other geometries** (shape targets: cylinders, lines, points) or matrices kept for them, decided by use.

## Baseline (this branch at 5d746d881; Chrome, RTX 4070 Ti, WebGPU, `--gpu-timing`)

MG_1184 (`sphere-raster-mg/*`, 24 cases) and PetWorld (`occlusion-lag/postsyn*-10deg-occlusion`) saved for
`--compare`; the yeast cell 40.6-41.6 s to the first frame, 7.56 GB of buffers (heap snapshot), 8.57 GB the
renderer at the first frame. Some cases (frame median, the opaque pass):

| case                                   | median   | opaque  |
| -------------------------------------- | -------- | ------- |
| MG whole, GPU culled, raster off       | 0.94 ms  | 0.64 ms |
| MG zoom 0.3, raster off                | 4.88 ms  | 4.49 ms |
| MG zoom 0.3, 2560×1440, raster off     | 17.5 ms  | 16.4 ms |
| MG whole rotating, raster off          | 2.05 ms  | 0.64 ms |
| MG zoom 0.3 rotating, raster auto      | 4.01 ms  | 0.34 ms |
| MG whole, culled on the CPU            | 13.2 ms  | 4.29 ms |
| PetWorld zoomed, rotating 10° a frame  | 1.73 ms  | 0.04 ms |

**Acceptance** at each step: no case's median or opaque pass more than 3% over the baseline beyond the runs' noise
(repeated where it is), the allocation guards unchanged, the visual suite unchanged, the yeast cell's load no slower;
the gains (heap, GPU storage, load) measured, not estimated. The vertex stages' work (a rotation by the quaternion for
the matrix product it replaces, two reads instead of three) is what the GPU cases measure.

## Decisions open

- **The carrier**: poses in cells of their own beside `aTransform` (empty then), read through `getTransformParts`,
  `getInstanceTransforms` and `getInstanceMatrices` (step 2, a smaller change), or `aTransform` holding
  `InstanceTransforms` so the compiler finds every reader: to revisit, a reader missed shows only at run time.

- **Scale an instance or a target**: an instance's (32 B, `scaleByRadius` per particle) or a target's (28 B a pose,
  the scale a segment constant): by which targets scale by radius per particle (shape targets do).
- **WebGL2 encoding**: instanced attributes or `tInstanceData` texels.
- **Particles not a range** (a target's particles scattered in the list): a particle index an instance (4 B), or the
  list sorted by target on load (the pack's order) so every target is a range.
