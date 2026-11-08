# Plan: a transient object transform on the GPU

## Status

Built: steps 1 to 4, then one transform per merged object (below). Next: a rigid release (below). Not done: the
heavy timing test and the frame time at rest (Perf below).

**As built**, where it differs from the plan below:

- `drawnBoundingSphere` is not a cell: `getDrawnBoundingSphere` (`renderable/util.ts`) derives it where read, cached by
  the versions of `boundingSphere` and `uObjectTransform`, `boundingSphere` itself at rest. Read by the scene's bounds,
  the renderer's visibility test and the raster's scissor, the draw table's rows, the clip caps, the CPU culling without
  a grid, merged objects' segment spheres, the culler's records. The camera's reset check keeps `boundingSphere` (a
  drag would reset the camera otherwise).
- CPU culling: `getCullGrid` (`renderable.ts`), a copy of the grid with its cell and batch spheres moved, a new object
  per transform (the clip classes, cached by the grid, classified again).
- A merged object moves as one (below), not by segment: no rows in `tSegment` or the culler's segment data.
- The particle target visual keeps its rebuild for `transform` (a placement); the preview uses `previewTransform`.

Tests: `cull-segment.test.ts` (a moved grid culls as its cells moved), `merged.test.ts` (a merged object's own
transform moves its drawn bounds, nothing else; a member's is not read), Prism's particles test (the preview moves the
drawn bounds, no instance touched), native WebGPU (`webgpu-preview-transform.test.ts`: mesh and impostors drawn as the
same placement; `webgpu-gpu-culling.test.ts`: a merged object moved drawn as the move baked, culled on the GPU and the
CPU; `webgpu-sphere-raster-split.test.ts`: a merged object moved drawn by the raster as the move baked), the visual
tests (WebGL against the commit before: unchanged; WebGPU against WebGL).

## One transform per merged object

**Built.** A merged object moves as one: its own `uObjectTransform`, never its members'.

**Why.** Members of a merged object always share the transform in practice: merging is within one representation
(`mergeRenderObjectsByType`, from units representations and the particle targets), and `previewTransform` is set per
representation (the move handles set the same delta on every representation of the object). The first build moved
each member by rows in its segment data; in a drag's profile of a mesoscale model (a 35 ms frame of `preview`) that
cost:

- 60 to 70% of the frame moving every culler cell on the CPU by its segment's transform and uploading them, O(cells);
- the merged `sync`: each member's transform version compared, its rows written, O(members);
- `Visual.setPreviewTransform` per visual (a few thousand targets), where only the drawn objects need it;
- on the GPU: three `tSegment` fetches per vertex while moved; merged objects left the sphere raster for the drag.

**What.**

- **The API**: the transient transform is a property of a drawn render object. A representation sets it on what it
  draws (its `renderObjects`: merged objects and the members left unmerged), not on every visual;
  `mergeRenderObjectsByType` takes the representation's `previewTransform` and sets it on each object it returns, so a
  regrouping mid-drag keeps it. A member's own `uObjectTransform` is not read while merged: no member can disagree.
- **Merged values**: `uObjectTransform` is the merged object's own cell, drawn as an unmerged object's;
  `boundingSphere` the members' union as placed. A transform change does not reach the merged `sync`. The CPU culling
  moves each member's grid by the object's transform (`memberCull`), its still-frame check and the marked lists' key
  read the transform's version.
- **Shaders**: the vertex stages read `uObjectTransform` (GLSL; `uObject.uObjectTransform` in WGSL), merged or not;
  `tSegment` two vec4s a segment.
- **GPU culler**: the cull uniforms `uObjectTransform` and its largest scale `uObjectScale`; the cell and copy kernels
  move each sphere by them (`moved` in `cell-cull.wgsl`) before the frustum, LOD, clip, Hi-Z and raster tests. A
  transform change writes the two uniforms and culls again, never `syncLayout`. The records' prediction moves the
  segments' bounds by it and scales the raster's largest sphere.
- **Sphere raster**: merged objects stay with the raster while moved (`uParams.uObjectTransform`, the object's own).
- **CPU-culled merged objects** (small by being CPU-culled): every segment's grid moved by the one transform
  (`getCullGrid`). Moving the camera's planes into the object's space instead (O(1)) only if it shows in a profile.

A drag's frame on the CPU: the representations' `setState`, the scene's bounds and the uniform writes, O(drawn
objects); nothing per cell, member or instance.

**Merging across representations** (proposed in `many-culled-objects.md`, direction B) would break the rule: two
structures merged together move apart under the handles. Proposed for it: the group key takes the transform's owner
(what the handles move as one: a structure or a particle list, with all its representations), so a merged object never
spans two owners. Per-segment rows (in `tSegment` and the culler's segment data, as first built) are the fallback if
that proves too coarse.

## Next: a rigid release moves what is built

**Built.** A release of the MG cellpack model (2.1M particles, 848 targets) from 1.1 s to about 560 ms after steps 1
and 2 (the browser); then step 3 and the hot loops inlined, in Node on Dawn (`--cpu-prof`-style, inspector sampling)
from about 800 to 310-490 ms. The release of a drag
commits the move to the data (`TransformParticleList`, `TransformStructureConformation`) and everything downstream is
built again as from new data, though a rigid move changes no topology. In a release of the MG cellpack model (about
1.1 s):

- the data: `ParticleList.transform` copies the list and resets its property data, so the hierarchy index, the target
  groups and the hierarchy colors are computed again (about 160 ms), the transforms and the copy (about 140 ms); the
  representation's signatures for precompiling again (`getParticleTargetSignatures`);
- the target visuals: `createTargetParticleTransform` (about 400 ms): the instance grid built again
  (`calcInstanceGrid`, 220 ms), the bounds (extrema, about 75 ms), the matrices; colors and sizes again (a new list
  counts as changed particles);
- the merged culler: a member's grid and `aTransform` change the layout key, every cell split again (`splitCell`,
  `setSegmentScales`, about 130 ms);
- the merged copy: `rebuildInstances`, `packAffineRows` (about 65 ms).

Position independent or moved rigidly: the hierarchy, the groups, the colors, the signatures, which instances share a
cell and the culler's split of them. Only the instances' matrices, the cells' spheres and the bounds move.

**What.** A placement change known to be a move by `M` (the decorator's new matrix times the inverse of its old) moves
what is built:

1. **The culler** keeps its layout for a moved grid. `moveInstanceGrid(grid, M)` (core, `instance-grid.ts`): the
   offsets, the instance order and the batches' cells shared, the cell and batch spheres and `cellTransform` moved;
   the grid records its layout (`layout`: the grid first built, shared by every grid moved from it) and its placement
   from it (`placement`). A member whose new grid has the layout the culler split, every other cell of the layout key
   the same but the grid, `aTransform` and the bounds: its cells' spheres moved by
   `placement · placement_split⁻¹`, its segment scales by the move's scale; no `splitCell`, the records kept. The
   instances are uploaded again (O(N)). As built: the culler keeps per member the grid, layout, placement, transforms'
   and bounds' versions it split, per cell its segment and sphere; the layout key is the fixed part (counts, invariant
   spheres, draw counts, levels) and the moved part (grids, transforms, bounds); a grid as split with its transforms or
   bounds changed is not a move. `MergedCullerCounts` counts splits and moves (tests).
2. **The target visuals** move instead of rebuilding where the new list is the old one moved by `M` (the same targets
   and indices): `aTransform` and the grid by `moveInstanceGrid`, the bounds moved (no extrema); the colors and sizes
   kept where the theme does not read positions. As built: `Particle.transform` records the move (`movedFrom`: the
   first list of a line of moves and the moves composed; `withTargets` keeps the line), `Particle.getMove` gives the
   move between two lists of a line; `moveTransformData` (`transform-data.ts`) moves `transform`, the grid and
   `aTransform` where no representation transform applies around them (else the rebuild); the visual keeps its
   index set where the groups computed again are equal, rebuilds colors and sizes of a vertex or volume granularity
   only.
3. **The data**: `Particle.transform` records its source and `M` and keeps the property data that does not depend
   on position (the hierarchy index, the target groups); the signatures cached by the targets (`targetMapping`) and the
   props. As built: `transform` keeps the target groups, the fiber mask and (without a scale) the spacing; the
   hierarchy colors' index is cached by the compartment array (and the entities), which a moved list shares; the
   target representation's signatures cached by its target mapping, groups, graphics and props (80 ms for 848
   targets). Then: `Particle.transform`'s loops and `moveTransformData`'s products inlined (`moveAffine`, in place;
   `moveTransforms` a column at a time). Left: the O(N) products and the instances' upload (`packAffineRows`,
   `rebuildInstances`), the first frame drawn.

What remains: the O(N) copies and uploads of the matrices (about 50-100 ms there).

**Not**: keeping the placement on the GPU (the decorator's matrix a lasting object transform, the representations on
the list as loaded): a release O(1), but every CPU reader of instances would apply it (export, themes by position,
the overlay's depth range, helpers, picked positions), against the transient scope below.

Tests: `moveInstanceGrid` (`instance-grid.test.ts`: the spheres and transforms moved, the cells and order shared, moves
composed); the culler keeps its split for moved grids and draws as members built there (`webgpu-gpu-culling.test.ts`);
a list's move recorded and given between two lists of a line (`particle-list.test.ts`); a list moved through
`particles.transform` (twice) leaves its target's grid moved from the one first built, its bounds and cells moved
(Prism's particles test); a list moved keeps its groups and spacing, its transforms computed again
(`particle-list.test.ts`).

## Scope: transient only

This is for transient visualization only: an object shown moved while something is under way, such as the drag of
the move handles before the release. It is not a way to place objects:

- where an object is lives in its data (the transform decorators, `TransformStructureConformation`,
  `TransformParticleList`) and in `Representation.State.transform`, baked into the instances as now;
- the transient transform is the identity at rest; whoever sets it puts it back (the move handles on release, once the
  state has the new placement);
- what reads instances on the CPU (export, themes by position, color smoothing, the overlay's depth range, helpers)
  ignores it: a transient view, not a placement.

## Why

`Representation.State.transform` (`TransformData.matrix`) is baked into the instances on the CPU: `updateTransformData`
writes `aTransform[i] = matrix · transform[i] · extraTransform[i]`, then rebuilds the instance grid from it;
`Visual.setTransform` then recomputes the bounding sphere from every instance. A merged object's GPU culler sees a
member's `aTransform` change as a new layout and splits every grid cell again (`syncLayout`, `splitCell`). The move
handles' preview sets it every frame, so a drag costs a full rebuild per frame, linear in the instances with large
constants:

- the profile of a drag of a mesoscale model: `createTargetParticleTransform` (the target visual rebuilt from the
  particles, the stopgap for the preview), `calcInstanceGrid`, `computeSphereExtrema` / `includePositionRadius` (the
  bounds), `splitCell` (the merged culler), then the uploads;
- structures take the same path (`Visual.setTransform`), cheap only because they have few instances.

That is right for a placement that stays (paid once), wrong for one that changes every frame.

## What

A second, transient transform: `Representation.State.previewTransform` (a `Mat4`, the identity by default), applied on
the GPU as an object uniform, `uObjectTransform` (`m4`), after everything baked:

- `aTransform`, the instance grid, the color grids and `boundingSphere` stay as baked (with `state.transform`);
- the world as drawn is `previewTransform · (baked)`: the vertex stages and the culling apply it;
- what culls and bounds the object while it is set is derived in O(1): `drawnBoundingSphere`, `boundingSphere`
  transformed by it (`Sphere3D.transform`: the largest scale for a scale or shear), and the grid's cell and batch
  spheres transformed where they are tested (cell counts, not instance counts).

Setting it writes one uniform and one sphere: no instance data, no grid, no culler layout, no upload but the uniform.

## Where it applies

### Vertex stages (WebGL and WebGPU)

- `uModel` is no place for it: one per scene group (`group.view`, never set: the identity) scaled by `camera.scale` (not
  1 only in XR), the same for every object of a pass. The scale stays outermost, so the transient transform is in Å.
- `assign_position` and the WGSL programs' equivalents: `model = uModel · objectTransform · applyTumble(aTransform)`.
  `objectTransform` is `uObjectTransform`, or the segment's (merged, below). No define: always applied, so a drag never
  switches programs; at rest the identity costs one matrix product per vertex (measure: below).
- Every program reading `aTransform`: mesh, spheres, cylinders, lines, points, text, image, direct volume (its
  `vTransform` carries the transient transform, so `inverse4(vTransform)` in the fragment stage stays right).
- Instance clipping (`clip-instance.glsl`, `clip.wgsl` `clipCenter`): the instance's center as drawn.
- Grid lookups in the vertex stage (`assign-color-varying`: the color grid by position, overpaint, emissive, substance,
  transparency volumes): the position before `objectTransform`, as the grids are built from `aTransform`. A local in
  the vertex stage, not a new varying.
- `vModelPosition` is the position as drawn (pixel clipping, LOD fade, camera-plane distances). Bump noise (`fbm` of
  `vModelPosition`) slides over a moving object during a drag: acceptable for a transient view.

### CPU culling (`renderable.ts` `cull`, `cullSegment`)

The cell and batch spheres transformed by the object's (or segment's) transient transform before the frustum, LOD,
clip and Hi-Z tests, skipped for the identity. The gridless path tests `drawnBoundingSphere`.

### Bounds

- `boundingSphere` stays the baked bounds; `drawnBoundingSphere` is it transformed, the same cell while the transient
  transform is the identity (no copy at rest).
- What the scene bounds and orders by while drawing (the scene's bounds, the camera's near and far, the transparency
  order, the gridless culling) reads `drawnBoundingSphere`. What fits the camera to the data, picks and focuses reads
  `boundingSphere` as now.

### Merged objects (`merged.ts`)

A segment's transient transform is three more rows in `tSegment` (two texels a segment now, five then; identity where
unset): `objectTransform` reads them where `dSegmented`, as the invariant sphere is read now (`segment.wgsl`,
`common-vert-params.glsl`). `merged.ts` watches each member's `uObjectTransform` version: its rows into `tSegment` and
the union of the drawn bounds again, both O(1) per member. A uniform flag (`uSegmentTransforms`: any member not the
identity) skips the three fetches at rest.

Alternative, for the GPU-culled path only: the instance copy kernel writes `M · T` into the compacted transforms (free
on the GPU), and the draw skips the segment's matrix. Two meanings of `aTransform` for one program; only if the fetches
measure.

### GPU culling (`merged-culler.ts`, `cell-cull.ts`) and the sphere raster

- `segmentData` (12 floats a segment) gets the segment's three rows and its largest scale; the cell kernel transforms
  each cell's sphere by its segment's, the copy kernel each instance's center (frustum, LOD, clip, Hi-Z, `rasterized`).
- The layout key (`getLayoutKey`) does not take the transient transform: a change writes `segmentData` only, never
  `syncLayout`. Gridless members' cells are their baked bounds, transformed in the kernel like the others.
- Sphere raster: unmerged objects' `modelView` times the transient transform (`uInvModelView` with it); the culled path
  reads the segment's rows with its other segment data.

### Representations and their users

- `Representation.State.previewTransform`, through `updateState` and every representation's `setState` (structure
  units, complex, volume, shape, particles, the particle targets) to `Visual.setPreviewTransform`: the uniform and
  `drawnBoundingSphere`.
- Prism's move handles (`sources/move/gizmos.ts`) preview with `previewTransform`, not `transform`, and put
  it back to the identity once the state has the new placement.
- The particle target visual (`target/common.ts`): the stopgap `setTransform` (a rebuild from the particles per frame)
  goes; it takes `previewTransform` as any visual. Its bounds after an in-place update of the particles (computed from
  the transforms) stay.
- `state.transform`, `Visual.setTransform` and their users (volume streaming, plugin code placing a representation)
  stay as they are.

## Steps

1. **Unmerged objects.** `previewTransform` in the representation state, the uniform and its schema entry (the object
   uniform block, both backends); every vertex stage and instance clipping; grid lookups before it;
   `drawnBoundingSphere`; CPU culling with transformed spheres; the sphere raster's unmerged path. Merged objects with a
   member not the identity unmerge until step 2 (or bake it for them: whichever is less code).
2. **Merged, CPU-culled.** `tSegment` rows, the member watch, `cullSegment` by the segment's transform.
3. **Merged, GPU-culled.** `segmentData` rows, the cell and copy kernels, the layout key, the raster's culled path.
4. **The move handles and the target visual** on `previewTransform`; the stopgap goes.

## Tests

- Quick: setting `previewTransform` leaves `aTransform`, `instanceGrid`, `aInstance` and `boundingSphere` at their
  versions and moves `drawnBoundingSphere`; back to the identity, the drawn bounds are the baked ones again;
  `cull-segment.test.ts`: an object with a transient transform culls as the same object with it baked (the visible
  cells, at each LOD level, clipped, occluded); the uniform layout tests with the new member; a merged member's
  transient transform changes `tSegment` and not the instances.
- Medium: the merged culler, a member moved: no `syncLayout`, the culled set as baked; the native GL suite; Prism's
  move test: a drag's frames leave the instances' versions alone.
- Visual: per program, an object with a transient transform against the same transform baked (both backends): mesh,
  spheres, cylinders, lines, points, text, image, direct volume; merged and unmerged; GPU-culled; with a smoothed
  color, an overpaint volume, instance clipping, LOD levels, the sphere raster. This reaches every program: the full
  suite.
- Heavy: setting the transient transform of an object of a million instances within a bound far below the rebuild's
  (the rebuild timed next to it).
- Perf: frame time at rest on the large scenes (PetWorld, the mesoscale models) within noise of before: the extra
  product per vertex, the segment fetches behind their flag.

## Open

- The merged draw's segment rows: gone with one transform per merged object.
- Non-rigid matrices (a scale or shear): the drawn bounds take the largest scale; impostors (spheres, cylinders)
  assume a similarity, as they do with it baked.
- `hasReflection` comes from the baked transforms only: a mirroring transient transform would flip faces wrongly for
  its duration; the move handles never mirror.
