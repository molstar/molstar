# Plan: Merged objects of many single-instance segments

## Status

- **Proposed**: a merged object of many segments of one instance each (unique units: a capsid, a cell) costs CPU and
  GPU time per segment. Step 1 (the cull view by one version) done; step 2 (item-range draws) measured and rejected,
  and so is flattening (below): both conflict with clipping by the original unit or add special cases for it. Step 3
  (culling by batches) open, with less to gain without step 2.
- **Constraint**: clipping by the original unit (the clip variant `instance`, against the variant `pixel`) must keep
  working for merged objects: a unit is clipped by its own sphere's center (`clip_instance_wgsl`, the segment's
  sphere and the instance's transform). Whatever joins segments' draws keeps a unit's identity per item for it, in
  every mode, without a path of its own for clipping.

## Question

The NPC capsid's per-unit spacefill, merged (94k segments of one instance and about 6 spheres each, WebGPU, rotating):
about 50 ms a frame of CPU time (70 before step 1) and, measured by hand, about 15 ms of GPU time; its structure-level
spacefill (one geometry, one instance) 0.9 ms of CPU and 1.5 ms of GPU. The structure-level visuals stay the choice for
such scenes (above 10,000 unit groups); this narrows the gap of the merged path, which keeps per-unit culling, levels,
updates and instancing.

Where a moving frame goes (CPU profile, 60 frames):

- **Recording the cull view**: the versions of six culling cells per member (565k) every frame the view changed
  (`createCullView.set`): about 18 ms. Step 1.
- **Culling per segment** (`cullSegment`): a frustum and level test and a draw list entry for each of 94k segments,
  about 0.4 µs each: about 36 ms.
- **A draw per visible segment and level**: an instanced draw covers a vertex range times an instance range, so
  segments of one instance each, with their own vertex ranges, cannot share one; tens of thousands of draws a pass,
  most of the GPU time. GPU culling does not help: it issues a draw per segment and level as well (`gpuCulling:
  'auto'` rightly leaves it off, under 30 instances a draw).

## Proposal

1. **The cull view by one version** (done): `Merged.cullVersion`, advanced by `sync` where what culling reads of the
   members changed (layout, transforms and grids, spheres, levels); the merged renderable's cull view compares it
   instead of the members' cells. 70 → 50 ms a moving frame (production build).
2. **Item-range draws** (rejected): where every segment has one instance, each item reads its instance's data by index
   (per item its instance, `tItemInstance`; per instance its transform and segment data, `tItemInstanceData`), so a
   draw is a range of items and consecutive visible segments join into one. Built for spheres on both backends: the
   capsid's GPU time 1.31 → 0.23 ms a frame (RTX 4080 class), the visual suite identical. As a define it is a program
   variant known only from the data; as an object uniform (a branch per vertex in every merged sphere program, both
   layers bound) it measured the same: `render/assembly-60-spacefill-rotating` (60 instances a segment, the attribute
   branch) WebGPU 0.777 ms before, 0.775 define, 0.810 and 0.777 uniform, WebGL2 0.778, 0.778, 0.748 and 0.791;
   `render/units-1000-spacefill` equal. Rejected for its complexity: a second way to read instance data in the sphere
   shaders of both backends (macros over the attributes, item reads, the WebGL2 instance texture), two layers, a
   uniform, item joins in the draw lists.
3. **Flattening single-instance segments** (rejected): members of one instance copied into one segment of one identity
   instance (positions transformed, groups offset by the member's group base), culled and picked per member by item
   range and group base, no shader changes. It loses the unit per item that per-instance clipping needs (one instance,
   one sphere for all units); keeping clipping would need either a fallback to instanced segments while the `instance`
   clip variant is active, a CPU clip test per member (inexact with the clipping layer's per-group flags), or a unit
   index per item read by the clip chunk: each a special case of clipping against not clipping.
4. **Culling by batches of segments** (open): the segments' spheres in a grid of batches, as the instance grid's
   batches (`calcTopGrid`); a batch outside the frustum or its levels skips its segments, only batches on the edge
   test per segment. Without joined draws a visible batch still emits a draw entry per segment, so this saves the
   tests of what is out of view, not the entries (most of a view of the whole capsid).
5. **Measured**: the capsid's per-unit spacefill (CPU and GPU time a moving and a still frame, draws a pass), against
   its structure-level spacefill; `merging/*` and `render/units-1000-spacefill` unchanged or faster; picking, marking
   and clipping cases identical, both backends.

## Open questions

- Whether joining draws of single-instance units belongs before merging: building such units' geometry together while
  keeping their identity per item for clipping, picking and marking (as the structure-level visuals do per element),
  see [merged-from-units.md](merged-from-units.md).
- The instances of one object may be scattered: batches and segment order by the segments' spheres, not by object.
- Whether GPU culling should cull by batch and segment for such objects.
