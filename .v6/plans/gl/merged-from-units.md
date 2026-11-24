# Plan: Merged objects built from units, without member render objects

## Status

- **Done** (steps 1 to 6): a units representation that merges builds its merged render object from the symmetry
  groups directly, through one visual over all of them, not from a visual, geometry values and render object per group
  that `createMergedValues` then joins. Follows [many-units.md](many-units.md): after its steps a member costs about
  19.5 KB, spread over its own objects; this removes most of them. What it left (mergeability found by building, the
  implicit per-segment keys, the merged visual's own update paths) continues in
  [merged-units-core.md](merged-units-core.md).
- **Step 2 done**: merged objects built from records and shared values; the particle pool on records (its first
  user).
- **Steps 1, 3 to 6 done**: `UnitsMergedVisual` (`repr/structure/units-merged-visual.ts`) used by
  `UnitsRepresentation` wherever `renderMerged` merges two or more groups and the factory chooses by the props alone
  (the element spheres, the bonds' cylinders and lines, the meshes on the CPU path); a part a merge key (step 4);
  per-group visuals where it does not merge. `repr/units-100k-spacefill` (WebGPU): 5.4 to 3.6 s, 19.6 to 7.9 KB a unit
  (the heap read after a second collection, see step 6); the marking and layer snapshots and the visual suite's images
  as the per-group visuals'. At 1,000 units (`merging/units-1000-*-merged*`) still, hover and rotating frames alike or
  faster, a recolor 3.8 to 2.1 ms, the scene's heap 31 to 13 MB.
- **What a unit holds** (a Chrome snapshot at 10k units, 8.2 KB): its record's cells (some 36 of its own, two objects
  each), its color, size, marker and empty layers' textures, two bounding spheres with their extrema (`Vec3` arrays),
  its geometry, the merged object's bookkeeping a member (its attributes' and textures' cells and versions in flat
  arrays; the index orders and marked positions made where used), the unit's own (model side). Records are fast objects
  from a template; the extra transforms shared until unit transforms are set.
- **Measured and not taken**: the empty layers' cells shared by the records until a layer is first set (0.2 KB a unit
  in Chrome, whose compressed pointers halve a cell; 0.9 in Node): a layer first set then replaces records' cells, so
  `setState` makes the merged object again and the canvas and the multi representations must follow a changed object
  list, new paths for a rare event.
- **The merged renderable's frames at 100k members** (profiled, `merging/units-*` at that count with the per-unit
  visual, by hand): a mark updates the object through its observed cells alone (`Scene.updateChanged`,
  `Renderable.update(observedOnly)`), so `sync` reads the members whose marker cells changed (two versions compared a
  member, not all of a member's cells), the culler's layout key is kept while `cullVersion` stands, the marker
  aggregates and the marked lists follow the members marked, and the CPU culling's lists stay; the culler's predictions
  read the segments' drawn spheres packed (refreshed with `boundsVersion`), and WebGL uploads a segmented layer's
  changed parts as WebGPU does. WebGPU: a hover 146 to 8.2 ms, a rotating frame 35.7 to 6.0. WebGL2: a hover 119 to
  39.6, a rotating frame 84 (its CPU cull of 100k segments and their draws, a uniform each). At 1,000 units no slower
  (a hover 1.1 to 0.8 ms). A clip moved (`merging/units-*-clip`): the parts' shared values' clip written once (the clip
  cells no record's), and the merged objects told no record was written (`recordsVersion`, `Merged.recordsUnchanged`:
  their sync as a mark's, not over every member) nor made again: 337 to 10.6 ms (WebGPU), 396 to 99.5 (WebGL2, its CPU
  culling's lists culled again for the clip). Other props changed alone (`merging/units-*-props`, the emissive): a part's first segment
  updated, the others not where that wrote its shared values alone (no cell of its record, nothing made again nor
  recolored; a size factor's levels and padding write every record): 300 to 17.4 ms (WebGPU), 337 to 84.5 (WebGL2).
- **Not pursued**: WebGL's frames at 100k members (the CPU cull a segment each, a draw each; segments coarser or
  drawn together would be structural changes for WebGL alone: WebGPU is the backend for such scenes); the 5 KB a unit
  target (7.9 KB is enough).

## Question

A merged member today is a full `UnitsVisual` (closures, props, theme state, location iterator), its geometry, its
render values (about 130 entries, as many cells: defines, uniforms, layers, transform data, bounds) and its render
object, kept so that each can update, mark, pick and be drawn alone, which merged members never are. The merged object
reads about 30 of those cells per member (`getMemberCells`) and copies or references their arrays as segment parts.

Measured (live heap a member, the NPC capsid's per-unit spacefill, 94k members): values object 2.3 KB, its cells 1.8
KB, visual closures and props 2.8 KB, color, size, marker and empty layers 2.5 KB, transform data 1.2 KB, merge
bookkeeping 1 KB; the unit's own (model side) 1.3 KB. The geometry itself (6 spheres) is under 0.5 KB.

A second user, the particle target's paged pool (`repr/particles/representation/target/paged.ts`): a member a target
of the list's catalog made up front (the yeast's 5,641), each a copy of all of a target member's cells (162 for plain
spheres, about 24 KB with its render object), never drawn, picked or updated alone, of whose cells the merged object
reads about 30. On records (step 2): 2.1 KB a member (38 keys, 17 cells its own), 132.8 to 11.5 MB for 5,641 plain
spheres members; the yeast's JS heap at the tiles benchmark's end 225.8 to 172-174 MiB (its target members lighter
than plain spheres').

What the per-group objects give that a merged object needs, and where it comes from without them:

- **Geometry**: `builder.createGeometry(ctx, unit, structure, theme, props)` per group, as now (the builder is the
  visual's, not the visual); kept as the segment's parts.
- **Instances**: the group's transforms (`createUnitsTransform`), the segment's instance range.
- **Bounds**: the group's invariant sphere (from the unit's boundary, `Sphere3D.expand` by the size), transformed per
  instance into the segment sphere.
- **Theme data**: color, size, marker and the layers per group, from the group's location iterator
  (`builder.createLocationIterator`), written into the segment's part of the merged data arrays.
- **Levels of detail**: per segment, as now (`aSegmentLod`, from the group's sphere count).
- **Picking**: instance → segment → group (`pickMergedLoci`), then `builder.getLoci` with the segment's group.
- **Marking and layers**: `builder.eachLocation(loci, structureGroup, apply)` per segment the loci reaches (by unit id,
  as `visualsByUnitId`), applied at the segment's group base in the merged marker and layer data.
- **Updates**: the visual's update state per prop and theme change (`builder.setUpdateState`), decided once for all
  segments (they share props and theme), then applied per segment.

## Proposal

1. **Segment records**: per symmetry group a record (`MergeRecord`) of cells its own: its geometry's (the type's
   vertex attributes and textures), counts, instances (`aTransform`, `aInstance`), bounds, grid, levels, its theme
   layers of a granularity (`tColor`, `tSize`, `tMarker`; the layers once set) and their aggregates. The rest the
   shared values' (one set: defines, uniforms, props, the clip, the layers of no granularity, the `*Type` cells). Beside
   the record, the segment: its group, its geometry (the builder's reuse argument), its transform's base and extra
   matrices (`transform`, `extraTransform`: the record's `aTransform` computed of them), its previous mark where one was
   made; the location iterator made per use (theme data, marking) and dropped. Measured at 100k units
   (`repr/units-100k-spacefill`, 22.6 KB a unit and 5.4 s now, WebGPU): the target under 5 KB a unit.
2. **Merged values from records** (done): `createMergedValues(type, shared, records, ...)` and its `sync` read
   members through records alone; the merged renderable, its culling (CPU and GPU) and drawing unchanged.
   - **A record** (`MergeRecord`, `merged.ts`): the cells of a member the merged object, its renderable and cullers
     read, so `sync` reads changes by their versions: counts, instances (`aTransform`, `aInstance`, `instancePoses`),
     bounds, grid (`instanceGrid`, `instanceCellSize`), levels (`lodLevels`, `padding`), sites (`uSiteOrigin`,
     `uSiteScale`), the markers' and layers' aggregates and types (`markerType`, `sizeType`: read per member, kept
     there rather than moved into the merge key), the type's vertex attributes, `elements` and the merged textures by
     key (`MergeRecordKeys`, `getMergeRecordKeys`). Typed so a read of another key does not compile; the quick and
     WebGPU suites passed with every member a picked record (`pickMergeRecord`: its keys alone).
   - **The shared values**: of the object's type and merge key (defines, uniforms, props, the clip, `uLod`, the layers
     of no granularity), an argument of its own: the cull's clip and level of detail, the second color's type and the
     dual color factor read from them, not from each member.
   - **The calls**: `createMergedRenderObject(members)` passes the members' values as their records (no copies) and the
     first's as the shared values; `createMergedRecordsRenderObject(type, shared, records, state, materialId)` takes
     records. Either gives one form, `Merged.records` and `shared` (no `members`: views of the records,
     `getMergedMembers`).
   - **A source's records** (`getSourceRecordFillKeys`, `fillMergeRecord`): what a fill writes its own cells
     (geometry, counts, bounds, levels, sites, aggregates, the layers of its own granularity), the rest the empty
     record's (what the source holds: instances, grid, markers; the layers of no granularity). A record's cells are
     resolved once by the merged object and its culler: a record is filled in place, not replaced.
   - **Measured** (the particle pool): see the question; the yeast's images at the tiles benchmark's still segments
     the same as before (outside 0 pixels, inside 113 of 1.66 M: the load's timing), its frames alike.
3. **A merged units visual** (`UnitsMergedVisual`, `repr/structure/units-merged-visual.ts`): one visual over all
   groups of a representation, for the units visuals that merge (`MergedTypes`), built from the per-group visuals'
   builders: creates and updates the records, the merged object (`createMergedRecordsRenderObject`), its picking,
   marking and layers. First the element spheres (spacefill, ball-and-stick atoms: impostors and meshes), then the
   bonds (cylinders, lines), then the meshes (surfaces per unit merge only on the CPU path). The contract:
   - **The builder**: a units visual carries its builder and geometry utils (`UnitsVisual.builder`, a reference); a
     factory choosing by the props alone (`chooseVisual`) or not at all gives one for all groups. One choosing by the
     data (`selectVisual`: GPU or CPU surfaces by a unit's size) gives a merged visual for each kind it selects that
     merges, per-group visuals for the others ([merged-units-core.md](merged-units-core.md)).
   - **The writers unchanged**: `createValues`, `createColors`, `createSizes`, `createMarkers`, `createUnitsTransform`,
     `updateValues` and `Visual.mark`, `setOverpaint` and the other layer setters write cells by name. A segment is
     written through its record, the shared values its prototype (one view set per segment was most of a recolor):
     the record's cells change, the shared ones are written alike by every segment (`updateIfChanged`). A segment's
     values made new (`createValues`) are picked to its record (a fast object from a template, its prototype set once)
     and dropped; the first segment's become the shared values, their segment-owned cells (`getMergeSegmentKeys`)
     following the part's first live segment ([merged-units-core.md](merged-units-core.md)). Who walks a record's keys reads its own alone
     (`Object.keys`). The layers' known-clear state and the previous mark keyed
     by the segment, not a render object.
   - **Updates**: the props' and theme's part of the update state decided once (they are shared), a group's part
     (`setUpdateState`, transforms, conformations, the hierarchy) per segment; the geometry made again (the builder's
     `createGeometry`, its `includesUnitKind` check) where its state says. The merged object follows the records' cell
     versions (`sync`), made again only where its layout or merge key changes.
   - **Picking**: the segment by the instance (`pickMergedLoci`'s search), then `builder.getLoci` with its group and
     the merged object's id. **Marking**: the segments a loci reaches by unit id (as `visualsByUnitId`), all for a
     loci of more units than segments; `lociApply` at the segment's group. **Layers**: per segment, never smoothed (a
     merged object keeps one layer type); clipping and wiggle alike. **State**: visibility, pickability, color only,
     theme strength, transform and preview transform the merged object's; unit transforms per segment (its extra
     matrices, its `aTransform`).
   - **Culling and merging**: on the GPU or CPU as `mergeRenderObjectsByType` decides for the members (its counts
     taken of the records: one function both call, `resolveMergedCulling`), per part; the merge check (`canMergeValues`)
     on each part's shared values, a segment in the part of its merge key (`getMergeKey` of its values made, step 4).
   - **Measured**: as step 6, against the per-group visuals merged (images identical, `merging/*`, both backends).
4. **The representation chooses** (done): a units representation whose merging resolves on (`renderMerged`, `auto`
   above `AutoRenderMergedThreshold` units, two groups or more) uses the merged visual where its factory makes one
   visual for all groups (`choosesByProps`) of a type that merges, per-group visuals otherwise (unmerged, a factory
   choosing by the data, palettes and the other cases `canMergeValues` rejects: all parts alike, their props and theme
   shared). The groups of another merge key are a part of their own (`Part`: its shared values, records and merged
   object; one or a few a representation), not visuals of their own: `auto` instance granularity resolves by a
   group's counts, and the first group made, whose values are shared, may be the one apart. Before, any group of
   another key left the whole representation unmerged. A new instance granularity makes the parts again. Whether the
   groups merge is found by building, not decided up front: a probe visual gives the builder, the merged visual's update
   ends at a part that does not merge (destroyed, the per-group visuals used), the failure kept by the props, theme and
   structure (`unmergeableAt`), and merging tried again where the per-group visuals' values merge after an update
   (`mergesAgain`, without their layers' types: smoothed there, not merged). Replaced by
   [merged-units-core.md](merged-units-core.md) step 4.
5. **The structure-level switch** (measured again after the frames' work, kept): representations still switch to
   their structure-level visuals above 10,000 unit groups (`getSpacefillParams`, `getBallAndStickParams`, …). At 100k
   units (`merging/units-*` at that count, WebGPU, by hand, a case a process), the structure-level spacefill against the
   per-unit one merged: built in 0.85 s against 5.4, the scene's heap 170 MB against 1.3 GB, a still frame 0.55 ms
   against 0.90, a hover 0.83 against 7.6 (215 before the frames' work), a rotating frame 0.59 against 4.7 (33.6), a
   clip moved 0.58 against 9.8, a prop changed 0.80 against 16.4. The per-unit visuals stay an option there (per-unit
   culling, levels of detail, the `instance` clip variant), now interactive on WebGPU.
6. **Measured** (done, without the 5 KB target: chased no further): `repr/units-100k-spacefill` 7.9 KB a unit
   (19.6 per-group merged) and 3.6 s; its ball-and-stick twin (element spheres and intra bonds, by hand) 13.6 KB a unit
   and 11 to 15 s, per-group merged not finished in 15 minutes. The merged path at 1,000 units as step 3's status
   says (hover, recolor, rotating no slower). Images identical to the per-group visuals merged on WebGL2 and WebGL1
   (the visual suite saved with the merged visual off, compared with it on: none changed) and WebGPU against WebGL2
   (all pass). The NPC capsid by hand in Prism (its per-unit spheres, 94k units): drawn, a backend switch and a clip drag
   fixed by it (its objects made with a disposed state; a clip moved through every segment).

## Open questions

- What a group's location iterator costs kept (theme data and marking need it): created per use and dropped, or kept
  for the groups a loci reaches often.
- Whether a structure-level visual (one geometry for all units, no instancing) and a merged visual (one object of
  segments, instanced per group) should become one: the merged visual keeps per-unit culling, levels of detail and
  partial updates, which the structure-level visuals lack.
- The visuals' state per group (visibility of a unit, `setVisualState`): kept per segment (a draw range dropped) or
  not supported merged.
