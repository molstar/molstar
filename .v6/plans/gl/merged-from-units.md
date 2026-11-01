# Plan: Merged objects built from units, without member render objects

## Status

- **Proposed**: a units representation that merges builds its merged render object from the symmetry groups directly,
  through one visual over all of them, not from a visual, geometry values and render object per group that
  `createMergedValues` then joins. Follows [many-units.md](many-units.md): after its steps a member costs about 19.5 KB,
  spread over its own objects; this removes most of them.

## Question

A merged member today is a full `UnitsVisual` (closures, props, theme state, location iterator), its geometry, its
render values (about 130 entries, as many cells: defines, uniforms, layers, transform data, bounds) and its render
object, kept so that each can update, mark, pick and be drawn alone, which merged members never are. The merged object
reads about 30 of those cells per member (`getMemberCells`) and copies or references their arrays as segment parts.

Measured (live heap a member, the NPC capsid's per-unit spacefill, 94k members): values object 2.3 KB, its cells 1.8
KB, visual closures and props 2.8 KB, color, size, marker and empty layers 2.5 KB, transform data 1.2 KB, merge
bookkeeping 1 KB; the unit's own (model side) 1.3 KB. The geometry itself (6 spheres) is under 0.5 KB.

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

1. **Segment records**: per symmetry group the least the merged values need: geometry, transform data, invariant
   sphere, group count, location iterator (created when needed: theme data, marking), the data versions. Measured
   alone at 100k units (`repr/units-100k-spacefill`): the target is under 5 KB a unit.
2. **Merged values from records**: `createMergedValues` and its `sync` take records (a member's arrays and counts, the
   fields `getMemberCells` reads now) instead of members' render values; the merged renderable, its culling (CPU and
   GPU) and drawing unchanged. The member-render-object path keeps working through a record view of a member's values,
   so both run until the switch.
3. **A merged units visual**: one visual over all groups of a representation, for the units visuals that merge
   (spheres, cylinders, meshes, lines, points; `MergedTypes`), built from the existing builders (`UnitsVisualBuilder`):
   creates and updates the records, the merged object, its picking, marking and layers. First the element spheres
   (spacefill, ball-and-stick atoms), then the bonds (cylinders, lines), then the meshes (surfaces per unit merge only
   on the CPU path, `MergedTypes`).
4. **The representation chooses**: a units representation whose merging resolves on (`renderMerged`, `auto` above
   `AutoRenderMergedThreshold` units) uses the merged visual where its builder supports it, per-group visuals
   otherwise (unmerged, smoothed layers, palettes and the other cases `canMergeValues` rejects: a group not mergeable
   falls back to a visual of its own beside the merged one).
5. **The structure-level switch**: whether representations still switch to their structure-level visuals above
   10,000 unit groups (`getSpacefillParams`, `getBallAndStickParams`, …) once a merged unit costs a few KB: measured
   (build, heap, frame, hover) at 100k units both ways. They stay an option: drawn, a merged object of many
   single-instance units costs a draw per visible unit ([merged-segment-draws.md](merged-segment-draws.md): the
   capsid's per-unit spacefill about 50 ms of CPU and 1.3 ms of GPU a moving frame here, structure-level 0.9 and 1.5 ms
   by hand), and clipping by the original unit must keep working (the `instance` clip variant).
6. **Measured**: `repr/units-100k-spacefill` and a ball-and-stick twin under 5 KB a unit; the capsid's per-unit
   spacefill and ball-and-stick within a tab's heap; updates as the merged benchmarks measure them (hover marking, a
   recolor, a trajectory frame) no slower; the merged render cases unchanged per frame; images identical to the
   member-merged path (`merging/*`, the marking cases), both backends.

## Open questions

- What a group's location iterator costs kept (theme data and marking need it): created per use and dropped, or kept
  for the groups a loci reaches often.
- Whether a structure-level visual (one geometry for all units, no instancing) and a merged visual (one object of
  segments, instanced per group) should become one: the merged visual keeps per-unit culling, levels of detail and
  partial updates, which the structure-level visuals lack.
- The visuals' state per group (visibility of a unit, `setVisualState`): kept per segment (a draw range dropped) or
  not supported merged.
