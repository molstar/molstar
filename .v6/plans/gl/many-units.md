# Plan: 100k units, the costs per member

## Status

- **Proposed**: a representation of some 100k small units, merged ([merged-render-objects.md](merged-render-objects.md)),
  within the JS heap of a browser tab: the per-unit visuals' memory per member from about 41 KB to under 10 KB, most
  of it bounding-sphere extrema. Steps 1 and 2 done: about 21 KB a member, the capsid's spacefill in 9.1 s (17.5) and
  2.8 GB of heap (4.6). Steps 3 to 5 done: 100k units 6.7 → 3.8 s but 20.8 → 19.5 KB a unit, the capsid 6.7 s and
  2.5 GB. What stays (live heap, a member): its values object (about 130 entries, 2.3 KB) and their cells (1.8 KB),
  the visual's closures and props (2.2 + 0.6 KB), its transform data (1.2 KB), color, size, marker and the empty layers'
  cells and objects (about 2.5 KB), the merge bookkeeping (cached keys and member versions, about 1 KB), the unit's
  symmetry group and transform groups (model side, about 1.3 KB): spread over the member's objects, under 10 KB only
  with fewer of them ([merged-from-units.md](merged-from-units.md)). Found with an NPC capsid (`NPC-NUP98_Capsid.data`, LAMMPS, not in the repository):
  its two scaling bugs are fixed (merged layers written part by part with a scan of all parts each, 72 s;
  `Renderable.memberCells` observing every member cell, a hang), see
  [draw-table-observing.md](draw-table-observing.md) for the latter.

## Question

The capsid is 525k atoms in 94,105 units (one per molecule, about 6 atoms each). Spacefill with its per-unit visual
(`element-sphere`, one render object a unit), merged, WebGPU: the visuals and their render objects take 17.5 s and the
heap grows from 0.75 to 4.6 GB, about 41 KB a member; the scene commits in seconds and draws at vsync. 4.6 GB is at
the edge of a tab's heap: ball-and-stick (two objects a unit) does not fit, which is why representations switch to
their structure-level visuals above 10,000 unit groups (`structure-element-sphere`; ball-and-stick so: 3.6 s, 2.3 GB).

Live heap by origin (sampled heap profile, all members):

| Where                                                                                                    | Heap    | A member |
| -------------------------------------------------------------------------------------------------------- | ------- | -------- |
| bounding-sphere extrema: up to 98 `Vec3` a sphere, a unit's sphere and those derived from it            | ~1.9 GB | ~20 KB   |
| absent layers: a null texture a grid layer (`createNullTexture`), grid uniforms, the empty layers' cells | ~0.4 GB | ~4 KB    |
| instance transforms and grids (`createTransform`, `calcInstanceGrid`, `calcTopGrid`) for one instance    | ~0.3 GB | ~3 KB    |
| the visual's closures and state (`UnitsVisual`)                                                          | ~0.2 GB | ~2 KB    |
| value cells (about 129 an object)                                                                        | ~0.2 GB | ~2 KB    |
| merge keys (strings, `getMergeKey`)                                                                      | ~0.1 GB | ~1 KB    |

Time: about 7.5 s of the 17.5 s is garbage collection; merging takes about 6 s, 4 s of it the members' bounding
spheres (`updateBoundingSpheres`: the union includes every member's extrema, 49 directions each).

The extrema: a unit's boundary (`Unit.boundary`, `getBoundary` with atomic radii) keeps the helper's 98 extremal
points (49 directions, two each, `BoundaryHelper('98')`) however few the atoms: with radii they are points on the atoms' spheres,
distinct for 6 atoms. Without radii and at most 98 points `getBoundary` keeps the points themselves; the renderable's
invariant sphere (`calculateInvariantBoundingSphere`) does so up to 14 positions, otherwise keeps the helper's 98 with
duplicates. Each derived sphere copies them: the geometry's sphere (`Sphere3D.expand` of the unit's, `element.ts`),
the transformed one (`calculateTransformBoundingSphere`: `clone` and `transform`), the merged union's inputs. Extrema
matter for quality on larger point sets (protein chains): a sphere from extrema is tighter than from centers and radii.

## Proposal

1. **A case of many small units** (done: `repr/units-100k-spacefill`, retained heap a unit in its notes): a synthetic benchmark like `render/units-1000-spacefill` with 100k units of a few
   atoms (setup time and heap, not only frames), and the capsid by hand. Measured per step: heap after the
   representation, its build time, the merge time.
2. **Never more extrema than points** (done; `@molstar/core`, `boundary.ts`, `boundary-helper.ts`, `sphere3d.ts`):
   - a sphere's extrema are at most its distinct input points: `BoundaryHelper.getSphere` drops duplicates (the
     extrema are copies of input points, exact equality); point sets of at most as many points as extremal directions
     keep the points (as `getBoundary` does without radii), the same limit in `calculateInvariantBoundingSphere`
     (98 for the fine helper, not 14);
   - sets with radii and few points: the centers as extrema with one radius for all (`Sphere3D.extremaRadius`, the
     largest radius; absent: 0), each extremum the ball of that radius around it: at most the points, a union or box
     still covering the atoms' spheres. Every reader takes the radius: `BoundaryHelper.includeSphere` and
     `radiusSphere` (`includePositionRadius`, `radiusPositionRadius`), `Sphere3D.includes`, `expand` (the radius
     grown by `delta`, the centers kept, instead of shifting the extrema along the principal axes), `transform` (the
     radius scaled as the sphere's), `clone`, `copy`, `Box3D.fromSphere3D` and `addSphere3D` (from 14 extrema: the box
     grown by the radius), `structure/util/boundary.ts`, the assembly-symmetry axes' sizes, the debug helper's extrema spheres;
     `setExtrema` sets it (0 unless given), so that a reused sphere keeps no stale radius. A field of the public
     `Sphere3D`: a reader outside missing it under-sizes by up to the radius, recorded in
     `.v6/plans/breaking-v6-changes.md`. Rejected: the points pushed outward by their radius (exact only radially:
     up to a radius short across, e.g. a two-atom unit has no thickness), no extrema (never too small, but loose for
     elongated small units);
   - larger sets keep the directions' extrema (quality: chains), copies as now;
   - derived spheres keep the count and the radius (`expand`, `transform`, `clone` map one to one).
   - measured (`repr/units-100k-spacefill`, WebGPU): at 10k units 1538 → 445 ms and 43.8 → 23.0 KB a unit; at 100k
     6.7 s and 20.8 KB a unit (before: past the heap, not finished). Bounds only grew: units with hydrogens by the
     largest radius for all (ace2 split: 7.75 → 7.96 Å), padded spheres' boxes to the exact padding (the smoothing
     grids, about 0.4 Å); every unit sphere encloses its atoms with their radii (`structure-building.test.ts`).
3. **Absent layers shared** (done): one null texture and frozen grid vectors for every layer without a grid
   (`NoLayerGrid`, `layer-type.ts`): their writers replace them (`ValueCell.update`), geo-export only reads them.
4. **One instance, a one-cell grid** (done): not without a grid, as gridless objects draw without levels of detail
   (`cellSize <= 1` in `cullSegment`, merged LOD scales); `calcInstanceGrid` of one instance builds its cell and batch
   directly (`calcSingleInstanceGrid`: three arrays, the offsets shared), not by two grid lookups.
5. **Merge keys interned** (done): equal keys and define key lists shared as one string and one array (strings compare
   by value: no result changes; a hash could merge unequal members on a collision).
6. **Value cells**: the single observer without a `Set` ([draw-table-observing.md](draw-table-observing.md) step 5).
7. **Measured**: the capsid's spacefill merged under 2 GB of heap and 8 s; its ball-and-stick with the per-unit visuals
   (about 190k members) within a tab's heap, drawn and marked; `render/units-1000-spacefill` and the merged benchmarks
   unchanged per frame; bounding spheres: images of the culling and clipping cases identical, and the spheres of the
   existing structures (`Unit.boundary`, the geometry's) no larger than now by more than a small margin (a test over
   the test structures' units, both helpers).

## Open questions

- Whether the visuals' structure-level switch (above 10,000 unit groups) should move once a member costs under 10 KB:
  the per-unit visuals mark, pick and update by unit, the structure-level ones rebuild all.
- Whether a merged object needs its members as full render objects (values, cells, spheres each): the larger change, a
  merged object built from the units directly. Only if the steps above leave a member above some 10 KB.
- The visuals' own state (`UnitsVisual`, about 2 KB a member): measured after the steps above.
