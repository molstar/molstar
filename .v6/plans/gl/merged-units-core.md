# Plan: One segment core for merged and unmerged units

## Status

- **Proposed**: units representations keep one set of segments (a record per symmetry group, its part's shared values)
  whether or not their groups merge; merging decides how the parts are drawn, not which visuals exist. Follows
  [merged-from-units.md](merged-from-units.md) (the merged visual over records) and the review of the merged path that
  found the issues below.
- **Fixed already, kept as they are**: a merged object's bounds leave out members that draw nothing; unit transforms
  found by unit for the representations' own groups (`StructureUnitTransforms.getSymmetryGroupTransforms`); the clip
  write skipped by the clip cell, not by the values object (records and their shared values one entry); the shifted
  cull lists the same object per source list (WebGPU's draw blocks a list object); a draw version per cull set on one
  counter; a new merged visual published once its first update merged; the per-group visuals' retry of merging judged
  without their layers' types and keyed by the theme too.

## Why

The merged path grew beside the per-group one instead of under it. Five choices make most of its defects:

1. **Mergeability found out by building.** `UnitsRepresentation` makes a probe visual for its builder, builds a merged
   visual, destroys it where a part does not merge, falls back to per-group visuals, remembers the failure
   (`unmergeableAt`) and retries by the per-group visuals' values (`mergesAgain`). Each switch destroys one side and
   builds the other: the markers go (a selection gone from the canvas while the selection manager keeps it), the work
   is done twice, and the result depends on history. The fixes since each narrowed a symptom of this.
2. **Which cells are a segment's is implicit.** A record's prototype is its part's shared values ("the writers
   unchanged"), and which keys are its own is spread over about ten lists (`MergeRecordKeys`, `getSourceRecordFillKeys`,
   `getMemberCells`, the texture specs' granularities, ...). Cells that differ by segment but are not in a record
   (`uInstanceCount`, `meta`) are written through to the shared values, the last writer winning: the member views
   (`getMergedMembers`, read by the geometry export) read another segment's. A part's shared values are its first
   segment's values object and outlive that segment, so record-key cells read through them (`sizeType`, `lodLevels`)
   freeze. A segment's part is chosen by its merge key once, when it is made: a later change of its key cells (its
   levels under a size theme) is never checked again.
3. **The merged visual decides updates its own way.** Its fast paths (`unchanged`, `clipOnly`, `propsOnly` and
   `sharedOnly`, the early return in `updateObject`) duplicate what a `UnitsVisual` decides per update and have drifted
   from it: `isCompatible` with a new context and the builders' update-info versions are not checked where nothing else
   changed. An update that fails disposes what it made but not what it rewrote in place (records written, versions not
   advanced, `currentProps` not moved), and `reset()` before it is not undone.
4. **Three producers of `MergedGraphicsRenderObject`** (`createMergedRenderObject` from members,
   `createMergedRecordsRenderObject` from records, the particle source's), each giving `members`, `records` and
   `shared` other meanings: every consumer has to know which made an object.
5. **`MergedRenderable` is a second renderable with its own culling state.** GPU culling exists only for merged
   objects, and a view not culled on the GPU for one draws nothing (`'none'`, no CPU fallback), so every view must call
   `prepareCullView`, and a side view (a ray pick) overwrites the frame's GPU state that a marking redraw then misses.

## Design

- **One visual over the groups** (`UnitsMergedVisual` generalized: segments, parts, records) for every factory that
  chooses by the props alone (`choosesByProps`) with a builder; per-group `UnitsVisual`s only for factories choosing by
  the data (`selectVisual`) or without a builder. Why: one place that holds geometry, theme data, markers and layers per
  group, so mergeability can change without a visual changing.
- **A part drawn merged or plain.** A part whose values merge (`canMergeValues` on its shared values) is one merged
  object; one that does not is one plain render object a segment, its values the segment's record over the part's
  shared values with the cells that must then differ by segment (smoothed layers, palettes) the record's own. A flip
  re-packages the part: the same records, geometry and markers; no visual made or destroyed, no geometry built to find
  out. `unmergeableAt`, `mergesAgain` and the probe go.
- **One table of segment keys** a type (which keys a record owns, which its part shares, which only a plain part's
  record owns), from which `MergeRecordKeys`, the fill keys, the member cells and the texture granularities are
  derived. A writer's write of a record-owned key through the prototype is an error in debug mode. Why: the
  per-segment state is then stated once and checked, not inferred from where writers happen to write.
- **Parts own their shared values**: a values object of the part's (as `createEmptyMergeMember` makes for the paged
  pool), not its first segment's. A segment whose merge key changes (its key cells' versions, read in `sync`) moves to
  the part of its new key.
- **Updates by the builders' update state**, decided once for the props and theme (shared) and per segment for its
  group, as `UnitsVisual` does: the fast paths become outcomes of that state (nothing to do, the clip alone, the shared
  values alone), not checks beside it; `isCompatible` checked per update. **Transactional**: an update builds what it
  changes beside the current segments and commits at its end (records, versions, `currentProps`, the parts), so a throw
  or an abort leaves the visual as it was.
- **One merged object form**: records and shared values; `createMergedRenderObject(members)` makes records of its
  members' values; member views are built from the key table (a record's own cells, the rest the part's). The export and
  the debug helpers read through them.
- **Culling not tied to merging** (last, larger): a view not culled on the GPU draws the CPU culling's lists of the
  identity region (as `merged-render-objects.md` describes the fallback) instead of nothing; the GPU state of a side
  view (`prepareCullView`) kept apart from the frame's. Then the question whether one segmented renderable with the
  culler chosen in `createRenderable` replaces `MergedRenderable`.

## Steps

1. **The key table** and parts owning their shared values; segments re-keyed by their key cells; member views from the
   table. Tests: a part outliving its first segment keeps following the size theme; a segment's levels changed moves it;
   the export of a part of groups with different unit counts. Measure: `repr/units-100k-spacefill` (heap a unit, build
   time) no worse.
2. **One merged object form**: the members' producer on records; the consumers of `members` moved to the views.
3. **Updates by the update state, transactional**: the fast paths folded in; a failed or aborted update leaves records,
   versions and parts as they were (tests for a throw and an abort mid-update, a backend switch alone).
4. **Plain parts**: a part that does not merge drawn as plain objects of its records; the representation's switch, the
   probe and the retry removed; the per-group visuals kept for `selectVisual` factories. Tests: markers kept through a
   flip both ways; a smoothed layer, a palette and a theme change alone flip without a rebuild of geometry. Measure: the
   visual suite identical (`merging/*`, both backends), a flip's time against a rebuild's.
5. **Culling without merging's coupling**: the CPU fallback for views not culled on the GPU, the side view's GPU state
   apart; then decide on one segmented renderable.

## Open questions

- Whether a plain part's objects are worth keeping at all where the merged object could draw a part whose values do
  not merge (per-segment palettes and smoothed layers as segment data): fewer object kinds, more segment data.
- Direction B of [many-culled-objects.md](many-culled-objects.md) (merging across representations) would add a producer
  of merged objects: after step 2, on the one form.
