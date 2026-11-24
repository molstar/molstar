# Plan: Every units representation merged, on one segment core

## Status

- **Proposed**: units representations keep one set of segments (a record per symmetry group, its part's shared values)
  and draw every part merged: what does not merge today gains the segment data it lacks, so mergeability stops being a
  question an update answers. Follows
  [merged-from-units.md](merged-from-units.md) (the merged visual over records) and the review of the merged path that
  found the issues below.
- **Fixed already, kept as they are**: a merged object's bounds leave out members that draw nothing; unit transforms
  found by unit for the representations' own groups (`StructureUnitTransforms.getSymmetryGroupTransforms`); the clip
  write skipped by the clip cell, not by the values object (records and their shared values one entry); the shifted
  cull lists the same object per source list (WebGPU's draw blocks a list object); a draw version per cull set on one
  counter; a new merged visual published once its first update merged; the per-group visuals' retry of merging judged
  without their layers' types and keyed by the theme too; palettes merged (step 4's first item); a marking redraw
  after a side view's GPU culling (a ray's pick) culls its frame's view again first (`Renderer.isGpuCullStale`); the
  CPU culling's key holds the clip objects (a clip moved alone culls again where a culler culls for the raster alone).

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

- **One visual over the groups** (`UnitsMergedVisual` generalized: segments, parts, records) for every units visual
  with a builder, whatever its group count; a part a geometry type and merge key, so a factory choosing by the data
  (`selectVisual`: GPU or CPU surfaces by a unit's size) gives a part of each type, not visuals per group. Why: one
  place that holds geometry, theme data, markers and layers per group; no switch between visuals.
- **Every part merged.** What `canMergeValues` rejects today, and what it takes:
  - **Palettes** (`dUsePalette`, done): no built-in structure color theme sets one (only `volume-value`, a volume
    theme); the palette, its domain and default color are the theme's, so shared by all segments (the merged values
    spread the shared ones), and the shaders look the palette up after reading `tColor`, which the merged path reads
    by segment.
  - **`vertexInstance` colors and layers** (per vertex and instance): a segment-aware index, its first texel plus its
    local instance times its vertex count plus its local vertex (two values more a segment, as `aSegment` holds the
    group-instance offset), in GLSL and WGSL. The texture split per segment (`getTexelCount`) a vertex-instance count.
  - **Smoothed colors and layers** (`smoothColors`, a theme's `preferSmoothing`; layers smoothed on surfaces): with a
    context `volume`/`volumeInstance` grid textures, a 3D grid in 2D tiles per object with its own transform and size;
    without one `vertex`/`vertexInstance`, interpolated on the CPU. Merged either as the vertex types (`vertex` merges
    today; CPU work and four bytes a vertex and instance, against the grid's smaller texture) or as an atlas of the
    segments' grids with each segment's grid offset, size and transform in its segment data. The vertex types first
    (they need only the index above), the atlas where their memory matters (many instances of large surfaces).
  - **GPU surfaces** (`texture-mesh`, `direct` colors): geometry and colors written on the GPU, per object textures
    on WebGL. On WebGPU the surfaces' vertices are chunks of one arena already (`gpu-surface-buffers.md`), which
    segments could address; on WebGL they stay plain objects, a render object a segment: the only plain kind left.
- **`renderMerged` retired** for units representations: the small counts measured no worse merged
  (`merging/units-{2,5,10,30}-{spacefill,cartoon}-*`, WebGL2 and WebGPU, RTX 4070 Ti): frames still, hovered, rotating,
  recolored, clipped and with a prop changed within the noise (0.3 to 1.2 ms either way, a rerun's order flipping);
  updates at 30 units cheaper merged (a prop 0.10 against 0.44 ms, a clip 0.09 against 0.20); setup alike (2 to 13 ms).
  `auto`'s threshold of 30 came from per-group merging, before records. The JS heap a scene retained is too noisy at
  these sizes to compare (0 to 6.5 MB for one scene).
  `unmergeableAt`, `mergesAgain`, the probe and `mergeRenderObjectsByType` for units go.
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
  view (`prepareCullView`) kept apart from the frame's, so a pass of the frame after it needs no culling again (now
  `Renderer.isGpuCullStale` and a new `prepareCull`, the marking redraw's). Then the question whether one segmented renderable with the
  culler chosen in `createRenderable` replaces `MergedRenderable`.

## Steps

1. **The key table** and parts owning their shared values; segments re-keyed by their key cells; member views from the
   table. Tests: a part outliving its first segment keeps following the size theme; a segment's levels changed moves it;
   the export of a part of groups with different unit counts. Measure: `repr/units-100k-spacefill` (heap a unit, build
   time) no worse.
2. **One merged object form**: the members' producer on records; the consumers of `members` moved to the views.
3. **Updates by the update state, transactional**: the fast paths folded in; a failed or aborted update leaves records,
   versions and parts as they were (tests for a throw and an abort mid-update, a backend switch alone).
4. **Every part merged**: palettes; the vertex-instance index (colors, sizes, layers); smoothing as the vertex types
   for merged parts; parts by geometry type (`selectVisual`), GPU surfaces as plain objects on WebGL; one group merged
   as well. Then the representation's switch, the probe, the retry and `renderMerged` for units removed. Tests: images
   as unmerged for each case lifted (palette, `vertexInstance`, smoothed colors and layers), markers kept through a
   theme change that used to unmerge. Measure: the visual suite identical (`merging/*`, both backends), the small
   counts (`merging/units-{2,5,10,30}-*`) against unmerged, a smoothed surface's memory as vertex types against grids.
   The grid atlas only where that memory says.
5. **Culling without merging's coupling**: the CPU fallback for views not culled on the GPU, the side view's GPU state
   apart; then decide on one segmented renderable.

## Open questions

- Whether WebGPU's GPU surfaces merge through the vertex arena (segments of its chunks), or stay plain objects as on
  WebGL.
- Per-unit visibility (`setVisualState`): a segment left out of the draws (culling skips segments already), needed
  once no unit has an object of its own.
