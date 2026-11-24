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

- **A merged visual a kind** (done for the kinds): a representation's groups by the kind of visual its factory
  selects for each (`getVisualKind`, the factory's `select`; one kind for a factory choosing by the props alone), a
  merged visual over a kind's groups where its builder's type merges, per-group visuals for the others (a GPU
  surface's texture mesh, a kind found not to merge, `unmergeableAt` and `mergesAgain` a kind each). A kind's groups the
  same array while they are (the merged visual's fast paths); the per-group visuals kept by group. Why: one place
  that holds geometry, theme data, markers and layers per group; a factory choosing by the data merged too (the
  Gaussian surface's CPU meshes beside its GPU ones). One group merged too (done): a kind's groups always in its
  merged visual.
- **Every part merged.** What `canMergeValues` rejects today, and what it takes:
  - **Palettes** (`dUsePalette`, done): no built-in structure color theme sets one (only `volume-value`, a volume
    theme); the palette, its domain and default color are the theme's, so shared by all segments (the merged values
    spread the shared ones), and the shaders look the palette up after reading `tColor`, which the merged path reads
    by segment.
  - **`vertexInstance` colors and layers** (done): the index of data per vertex and instance the segment's offset (its
    first texel less its first instance times its vertex count and its first vertex) plus the merged instance times
    its vertex count plus the merged vertex: both in a sixth segment element (`SegmentTexels`, high and low 16 bits,
    exact as floats), read where used (WGSL: a color of that type, the layers of it; GLSL: `vertexInstanceIndex`).
    Sizes of it not (none made; the sphere raster's index not segment-aware), nor a source's (its instances on the
    GPU). Tested on WebGPU: two meshes of other vertex and instance counts, colors and an overpaint layer of it, drawn
    merged as alone.
  - **Smoothed colors** (`smoothColors`, a theme's `preferSmoothing`; done): a merged visual's segments smooth on the
    CPU (`processValues` without a context) into data per vertex (`vertex`, or `vertexInstance` of a theme per
    instance), not a grid per object (`volume`/`volumeInstance`). Images as the unmerged grids' on the visual cases
    (`geometry/merged/molecular-surface-smooth*`: the diff noise, a pixel's largest an edge's). Measured
    (`merging/assembly-gaussian-smooth-*`, 60 copies of 7qpd colored by element-index, CPU meshes): frames alike or
    faster merged, the representation built in 87 against 55 ms (WebGL2) and 70 against 60 (WebGPU), its JS heap
    12.5 against 5.7 MB and 15.1 against 6.0 (the colors per vertex and copy). A grid atlas where that memory matters.
  - **Smoothed layers** (planned, step 4b): a layer's type is its part's (`overpaintType`, `dOverpaintType` and the
    others' are shared cells, `dOverpaintType` in the merge key), so a part's segments hold one type of each layer.
    Segments by their smoothed layers: a segment a smoothed layer touches (`smoothColors` on, a mesh's
    `meta.resolution`, `reachesAny` of its loci) is in a part of that layer's `vertexInstance` key, the others stay in
    their `groupInstance` part. The keys first, the writes after: `setState` decides each segment's layer types by the
    probe, moves the segments whose key changed (`regroupByKey`: their records kept, a key of no part a new one), then
    writes the layers, smoothing on the CPU per vertex and instance (no context: no grid, as the merged colors), every
    segment of a part writing the same types into its shared cells. A layer cleared or no longer touching: back to the
    plain part, parts of one key joined. Why: memory only where a layer is (vertex count times instances per touched
    segment and smoothed layer, as the smoothed colors of a theme per instance), not every segment of a part; a part
    more per combination of smoothed layers (in practice one). Not all or none per part: a layer on one chain of a
    60-copy assembly would give every segment vertex data of zeros. Also: a segment made later or whose geometry is
    made again (its vertex data the old geometry's) smoothed again as it is written; GPU surfaces (texture meshes) are
    per-group visuals and keep their grids.
  - **A segment's update of another merge key** (done): of a builder that writes its values by the segment
    (`processValues`: smoothing by a surface's resolution), a part's key after each segment's update against its first
    segment's: another (segments wrote their own keys into the shared values) makes the visual again, each segment in
    the part of its key.
  - **GPU surfaces** (`texture-mesh`, `direct` colors): geometry and colors written on the GPU, per object textures
    on WebGL. On WebGPU the surfaces' vertices are chunks of one arena already (`gpu-surface-buffers.md`), which
    segments could address; on WebGL they stay plain objects, a render object a segment: the only plain kind left.
- **`renderMerged` retired** (done, removed from the geometry params): units representations always merge; the
  particle targets merge two or more objects of a type, or one culled on the GPU (`resolveMergedCulling`). The small
  counts measured no worse merged
  (`merging/units-{2,5,10,30}-{spacefill,cartoon}-*`, WebGL2 and WebGPU, RTX 4070 Ti): frames still, hovered, rotating,
  recolored, clipped and with a prop changed within the noise (0.3 to 1.2 ms either way, a rerun's order flipping);
  updates at 30 units cheaper merged (a prop 0.10 against 0.44 ms, a clip 0.09 against 0.20); setup alike (2 to 13 ms).
  `auto`'s threshold of 30 came from per-group merging, before records. The JS heap a scene retained is too noisy at
  these sizes to compare (0 to 6.5 MB for one scene).
  The representation no longer merges its per-group objects by type (`mergeRenderObjectsByType` the particle targets'
  alone). What merging everything exposed, fixed with it: a merged object's instances split for the sphere raster by
  their whole spheres, one through the near plane drawn as impostors, and `'only'` skips the impostors of the objects
  it draws whole alone (`isDrawnWhole`); the overlay's screen-sized planes read a merged object's transforms and
  segment spheres; a new context keeps the merged objects unless they hold foreign resources, Canvas3D's image waiting
  on the representations' updates; records written before a scene adds the object (a layer by `setState`) synced as
  the renderable is made, else WebGPU's program lacks the layer. Left: a merged object of one instance through the near
  plane never rasterized; layers on merged surfaces unsmoothed (above).
- **One table of segment keys** a type (done, `getMergeSegmentKeys`): a segment owns its record's keys (what the merged
  object reads, `MergeRecordKeys`) and `MergeSegmentKeys` (what differs by segment that the merged object has its own
  of: `uInstanceCount`, `uInvariantBoundingSphere`, a mesh's `meta`, the transform cells), its part's shared values
  the rest. Found by building per-group visuals of groups of other sizes and instance counts and comparing every
  cell: no other cell differs by group for the merged types. Debug mode tells a shared cell its segments write with
  different values in one update (`checkSharedWrites`, once a key): a cell missing from the table. Why: the
  per-segment state stated once and checked, not inferred from where writers happen to write.
- **Parts' shared values follow a live segment** (done): a part keeps the values its first segment was made with as
  its shared values (no copy a part), their segment-owned cells re-pointed to its first segment where that is another
  (`bindShared`: the one gone, or moved), so what is read of a segment through them (a layer's granularity,
  `sizeType`, the merge key's levels) is a live one's and a gone segment's geometry is not held. An update that changes
  segments' levels puts every segment in the part of its key (`regroupByKey`, `withMergeLodKey`): its record the same
  (cells, markers, layers kept) over another part's shared values; parts of one key become one; a key of no part a new
  part, its shared values a copy (`createEmptyMergeMember`). Measured: `repr/units-100k-spacefill` (WebGPU) 8.06 KB a
  unit against 7.91 (the two cells a segment owns more), its build 4.04 s alike.
- **Updates by the builders' update state** (done): a segment's update from the props, theme and structure it was last
  written with (its own: `Segment.props`, `theme`, `structure`), by its update state (`setUnitsUpdateState`: the
  context's compatibility, the builder's versions a group), nothing written where the state is none and nothing it was
  written with changed. The fast paths (nothing changed, the clip alone, the shared values alone) only where every
  segment is at the last update's props (it ended), in the same context; a builder with versions a group
  (`initUpdateState`: the cartoon's secondary structure, the interactions) reads each segment's update state on a
  props change too, written where it says. Measured (`merging/units-1000-*-merged-*`, WebGPU): spacefill alike; the
  cartoon's clip move and prop change 0.23 to 0.65 and 0.31 to 0.72 ms of update (the update state of 1000 segments,
  their versions read), its frames alike.
- **An update that does not end** (a throw, an abort; done): no record versioning lost (the records read again,
  `recordsWritten`), the next update without the fast paths, each segment brought from its own; a rebuild (a backend
  switch, another instance granularity) made beside the current segments and parts, which stay drawn and are disposed
  once it ends; the loci's structure and granularity as before. Not a rollback of the segments it wrote: they keep the
  new props, which the next update reads as theirs.
- **One merged object form** (done): records and shared values (`MergedGraphicsRenderObject.shared`, always), no
  `members`; `createMergedRenderObject(members)` takes its members' values as its records, the first's as the shared
  values. Member views (`getMergedMembers`: a record's own cells over the shared values; none of a source's) for the
  export and the debug helpers; the per-group merging's entry keeps its members' ids for picking (`pickMergedLoci`).
- **Culling not tied to merging** (last, larger): a view not culled on the GPU draws the CPU culling's lists of the
  identity region (as `merged-render-objects.md` describes the fallback) instead of nothing; the GPU state of a side
  view (`prepareCullView`) kept apart from the frame's, so a pass of the frame after it needs no culling again (now
  `Renderer.isGpuCullStale` and a new `prepareCull`, the marking redraw's). Then the question whether one segmented renderable with the
  culler chosen in `createRenderable` replaces `MergedRenderable`.

## Steps

1. **The key table** (done) and parts' shared values following a live segment; segments re-keyed by their levels;
   member views of the segments' own cells. Tested: a part outliving its first segment follows a new size theme as one
   made with it; a size theme's levels split and join parts as made; the member views of a part of groups with other
   instance counts each with its own; debug mode's check of shared writes.
2. **One merged object form** (done): the members' producer on records; the consumers of `members` moved to the views
   and the records (the per-group picking to the entry's ids, the particle targets' clip to `shared`).
3. **Updates by the update state** (done): the fast paths behind it; an update that does not end leaves every segment
   updatable from what it holds and a rebuild's current objects drawn. Tested: a throw mid-update, then the props it
   was made with again, drawn as made; a rebuild that throws keeps its objects; a backend switch alone checks each
   group against the new context; a builder's version moved with nothing else changed.
4. **Every part merged** (done but smoothed layers): palettes; the vertex-instance index (colors, sizes, layers); smoothing as the vertex types
   for merged parts; parts by geometry type (`selectVisual`), GPU surfaces as plain objects on WebGL; one group merged
   as well. Then the representation's switch, the probe, the retry and `renderMerged` for units removed. Tests: images
   as unmerged for each case lifted (palette, `vertexInstance`, smoothed colors and layers), markers kept through a
   theme change that used to unmerge. Measure: the visual suite identical (`merging/*`, both backends), the small
   counts (`merging/units-{2,5,10,30}-*`) against unmerged, a smoothed surface's memory as vertex types against grids.
   The grid atlas only where that memory says.
4b. **Smoothed layers on merged parts**: segments by their smoothed layers' keys as above, in `setState`, in
    `applyStateTo` for a new segment and after a segment's geometry is made again. Tests: a layer on one chain of a
    two-chain Gaussian surface splits its part in two and draws as the per-group visuals' smoothed layers (the image
    within the smoothed colors' noise), clearing it joins the parts; markers, other layers and the records kept across
    the split; the four layer kinds at once; `marking.test.ts`'s surface snapshots back to `vertexInstance`; the visual
    cases `geometry/mesh/gaussian-surface-cpu-smooth-layers` smoothed on both backends and a merged case against the
    per-group grids (`geometry/merged/*-smooth-layers`). Measure: `setState` time and the JS heap of a layer on one
    chain and on all of a 60-copy assembly surface, merged against the per-group grids.
5. **Culling without merging's coupling**: the CPU fallback for views not culled on the GPU, the side view's GPU state
   apart; then decide on one segmented renderable.

## Open questions

- Whether a layer on most segments of a large assembly (vertex data per instance) needs a cap or the grid atlas: decided
  by step 4b's measurement.

- Whether WebGPU's GPU surfaces merge through the vertex arena (segments of its chunks), or stay plain objects as on
  WebGL.
- Per-unit visibility (`setVisualState`): a segment left out of the draws (culling skips segments already), needed
  once no unit has an object of its own.
