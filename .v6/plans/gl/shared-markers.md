# Plan: Shared marker data (one marker layer per structure group)

## Status

- **Current**: proposed, not started. Step 1 (measure) decides whether the rest is worth it.
- Related: [shared-textures.md](shared-textures.md) (one GPU copy of an image held by several render items, the model
  for the GPU side), [shared-draw-resources.md](shared-draw-resources.md) (theme layers in data pages, updated by
  range), [merged-render-objects.md](merged-render-objects.md) (members' layers as parts).

## Why

- Visuals of different representations over the same structure group often show the same atoms: ball-and-stick and
  a molecular surface of a ligand, spacefill and a gaussian surface. Each holds its own marker layer (`tMarker`,
  `createMarkers`), and a hover or select applies the loci to each (`lociApply`: finding the loci's groups in the
  visual), computes each one's status or average (`getMarkersAverage`, a pass over the layer when the status is mixed)
  and uploads each one's changed range.
- Visuals whose marker layers mean the same per index could hold one: the loci applied once, the average computed
  once, one upload, one copy in the data pages.
- What it saves is CPU per mark and per visual of a structure group shown more than once, plus uploads; memory hardly
  (1 byte per group and instance). Markers already upload by range ([shared-draw-resources.md](shared-draw-resources.md)
  step 3.2), and Canvas3D visits only the representations of the loci's root (`MarkIndex`). So the gain depends on how
  often groups are shown by more than one compatible visual, unmeasured: step 1.

## Design

1. **Compatibility**: two visuals' layers mean the same per index when they have the same location iterator kind over
   the same structure group (`StructureGroup`: the same `group` object, so the same units and elements) and the same
   instance granularity: index = instance × groupCount + group. Kinds today:
   - element (`ElementIterator.fromGroup`): element sphere, point, cross, ellipsoid, the gaussian, molecular and blob
     surfaces and wireframes, the gaussian density volume;
   - polymer element (`PolymerLocationIterator.fromGroup`): trace, tube, backbone, direction wedge;
   - nucleotide (`NucleotideLocationIterator.fromGroup`), bonds (`BondIterator.fromGroup`), and the structure-level
     element iterator (`ElementIterator.fromStructure`, per structure, not per group) each a kind of its own.

   Kinds never mix: a bond index is not an element index, and a cartoon (polymer elements) and a surface (elements)
   share nothing. The key is (kind, group object, granularity, `markerActions`), see 3.
2. **The shared data**: the whole `MarkerData` (`uMarker`, `tMarker`, `markerAverage`, `markerStatus`, `markerType`)
   by reference: the cells, not copies, so a mark through one visual is the state of all and the status uniform stays
   consistent. A registry per structure (`WeakMap` by group object, then by key) counts holders; a visual takes the
   data when its render object is created (`createMarkers` given the shared data) and releases it on destroy; the last
   holder frees it. A visual that joins a share sees its current marks (a new surface shows the selection at once;
   today it starts clear and the plugin selects it again, idempotent).
3. **Marks applied once**: today every visual applies the mark (`Visual.mark`). With shared data the first visual of a
   share applies it, the others skip it within the same mark: a share remembers the mark it applied (Canvas3D's mark
   sequence number with the loci and action). Not optional: `Toggle` (the plugin's `LociSelectManager.toggle`) applied
   twice undoes itself, and `getMarkersAverage` would run per holder. `markerActions` is part of the key: a
   representation that takes highlights only (validation reports, symmetry axes, measurement shapes) shares only with
   such.
4. **Marks aimed at one representation** (`Canvas3D.mark({ repr, loci })`, e.g. the mesoscale explorer's entity hover
   with `EveryLoci`): with shared data the other holders show it too. Either accepted (they show the same atoms), or a
   representation opts out (`state.shareMarkers: false`) where marks are meant per representation. Decision needed
   (Open).
5. **GPU, both backends** (as [shared-textures.md](shared-textures.md) §2 and §4): the layer's data-page range keyed by
   the `DataArray` object, counted by holders, written once per change (`TextureChanges` ranges, as now). Each layer
   already has its own word offset (WebGPU `getDataOffsetName`, WebGL object uniforms), so an item's `tMarker` offset
   may point into the shared range while its other layers stay in its own. The range must be on page 0 (or the item's
   page for that layer), else the item's `dDataLayerPages` define and so its program change. Byte counts divide by
   holders (shared-textures §5).
6. **Merged render objects**: members' layers are parts of the merged layer, written by the merged item
   ([merged-render-objects.md](merged-render-objects.md) step 4). Members of two merged objects sharing a marker array
   still write it into both merged layers: sharing saves the CPU (one apply), not those uploads. A merged layer reading
   a shared member range directly is out of scope.
7. **Draw lists** (WebGPU render bundles): a shared range that moves (a larger layer, a new group) must advance
   `drawResourcesVersion` as shared textures do (§7), so the other holders' recordings do not replay stale bindings.

## Steps

1. **Measure first**: count, in the default presets and common user setups (polymer cartoon with ligand ball-and-stick,
   plus a molecular surface of the ligand or pocket, spacefill with a surface), the visuals per structure group and
   kind that would share, and the time of `Visual.mark` per visual in hover and select (the many-structures bench,
   `marking.test.ts` fixtures, a browser profile on 1cbs and a large assembly). Continue only if marks of shared
   groups are a visible part of a hover's or a click's CPU.
2. **CPU sharing**: the registry, `createMarkers` taking shared data, marks applied once per share (the mark sequence
   number), `markerActions` in the key, the opt-out of Design 4. Tests: two visuals of a share marked through either
   show the same layer; `Toggle` applied once; a highlight-only representation not sharing with a full one; a visual
   joining a share sees its marks; the last release frees it.
3. **GPU sharing**: the data-page range by `DataArray` on WebGPU and WebGL, offsets per item, byte counts by holders,
   `drawResourcesVersion` on a moved range. Tests: Dawn (two items, one range, one write per mark, read back), WebGL
   data pages.
4. **Measure again**: hover and select per event, uploads per mark, data-page bytes; the many-structures bench and the
   step 1 scenes.

## Open

- Design 4: are marks aimed at one representation (`mark({ repr, ... })`) allowed to show on the other holders of its
  share, or does such a representation opt out?
- Is the structure-level element iterator (`ElementIterator.fromStructure`: labels, plane images, structure-level
  surfaces) worth a key of its own (per structure, not per group)?
- Visuals over equivalent but not identical groups (another structure object of the same model, `includeParent`
  proxies) do not share: a key by group object only. A key by equivalence would need the remapping marks do.
