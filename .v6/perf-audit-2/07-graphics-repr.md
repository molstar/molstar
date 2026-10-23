# Slice 07: graphics repr (v2)

Scope: [`packages/graphics/src/repr`](../../packages/graphics/src/repr): units, complex, volume, particle and shape
representations and their visuals, `Visual.*` (mark, layers, smoothing), loci iteration
(`structure/visual/util/{element,bond,polymer,loci-index}.ts`, `volume/util.ts`), link builders, `merged.ts`. Mesh
builders, marching cubes and smoothing internals are the geo slice; `gl/` and Canvas3D are frozen (WebGPU). MVS is in
[mvs-optimizations](../mvs-optimizations/research/README.md). HEAD 0944823d4 against v1 base 60fe4c2f1; fixtures and
method as in [v1 00-profile](../perf-audit/00-profile.md) (`flat-192`: 11.2k unshared units, 1.1M atoms). Measured
2026-10-22, Node 25, Apple M4.

Changed since v1 (10 commits on the slice): #1995 loci index, per-unit visual lookup, single-pass mark, run-merged
intervals, known-clear layers, no smoothing on merged members; #2005 / #2012 instance-bounded smoothing; #1993 scalar
cartoon curve; #2004 interactions as impostors (`cylinderCount`); #2007 `CylindersBuilder` counts cylinders; #2022
cached foreign-resource check; #2021 `setTransform` bounds from `aTransform`; #1987 isosurface id field. Volume code is
unchanged apart from #1987.

## Hot paths now

| Path | Fixture | Now | Where the time goes |
| --- | --- | ---: | --- |
| hover / select every copy | flat-192, all reprs | 1–3 ms | fixed (P3–P5) |
| one-chain overpaint / transparency | flat-192 | 9–23 / 11–27 ms | fixed (P1, 13/P5) |
| `repr:ball-and-stick` | flat-192 | 4.1 s | `createLinkCylinderImpostors` 39 % (P9, P15), merged copy 18 % (P7, frozen copy) |
| `repr:interactions` | flat-192 | 1.5 s | impostors 54 %; `_compact` 7 % + `create2` 4 % + GC 16 % from the 6× reservation (P11) |
| `repr:molecular-surface` | flat-192 | 20.5 s | core 92 % (`projectTorus` 63 %); repr: P6 |
| `repr:spacefill` / gaussian | flat-192 | 0.68 / ~0.7 s | merged copy 23–24 % (P7) |
| first recolor after build | flat-192 cartoon | 75 ms | `hasForeignResources` 39 % (13/P4, partly fixed), `getMergeKey` 6 % (P17) |
| CPU isosurface rebuild | 400³ | +35 ms, +256 MB | serial id field per build (P12) |
| hover a segment, slice of the segmentation shown | 256³ labels | ~70 ms per scan, ×2 per move | full-grid scan (P13) |

## v1 → now

32 v1 items: **8 fixed**, **2 partly fixed**, **22 still open**, 0 regressed, 0 dropped.

| v1 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| P1 | Layer smoothing redone for every visual | fixed | #1995 (+#2005, #2012) | `isKnownClear` + `reachesAny` early return ([visual.ts:343](../../packages/graphics/src/repr/visual.ts#L343)), no smoothing when untouched (`:378`) |
| P2 | Segment visuals size markers/layers to the whole grid | still-open | | [segment.ts:234](../../packages/graphics/src/repr/volume/segment.ts#L234), `:397` unchanged |
| P3 | Loci application O(visuals × loci.elements) | fixed | #1995 | `loci-index.ts` (WeakMap per loci), `getMarkedVisuals` by unit id |
| P4 | Run-merge never merges; bonds per directed edge | fixed | #1995 | `applyIndexRuns` ([element.ts:198](../../packages/graphics/src/repr/structure/visual/util/element.ts#L198)), `applyIntraBondRuns` ([bond.ts:290](../../packages/graphics/src/repr/structure/visual/util/bond.ts#L290)) |
| P5 | `Visual.mark` applies twice | fixed | #1995 | one pass summing `Interval.size` ([visual.ts:203](../../packages/graphics/src/repr/visual.ts#L203)); see B9 |
| P6 | Molecular-surface radii sized to the model | still-open | | [molecular-surface.ts:48](../../packages/graphics/src/repr/structure/visual/util/molecular-surface.ts#L48) |
| P7 | `renderMerged: 'auto'` by unit count; merged copy | still-open | | policy `base.ts:127`, call [merged.ts:50](../../packages/graphics/src/repr/merged.ts#L50) |
| P8 | Segment mask: 7 `includes` + 7 `dataOffset` per voxel | still-open | | [segment.ts:154](../../packages/graphics/src/repr/volume/segment.ts#L154) |
| P9 | Bond callbacks per half-edge | partly-fixed | #2004 | `radius()` only for multi styles now; position/style/reference/shift still per half-edge |
| P10 | `eachInterBond` linear `includes` per loci atom | still-open | | [bond.ts:430](../../packages/graphics/src/repr/structure/visual/util/bond.ts#L430) |
| B1 | Volume/particle visuals lose state on new data | still-open | | [volume/representation.ts:79](../../packages/graphics/src/repr/volume/representation.ts#L79) |
| B2 | New units visuals never get the repr transform | still-open | | re-measured, also on `mustRecreate` (B2 below) |
| B3 | Labels throw on coarse units | still-open | | [label-text.ts:118](../../packages/graphics/src/repr/structure/visual/label-text.ts#L118) |
| B4 | Segment voxel loci use bbox-local ids | still-open | | [segment.ts:111](../../packages/graphics/src/repr/volume/segment.ts#L111) |
| B5 | Impostor multi/aromatic offsets stale after size change | still-open | | [bond-intra-unit-cylinder.ts:339](../../packages/graphics/src/repr/structure/visual/bond-intra-unit-cylinder.ts#L339) |
| B6 | Disk style stale midpoint | still-open | | [link.ts:563](../../packages/graphics/src/repr/structure/visual/util/link.ts#L563) |
| B7 | Complex/shape position iterator from old geometry | still-open | | [complex-visual.ts:265](../../packages/graphics/src/repr/structure/complex-visual.ts#L265) |
| B8 | Highlight removal mutates `tMarker` only | still-open | | [units-representation.ts:169](../../packages/graphics/src/repr/structure/units-representation.ts#L169) |
| N1 | Segment mask reads out of grid at faces | still-open | | same lines as P8 |
| N2 | Smoothed layers on merged surfaces use member 0's grid | fixed | #1995 | `setLayers(smooth = !merged)` ([units-representation.ts:368](../../packages/graphics/src/repr/structure/units-representation.ts#L368)), re-applied when merging changes (`:291`); tests in `structure/_test/marking.test.ts` |
| N3 | `eachVolumeLoci` segment: loop index as instance; set accumulates | still-open | | [volume/util.ts:118](../../packages/graphics/src/repr/volume/util.ts#L118) |
| 13/P1 | Instanced-surface smoothing spans every instance | fixed | #2005, #2012 | 00-profile: asm-192-low one chain 4.7 s → 20 ms |
| 13/P2 | quality `auto` by root atoms | still-open | | [util.ts:148](../../packages/graphics/src/repr/util.ts#L148) |
| 13/P4 | Recolor pays a per-visual resource check | **partly-fixed** | #2022 | cached only after a first clean call: the first update after a build still walks every cell (39 % of the flat-192 cartoon recolor) |
| 13/P5 | Overpaint/transparency touch every visual | fixed | #1995 | known-clear early return |
| 14/N4 | Cartoon meshes reserve ~1.4× their vertices | still-open | | [polymer-trace-mesh.ts:85](../../packages/graphics/src/repr/structure/visual/polymer-trace-mesh.ts#L85) |
| 12/B1 | Isosurface id field read through the floodfilled `get` | fixed | #1987 | `getIsosurfaceFields` keeps ids on the unfilled space; `volume/_test/isosurface.test.ts` |
| 12/V1 | `createVolumeTexture3d`: 7 `dataOffset` per voxel | still-open | | [volume/util.ts:302](../../packages/graphics/src/repr/volume/util.ts#L302) |
| 12/V2 | `createVolumeTexture2d`: 1 or 7 `dataOffset` per voxel | still-open | | [volume/util.ts:199](../../packages/graphics/src/repr/volume/util.ts#L199) |
| 12/V3 | Segment texture / cells: 7 `dataOffset` + `includes` | still-open | | same code as P8 (tracked there) |
| 12/V4 | Volume dots/points: `get` + `dataOffset` per voxel | still-open | | [dot.ts:178](../../packages/graphics/src/repr/volume/dot.ts#L178) |
| 12/L1 | Slice: `get`/`dataOffset` per pixel | still-open | | [slice.ts:511](../../packages/graphics/src/repr/volume/slice.ts#L511) |

## Open performance findings (confirmed and plausible, ranked)

| # | Id | vs v1 | Title | File | Kind | Impact | Estimate | Effort | Verdict | Frozen |
| ---: | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | P2 | still-open | Segment visuals size markers and layers to the whole grid, per segment | [segment.ts:234](../../packages/graphics/src/repr/volume/segment.ts#L234) | memory | high | −152 MB for 20 segments in 200³; ~2.7 GB for 100 in 300³ (+4× per layer) | S | confirmed | no |
| 2 | P11 | new; **fixed** #2025 | Interactions impostors reserve 6× their cylinder buffers (×6 applied twice) | [link.ts:391](../../packages/graphics/src/repr/structure/visual/util/link.ts#L391) | memory | medium | Float32 alloc 1752 → 312 MB per 1M cylinders; peak RSS −25–30 % in `repr:interactions`; −10–15 % of its time; **measured**: `repr:interactions` alloc flat-32 148–163 → 76–89 MB, asm-192 275–292 → 72–81 MB; asm-192 build 429–505 → 268–290 ms | S | confirmed | no |
| 3 | P7 | still-open | `renderMerged: 'auto'` by unit count; merged copy doubles CPU memory | [merged.ts:50](../../packages/graphics/src/repr/merged.ts#L50) | memory | medium | retained ~2.2× → 1× for large members; −18–24 % of ball-and-stick/spacefill/gaussian builds | S | plausible | policy no / member release yes |
| 4 | P12 | new | Every CPU isosurface build allocates and fills a grid-sized serial id field | [isosurface.ts:91](../../packages/graphics/src/repr/volume/isosurface.ts#L91) | time+memory | medium | −4 B/voxel transient (256 MB at 400³), −35 ms (~5–10 %) per rebuild | S | confirmed | no |
| 5 | 13/P4 | partly-fixed | First update after a build still walks every value cell | [units-visual.ts:220](../../packages/graphics/src/repr/structure/units-visual.ts#L220) | time | low | first recolor at 11k visuals: −29 of 75 ms (cartoon) | S | confirmed | no |
| 6 | P17 | new | Every representation update re-derives every member's merge key | [merged.ts:51](../../packages/graphics/src/repr/merged.ts#L51) | time | low | −3–9 % of each update at 10k+ visuals | S | confirmed | no (repr cache) |
| 7 | P16 | new | `mustRecreate` with a smoothable layer leaves a merged surface unmerged | [units-representation.ts:368](../../packages/graphics/src/repr/structure/units-representation.ts#L368) | time | low | after GPU↔CPU/backend switch: N draw calls instead of 1 and N smoothings (seconds at 10k) | S | confirmed | no |
| 8 | P13 | new | Segment and isosurface loci scan the whole grid per mark on slice and direct-volume visuals | [volume/util.ts:118](../../packages/graphics/src/repr/volume/util.ts#L118) | time | low | ~140 ms per hover move on 256³ → ~1 ms; only for a slice/direct volume of the segmentation itself | M | confirmed | no |
| 9 | P14 | new | Clipping and wiggle layers clear and upload every visual on each change | [visual.ts:648](../../packages/graphics/src/repr/visual.ts#L648) | time | low | 10k visuals ~100+ ms → few ms per change | S | confirmed (code) | no |
| 10 | P9 | partly-fixed | Bond position, style, reference and shift evaluated per half-edge | [link.ts:405](../../packages/graphics/src/repr/structure/visual/util/link.ts#L405) | time | low | ~200–250 ms of 4.1 s ball-and-stick on flat-192 | M | confirmed | no |
| 11 | 14/N4 | still-open | Cartoon meshes reserve ~1.4× their vertices; `compactTrimmed` copies | [polymer-trace-mesh.ts:85](../../packages/graphics/src/repr/structure/visual/polymer-trace-mesh.ts#L85) | time+memory | low | −0.3 s of `repr:auto` at N = 300 | S | confirmed | no |
| 12 | P6 | still-open; **fixed** #2034 | Molecular-surface radii sized to the model's ElementIndex end, per unit | [molecular-surface.ts:48](../../packages/graphics/src/repr/structure/visual/util/molecular-surface.ts#L48) | memory | low | 25 GB cumulative zero-filled alloc on flat-192, ~120 ms; external-memory GC risk | S | confirmed | no |
| 13 | P15 | new | Bond impostor buffers reserved at 2 cylinders per half-edge (fill 0.63), then copied | [link.ts:392](../../packages/graphics/src/repr/structure/visual/util/link.ts#L392) | memory | low | −37 % transient bond-buffer alloc; ~0–4 % time; only with P9 | S | confirmed | no |
| 14 | P18 | new | Sphere impostors evaluate the size theme per atom only for the bounding radius | [element.ts:157](../../packages/graphics/src/repr/structure/visual/util/element.ts#L157) | time | low | ~40–60 ms (~1 %) of ball-and-stick on flat-192 | S | confirmed | no |
| 15 | P8 | still-open | Segment mask: 7 `includes` + 7 `dataOffset` per voxel (= 12/V3) | [segment.ts:154](../../packages/graphics/src/repr/volume/segment.ts#L154) | time | low | `getSegmentCells` 142 → ~20–30 ms per build | S | confirmed | no |
| 16 | 12/V1 | still-open | `createVolumeTexture3d`: 7 `dataOffset` per voxel | [volume/util.ts:302](../../packages/graphics/src/repr/volume/util.ts#L302) | time | low | 1.4–8× on the texture build | S | confirmed | no |
| 17 | 12/V2 | still-open | `createVolumeTexture2d`: 1 or 7 `dataOffset` per voxel | [volume/util.ts:199](../../packages/graphics/src/repr/volume/util.ts#L199) | time | low | `normals` as V1 | S | confirmed | no |
| 18 | 12/V4 | still-open | Volume dots/points: `get` + `dataOffset` per voxel | [dot.ts:178](../../packages/graphics/src/repr/volume/dot.ts#L178) | time | low | small | S | confirmed | no |
| 19 | 12/L1 | still-open | Slice: `get`/`dataOffset` per pixel | [slice.ts:511](../../packages/graphics/src/repr/volume/slice.ts#L511) | time | low | small (color per pixel dominates) | S | confirmed | no |
| 20 | P10 | still-open | `eachInterBond` scans `connectedIndices` per loci atom | [bond.ts:430](../../packages/graphics/src/repr/structure/visual/util/bond.ts#L430) | time | low | O(atoms × connected) → O(min); 2–3 ms today | S | confirmed | no |
| 21 | P19 | new | Scalar cartoon code deoptimized after full GC, boxes doubles until re-optimized (R1) | [curve-segment.ts:164](../../packages/graphics/src/repr/structure/visual/util/polymer/curve-segment.ts#L164) | time+memory | low | +5–8 ms, +40 MB per small cartoon build after a GC | M | plausible | no |
| 22 | 13/P2 | still-open | quality `auto` scores the root's atoms, not geometry per group | [util.ts:148](../../packages/graphics/src/repr/util.ts#L148) | architecture | low | look and cost of large assemblies (question) | S | confirmed | no |

### P2: segment group space = whole grid

- Mechanism: both segment visuals take `groupCount = volume.grid.cells.data.length` (`:234`, `:397`, and the location
  iterators), while mesh ids are `fillSerial` over the expanded bbox (`:188`). Markers and every layer are
  `instanceCount × groupCount` per segment; a whole-segment mark applies a grid-sized interval per visual.
- Evidence: unchanged since v1 (`segment.mjs`: 200³, 20 segments of 6³: 152.6 MB of markers for 4320 vertices; the
  volume is 30.5 MB).
- Fix: `'object'` granularity (the segmentation default, `formats/volume/segmentation.ts:67`): `groupCount = 1`, a
  constant group. `'voxel'`: bbox-local groups, mapped to global cells in `getSegmentLoci`/`eachSegment` via
  `segmentation.bounds[key]` (fixes B4 too).
- Risk: voxel picking and cross-visual highlight need tests; the texture-mesh path packs groups the same way.

### P11: interactions impostors reserve 6×

- Mechanism: #2004 calls `CylindersBuilder.create(linkBuilder.cylinderCount() * 6, linkCount, …)` ("sizes its arrays
  in vertices"), but #2007 had already made `create` count cylinders and multiply by 6 itself
  ([cylinders-builder.ts:66](../../packages/graphics/src/geo/geometry/cylinders/cylinders-builder.ts#L66)). A fresh
  build reserves 36 vertices per cylinder in six Float32 chunked arrays; fill 1/6 < 0.8, so `compactTrimmed` copies.
  The exact count meant to avoid growth forces a 6× reservation plus a copy. All interactions impostor visuals
  (intra, inter, bridge; `withDashedCylinderCount`, `props/computed/representations/shared.ts:40`). Updates that
  reuse buffers (`cylinders` given) are not affected.
- Evidence: `cylsize.ts` rerun, 1M cylinders: 1752 vs 312 MB Float32 allocated, 91–123 vs 61–71 ms. flat-192 profile:
  `_compact` 108 ms, `create2` 67 ms, GC 248 ms of 1.5 s. Auditor A/B (flat-32): peak RSS 729–796 → 555–575 MB.
- Fix: `CylindersBuilder.create(linkBuilder.cylinderCount(), Math.max(1, linkCount / 4), cylinders)`, fix the comment;
  test that a fresh build ends in one exactly filled chunk.
- Risk: none; the count is exact.

### P7: merge policy ignores size (plausible)

- Mechanism: `resolveRenderMerged('auto', n)` is `n > 100`; `createMergedValues` copies every member buffer, members
  keep theirs.
- Evidence: flat-192: merged copy 18 % of ball-and-stick (`copyVertexBlock` 360 ms self), 23–24 % of spacefill and
  gaussian; 13: ball-and-stick retained 2.3 GB at 1M.
- Fix: `'auto'` also requires small members (mean vertex bytes under a budget). Releasing member arrays after merging
  is the frozen alternative (`gl/renderable/merged.ts`).
- Risk: draw calls for hundreds of large members; needs a GPU frame-time measurement first.

### P12: serial id field per isosurface build

- Mechanism: `getIsosurfaceFields` builds `fillSerial(new Int32Array(cells.data.length))` on every CPU mesh and
  wireframe build (`:192`, `:484`); the wireframe is always CPU, the mesh is CPU on floodfill, without graphics or on
  a GPU fallback. In the strided case marching cubes reads `idf[li*idSx+…]`, which equals the cell's data offset. Same
  in v1 (not reported then).
- Evidence: `serial.mjs` rerun: 256³ 9 ms / 67 MB, 400³ 35 ms / 256 MB, 512³ 74 ms / 537 MB. Against marching cubes
  (~400 ms at 400³ by 12's 220³ numbers) that is ~5–10 % of time (the auditor's 10–20 % was high) and doubles the
  transient memory over the Float32 grid.
- Fix: marching cubes takes an implicit serial id (`idField: 'serial'`): unwrapped, the id is the offset it already
  computes (`algorithm.ts:208`); wrapped, the wrapped space's `dataOffset`. A WeakMap cache per grid is the smaller,
  resident alternative.
- Risk: low; test ids for plain, floodfilled and periodic cases (`volume/_test/isosurface.test.ts` has the first two).

### Low

- 13/P4: `createForeignResourcesCheck` remembers only a checked clean pair, and creation never checks, so the first
  `prepareUpdate` after a build walks every cell. flat-192 `theme-apply_cartoon`: `hasForeignResources` 29.3 of
  74.7 ms (39 %), ball-and-stick and surfaces 4–12 ms of `add` self. Fix: mark (renderObject, graphics) clean when the
  visual creates its render object with `graphics` (resources are created with that context, render-object.ts:233).
- P17: `mergeRenderObjectsByType` calls `canMergeRenderObjects` before comparing `prev.key`, and `canMergeValues` runs
  `getMergeKey` (Object.keys + filter + sort + template per define) for every member, on every update. Self time at
  flat-192: 4.3–9.9 ms per recolor (3–9 %). Fix at the repr: when `prev.key === key`, skip the check unless a member's
  define cells, `uInstanceGranularity`, `uLod` or `lodLevels` changed version (fingerprint per render object). Do not
  move the cache into `gl/renderable/merged.ts` (frozen).
- P16 (corrected; was "smoothing applied, then discarded"): on `mustRecreate`, new visuals get `setVisualState` while
  `mergedVisuals` holds the destroyed ones, so `setLayers` smooths them; smoothing switches `overpaintType` to
  `vertexInstance`/`volumeInstance`, `canMergeValues` rejects the group, and the surface stays unmerged with every
  member smoothed. `verify/newvis.ts` (1cbs gaussian, `renderMerged: true`, whole-structure overpaint, forced
  `mustRecreate`): 1 merged object → 3 objects, 0 merged, `vertexInstance` ×3, 3 smoothing calls. Fix: apply layers
  of new visuals after `updateRenderObjectList` decides merging (collect, then `setLayers` once), or predict merging
  from `renderMerged` and the group count. Triggers: GPU↔CPU or backend switch. A new group from a non-equivalent
  structure is not affected (its layers are on the old structure and do not reach it; the plugin re-applies them).
- P13: `eachVolumeLoci` without props (slice `:736`, direct volume `:213`) tests every cell for `Volume.Segment` loci
  (`s.includes`, one singleton per hit) and every cell by `equalEps` for `Volume.Isosurface` loci. `volscan.mjs`:
  256³, 3 labels: 67–73 ms and 519k singletons per scan, before the slice wrapper; highlight + remove doubles it. Only
  for a slice or direct volume of the segmentation volume itself (`Volume.areEquivalent` rejects the tomogram). Fix:
  slice iterates its own sampled cells with a label LUT; direct volume scans `segmentation.bounds` and applies runs;
  fix N3 there.
- P14: `setClipping` / `setWiggle` lack the known-clear and `touched` handling: `createClipping`, a full clear and
  `ValueCell.update` per visual per change. Prism cuts write a clipping layer when kept atoms change (user action, not
  per frame). Fix: the same `isKnownClear` + `reachesAny` early return; share a helper with the four layer setters.
  Not measured (13/P5 measured 106–665 ms for the same pattern at 11k visuals).
- P9: the loop runs over directed edges (`linkCount = 2 × edgeCount`); `position()` (2 `invariantPosition`, plus
  `radiusA/B` and a normalize with `adjustCylinderLength`), `style()`, and for multi styles `referencePosition()` and
  `calculateShiftDir` run for both halves. flat-192: position 158, style 104, reference 38, shift 163 ms total. Fix:
  iterate `a[t] < b[t]`, compute the symmetric parts once, emit both halves by the reverse edge index; radius stays per
  half (`uncertainty` theme). Cylinder order changes; picking ids must stay directed edge indices.
- 14/N4: estimate assumes a tube for every element; sheets use 8 vertices per segment, caps only at SS ends; fill ~0.70.
  Estimate per element from the SS type.
- P6: `new Float32Array(OrderedSet.end(indices))` per unit. `radii.mjs`: 25.2 GB cumulative, 121–134 ms vs 2–3 ms
  compact, ~170 mark-compacts per pass driven by external memory. With 4M live JS objects the same loop takes 19–22 s
  vs 11–91 ms (`verify/radii-heap.mjs`), but the real flat-192 build shows GC at 369 ms (2 %) of 20.5 s (scavenges
  free the young buffers), so that is a risk, not the observed cost. Fix: compact x/y/z/r per unit with
  `indices = Interval(0, n)`.
- P15: without `cylinderCount`, `create(linkCount * 2, …)` reserves 2 per half-edge; most add one (fill 0.63) and
  `compactTrimmed` copies. Auditor A/B with an exact style pre-pass: no time win. Do it with P9 (styles computed once).
- P18: with a non-uniform size theme, `themeSize(l)` per atom only feeds `maxSize`, which is also unused when the old
  bounding sphere is reused; `createSizes` evaluates the theme again. Use the theme's max or the max of the group size
  array.
- P8 (= 12/V3): one membership mask pass over bbox+1 with strides, then the stencil; do it with N1.
- 12/V1, 12/V2, 12/V4, 12/L1: take `Tensor.strides` once and carry offsets (#1987 pattern).
- P10: iterate the smaller side (`connectedIndices` with `OrderedSet.has`, or `e.indices` with `hasEdges`).
- P19 (plausible): after the bench's full GCs, optimized code is marked for deoptimization ("embedded weak objects
  cleared") for curve functions and also for the CIF tokenizer, iterators and theme callbacks, so it is the per-run
  closures/maps dying, not curve-segment itself; the scalar locals box while unoptimized. 7qpd cartoon: 34–41 ms,
  48–82 MB vs v1 32 ms, 38 MB. Check in Chrome before acting; drop if bench-only.
- 13/P2: question for the maintainers: score instanced representations by the largest symmetry group's atoms.

## Bugs

| Id | vs v1 | Severity | Bug | Where | Verdict |
| --- | --- | --- | --- | --- | --- |
| B2 | still-open; **fixed** #2032 | medium | New units visuals never get the representation transform (also on `mustRecreate`) | [units-representation.ts:394](../../packages/graphics/src/repr/structure/units-representation.ts#L394) | confirmed (measured) |
| B1 | still-open; **fixed** #2033 | medium | Volume and particle visuals lose state on new data or recreate | [volume/representation.ts:79](../../packages/graphics/src/repr/volume/representation.ts#L79) | confirmed (code; v1 measured) |
| N1 | still-open; **fixed** #2033 | medium | Segment mask reads wrapped or out-of-grid voxels at grid faces | [segment.ts:154](../../packages/graphics/src/repr/volume/segment.ts#L154) | confirmed (code; v1 measured) |
| B3 | still-open; **fixed** #2032 | medium | Label representation throws on coarse units | [label-text.ts:118](../../packages/graphics/src/repr/structure/visual/label-text.ts#L118) | confirmed (code; v1 measured) |
| B5 | still-open; **fixed** #2034 | low | Impostor multi-bond/aromatic offsets stale after a size-theme change | [bond-intra-unit-cylinder.ts:339](../../packages/graphics/src/repr/structure/visual/bond-intra-unit-cylinder.ts#L339) | confirmed (code; v1 measured) |
| B4 | still-open | low | Segment voxel loci use bbox-local ids as global cells | [segment.ts:111](../../packages/graphics/src/repr/volume/segment.ts#L111) | confirmed |
| B8 | still-open | low | Highlight removal on reused visuals writes `tMarker` without a cell update | [units-representation.ts:169](../../packages/graphics/src/repr/structure/units-representation.ts#L169) | confirmed |
| N3 | still-open | low | `eachVolumeLoci` segment branch: loop index as instance id; set accumulates | [volume/util.ts:129](../../packages/graphics/src/repr/volume/util.ts#L129) | confirmed |
| B9 | new | low | `Visual.mark` takes overlapping intervals for full coverage | [visual.ts:215](../../packages/graphics/src/repr/visual.ts#L215) | plausible |
| B6 | still-open | low | Disk style uses a stale midpoint in impostor and line builders | [link.ts:563](../../packages/graphics/src/repr/structure/visual/util/link.ts#L563) | confirmed (latent) |
| B7 | still-open | low | ComplexVisual and ShapeRepresentation build the position iterator from old geometry | [complex-visual.ts:265](../../packages/graphics/src/repr/structure/complex-visual.ts#L265) | confirmed (latent) |

### B2: transform not applied to new visuals

`setVisualState(visual, group, _state)` for new visuals compares `_state.transform` with itself, so `setTransform` is
skipped. `verify/newvis.ts` (1cbs gaussian): after `setState({transform: +100 Å x})`, a forced `mustRecreate` of every
visual gives tx 100,100,100 → 0,0,0 (the whole representation jumps back); a structure update adding a group gives
100,100,0. Real cases: GPU↔CPU or backend switch after Prism's structure-move gizmo (#2021) or volume streaming set a
repr transform. Fix: an `isNew` flag that bypasses the guard, or `setTransform(_state.transform)` for new visuals.

### B1: volume/particle state lost

`VolumeRepresentation.visual()` creates or recreates visuals for new keys and `mustRecreate`, and `VolumeVisual` sets
`createNew` for a non-equivalent grid (`volume/visual.ts:184`); none re-applies `_state`. Same in
`ParticleRepresentation` (`particles/representation.ts:97`). Scenario: volume streaming
(`plugin/core/src/behavior/dynamic/volume-streaming/transformers.ts:358` re-sets only the transform at apply): each new
box drops alpha, overpaint, transparency, clipping (v1 `volstate.mjs`: visible true, alpha 1 after
`setState({visible: false, alphaFactor: 0.25})`). Fix: call `setVisualState` for new or replaced render objects, with
the transform applied unconditionally.

### N1: segment mask at grid faces

`getSegmentCells` expands the bbox by 2 without clamping; the centre read is unclamped and
`Math.min(xn1 + maxx, x + 1 + minx)` is always ≥ its argument. Reads past a face wrap into the next row or give
`undefined`. v1 `segedge.mjs`, 24³: mesh x = [−2, 25] on a 0..23 grid. Membranes in tomograms touch faces. Fix: a
zero-outside accessor with the grid dimensions, in the P8 mask pass.

### B3: labels on coarse units

`createChainText` calls `auth_asym_id`/`label_asym_id`, and residue/element text read `atomicHierarchy`, on every unit;
`isApplicable` checks only `elementCount > 0`. 8zzh: all three levels throw "Property only available for atomic
models." Fix: coarse accessors at chain level, skip coarse units at residue/element level.

### B9: overlapping intervals taken as full coverage (new, plausible)

The single pass sums `Interval.size` and, when the sum equals `count`, takes the `EveryLoci` outcome and sets every
group. A loci listing a unit twice overcounts. In-tree producer: the superposition panel's atom entries
(`plugin/ui/src/structure/superposition.tsx:532`, one element per selected atom, highlighted with
`highlightOnly(…, false)`, which does not merge elements). Scenario: a cartoon of a 2-residue unit, two selected atoms
in residue 1: two residue-1 intervals, size 2 === count, residue 2 is highlighted too; with a toggle action the
overlap applies twice and cancels. v1's two-pass mark had the same assumption. Fix: take the full path only when the
group loci has unique units (`getPositions` defined, or a per-visual check); otherwise the partial path.

### Low

- B5: radius-scaled offsets are baked into impostor endpoints; geometry is rebuilt on a size change only with
  `adjustCylinderLength` (`:339`), and the structure intra (`:605`) and inter (`bond-inter-unit-cylinder.ts:290`)
  blocks never check size. v1: 12,744 of 38,088 `aStart` floats stale on 7qpd. Rebuild when
  `multipleBonds !== 'off' || aromaticBonds || adjustCylinderLength` and the size theme changed.
- B4: `'voxel'` granularity on a segmentation picks and highlights the wrong voxel; fix with P2.
- B8: after reusing a visual for a changed structure, `applyMarkerAction` on the array only: GPU texture, merged copy and
  `markerStatus` keep the old highlight. Use `visual.mark(EveryLoci, MarkerAction.RemoveHighlight)`.
- N3: `offset = j * cellCount` with `j` the loop counter, not `OrderedSet.getAt(instances, j)`; `set` is shared by all
  loci elements. Wrong instance for Segment loci naming a non-zero instance; fix with P13.
- B6: the Disk branches (`:563`, `:722`) read `vm` left by the previous edge. Mirror the mesh: offset `va` and `vm` by
  `0.475·(vb − va)`. No in-tree impostor/line user.
- B7: `createPositionIterator(geometry, …)` should take `newGeometry || geometry`; shape builds colors (`:200`) before
  refreshing `positionIt` (`:218`).

## Checked and dropped

- P16 as the auditor stated it (new visuals smoothed, then re-applied unsmoothed when merged): refuted by
  `verify/newvis.ts`: smoothing makes the group unmergeable, so the re-application never happens; restated above.
- 13/P4 "fixed": corrected to partly fixed (first update after a build still pays).
- P12 estimate cut from 10–20 % to ~5–10 % of a rebuild; P13 impact cut to low (only a slice/direct volume of the
  segmentation volume itself).

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/1f4e2fbb-da6c-4948-aacb-385a4dab465c/scratchpad/perf-audit-2/graphics-repr/`
(bundle with `node build.mjs <in.ts> <out.mjs>`, esbuild, `molstar-src` condition):

| Script | Covers |
| --- | --- |
| `cylsize.ts` | P11, P15 |
| `bench-fix6.mjs`, `bench-exact.mjs` | P11, P15 A/B (patched bench bundles) |
| `serial.mjs` | P12 |
| `volscan.mjs` | P13 |
| `radii.mjs`, `verify/radii-heap.mjs` | P6 |
| `curve.ts`, `deopt.log`, `opt.log` | P19 |
| `verify/newvis.ts` | B2, P16 (forced `mustRecreate` by patching the bundle's gaussian `mustRecreate`) |
| `../profile/prof-flat-192*` | hot paths, P7, P9, P17, 13/P4 |
