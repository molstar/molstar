# Slice 07: graphics repr

Scope: [`packages/graphics/src/repr`](../../packages/graphics/src/repr): `Representation.createMulti`, units, complex,
volume, particle and shape representations, their visuals (`units-visual.ts`, `complex-visual.ts`, `volume/visual.ts`,
`particles/visual.ts`), `Visual.*` (mark, layers, smoothing), loci iteration (`structure/visual/util/{element,bond,
polymer}.ts`, `volume/util.ts`), link builders and `merged.ts`. Mesh builders, marching cubes and color smoothing
internals belong to the geo slice; `gl/` and Canvas3D are frozen (WebGPU). MVS is covered by
[mvs-optimizations](../mvs-optimizations/research/README.md). Fixtures and method: [00-profile.md](00-profile.md);
`big-bcif` is `1tup-2x2x2` written as a single model (186k atoms, 416 units, no symmetry sharing). Measured
2026-10-19, Node 25.6, Apple M4.

## Hot paths

Build time of most representations is outside this slice: molecular surface spends 96 % in core
`molecular-surface`, gaussian in core density and marching cubes, cartoon in mesh builders plus the merged copy.
Ball-and-stick (big-bcif 456 ms) is ~100 ms `CylindersBuilder`, ~128 ms `createValues`/colors/sizes and ~64 ms of
repr-side link callbacks. The repr-owned costs that matter are in the interactive and update paths, which the
pipeline profile does not cover:

| Path | Fixture | Now | Where the time goes |
| --- | --- | ---: | --- |
| `setState({overpaint})`, layer on 1 of 416 units | big-bcif gaussian, merged | 695–719 ms | `applyMeshOverpaintSmoothing` on all 416 members, 99.7 % (P1) |
| `setState({transparency})`, same | big-bcif gaussian | 458–479 ms | transparency smoothing on all members (P1) |
| `mark` all carbons / every other residue | big-bcif ball-and-stick (577 objects) | 26.1 / 21.4 ms | two `lociApply` passes (P5), per-atom/per-bond singletons (P4), O(visuals × loci.elements) scans (P3) |
| `mark` all but one atom (Interval indices) | big-bcif ball-and-stick | 10.5 ms | ~390k per-bond singleton applies per pass (P4) |
| overpaint carbons, set + clear | big-bcif ball-and-stick | 9.3 ms | one pass of the above |
| segment representation build, 20 small segments | 200³ grid | 34 ms, **153 MB markers** | per-segment group data sized to the whole grid (P2) |

Already efficient: impostor paths, symmetry-group geometry sharing, merged-object reuse key, `EveryLoci` shortcut for
whole-structure marks (`units-representation.ts:276`), word-wise `applyMarkerAction`, instance-granularity fallback,
`Loci.remap` identity early-out.

## Performance findings (confirmed and plausible, ranked)

| Id | Title | File | Kind | Impact | Estimate | Effort | Verdict | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| P1 | Layer smoothing redone for every visual, also those the layers never touch | [visual.ts:303](../../packages/graphics/src/repr/visual.ts#L303) | time | high | 1-chain overpaint on 416-unit surface 700 → few ms; transparency 470 → few ms; restore with layers likewise | S (M with merged) | confirmed | no |
| P2 | Segment visuals size markers/layers to the whole grid, per segment | [segment.ts:234](../../packages/graphics/src/repr/volume/segment.ts#L234) | memory | high | −152 MB on 20 segments in 200³; ~2.7 GB for 100 segments in 300³ (+4× per overpaint layer) | S (object granularity) / M (voxel) | confirmed | no |
| P5 | `Visual.mark` runs every partial `lociApply` twice | [visual.ts:186](../../packages/graphics/src/repr/visual.ts#L186) | time | medium | ~2× on every partial hover/select | S | confirmed | no |
| P4 | Run-merge in `eachElement`/`eachPolymerElement` never merges; bonds applied per directed edge | [element.ts:217](../../packages/graphics/src/repr/structure/visual/util/element.ts#L217) | time | medium | 93k → 16k applies on even-residue loci (1.41 → 1.07 ms per pass); ball-and-stick interval loci ~10 → ~2 ms | S | confirmed | no |
| P3 | Loci application O(visuals × loci.elements); whole-structure test per visual | [units-representation.ts:286](../../packages/graphics/src/repr/structure/units-representation.ts#L286) | time | medium | indexed `eachElement` 1.41 → 0.63 ms per pass; ~10 of 26 ms of a big ball-and-stick mark | M | confirmed | no |
| P7 | `renderMerged: 'auto'` merges by unit count only; merged copy doubles CPU memory | [merged.ts:50](../../packages/graphics/src/repr/merged.ts#L50) | memory | medium | retained ~2.2× buffers (gaussian 73 → 165 MB, molsurf 63 → 141 MB on big-bcif) | S | plausible | no (policy) / yes (member release) |
| P6 | Molecular-surface radii sized to the model's ElementIndex end, per unit | [molecular-surface.ts:48](../../packages/graphics/src/repr/structure/visual/util/molecular-surface.ts#L48) | memory | low | 214 MB cumulative transient alloc on big-bcif (peak < 1 MB per unit); runtime unchanged | S | confirmed | no |
| P8 | Segment mask: 7 `Array.includes` + 7 `dataOffset` closures per voxel | [segment.ts:131](../../packages/graphics/src/repr/volume/segment.ts#L131) | time | low | `getSegmentCells` 142 of 339 ms → ~20–30 ms | S | confirmed | no |
| P10 | `eachInterBond` linear `connectedIndices.includes(v)` per loci atom | [bond.ts:382](../../packages/graphics/src/repr/structure/visual/util/bond.ts#L382) | time | low | O(atoms × connected) → O(min) per unit pair; small today (12 ms / 80 marks) | S | confirmed | no |
| P9 | Bond callbacks evaluated per half-edge | [bond-intra-unit-cylinder.ts:135](../../packages/graphics/src/repr/structure/visual/bond-intra-unit-cylinder.ts#L135) | time | low | ~15–20 ms (3–4 %) of 456 ms ball-and-stick build; radius must stay per half | M | plausible | no |

### P1: layer smoothing on untouched visuals

- Mechanism: `Visual.setOverpaint/setTransparency/setEmissive/setSubstance` clear the array (`clear = true` from
  units/complex visuals), apply layers through `lociApply` (return value ignored), then smooth whenever
  `layers.length > 0`. `UnitsRepresentation.setState` and `setVisualState` (new visuals, restore) do this for all G
  visuals; G−1 smooth an all-zero array (`calcMeshColorSmoothing`, 216 cells per vertex).
- Evidence (`overpaint.mjs`, rerun): 695–719 ms overpaint, 458–479 ms transparency on big-bcif gaussian;
  `applyMeshOverpaintSmoothing` 1959 of 1965 ms in `setOverpaint` (`prof-op`).
- Fix: `touched = lociApply(...) || touched` over layers; when untouched and cleared, skip smoothing. Optionally a hash
  of the pre-smoothing array per visual (next to `overpaintTexture` in `geometry.meta`) to skip unchanged ones when a
  second layer is added.
- Risk: **merged objects**. `getMergeKey` includes every define, but the merged layout key reads only the first
  member (`gl/renderable/merged.ts:418`) and textures take granularity from the first member. Setting `dOverpaint`
  false or leaving `overpaintType` at `groupInstance` on some members while others become `volumeInstance` makes
  members diverge silently. Keep defines uniform within a merged group: for untouched members keep `dOverpaint` true
  with the zeroed `groupInstance` array, and do not smooth members of a merged object at all (see N2). Unmerged case
  (≤ 100 units, `renderMerged: false`) is a clean skip. Test CPU mesh and GPU texture-mesh paths, then clear.

### P2: segment group space = whole grid

- Mechanism: both segment visuals use `groupCount = volume.grid.cells.data.length`
  ([segment.ts:234](../../packages/graphics/src/repr/volume/segment.ts#L234), `:419`), while mesh group ids are
  `fillSerial` over the expanded bbox (`:188`). `createMarkers`, and every `set*` layer, allocate
  `instanceCount × groupCount` per segment. Marking a whole segment applies a cell-count interval per visual.
- Evidence (`segment.mjs`, rerun): 200³ grid, 20 segments of 6³ voxels, 4320 vertices: markers 152.6 MB, all buffers
  154.1 MB (the volume itself is 30.5 MB).
- Fix: segmentations from `formats/volume/segmentation.ts` set picking granularity `'object'`; there use
  `groupCount = 1` and a constant group id (cheap, default case). For `'voxel'`, make groups bbox-local
  (`groupCount` = expanded bbox cells) and map local ↔ global in `getSegmentLoci` and `eachSegment` via
  `segmentation.bounds[key]` (also fixes B4).
- Risk: voxel-granular picking and cross-visual highlight need tests; the texture-mesh path packs groups the same way.

### P5: double `lociApply` in `Visual.mark`

- Mechanism: a size-only pass (`visual.ts:187`) detects empty / whole coverage, then a second pass applies. All
  `eachLocation` work runs twice per visual per hover/select.
- Evidence: `mark.mjs`: ball-and-stick mark carbons 26.1 ms vs overpaint carbons (one pass, no partner checks)
  9.3 ms; both passes show identical subtrees in `prof-mark`.
- Fix: apply while summing `Interval.size`; if 0 return false; if `== count`, take status/average from
  `getMarkerInfo(action, currentStatus)` (array already consistent); else the existing partial logic. Keep
  `previous.loci = EveryLoci` in the full case.
- Risk: low; test toggle/select/highlight sequences and the reverse-action shortcut.

### P4: singleton intervals

- Mechanism: `while (endI < _i && e.indices[endI] === start)` compares against the constant `start`; sorted indices
  are unique, so each index is its own interval (also [polymer.ts:236](../../packages/graphics/src/repr/structure/visual/util/polymer.ts#L236)).
  `eachIntraBond` (`bond.ts:318`) and `eachStructureGroupsBond` (`:460`) apply one `ofSingleton` per directed bond
  inside an `OrderedSet.forEach` closure with an `OrderedSet.has` binary search, even for `Interval` indices.
- Evidence (`eachel.mjs`, rerun): 93,296 → 16,272 apply calls, 1.41 → 1.07 ms per pass; ball-and-stick allButOne mark
  10.5 ms, driven by ~390k bond applies per pass.
- Fix: merge while `indices[endI] === indices[endI-1] + 1`. Bonds: for `Interval` indices and `!isMarking`, one range
  `[offset[min], offset[max+1])`; when marking, extend a pending interval while `t` is contiguous and the partner is
  in the set. Plain loops instead of closures.
- Risk: low; intervals must stay non-overlapping. Unit test against today's marked arrays.

### P3: O(visuals × loci.elements)

- Mechanism: per visual: `invariantId` scan of `loci.elements` (`units-visual.ts:526`, full scan for non-matching
  visuals), `lociIsSuperset` → `isWholeStructure` (sums all element sizes; redundant for `mark`, which already turned
  whole-structure loci into `EveryLoci`), then `eachElement`/`eachIntraBond`/`eachPolymerElement` loop all L elements
  with `unitIndexMap.get`. 416 × 416 × 2 passes ≈ 350k lookups per mark. `Loci.remap` per sub-representation costs
  only on child/component structures (identity early-out otherwise).
- Evidence: `prof-mark` (80 big marks): `eachIntraBond` self 505 ms, `eachElement` self 155 ms, `size` 106 ms,
  invariantId scan 76 ms; `eachel.mjs` indexed variant 0.63 vs 1.41 ms per pass.
- Fix: per call, a `Map<unitId, element>` plus whole-structure flag, cached in a `WeakMap` keyed by the loci (same
  object reaches every sub-representation); visuals iterate `group.units` and look up. Cache `Loci.remap` per
  (loci, structure) there too. Keep the old `eachLocation` signature as fallback (extensions use it).
- Risk: medium-low; test hover, select, overpaint on multi-unit, merged and child structures.

### P7: merge policy ignores size

- Mechanism: `resolveRenderMerged('auto', n)` is `n > 100` ([base.ts:127](../../packages/graphics/src/geo/geometry/base.ts#L127));
  `createMergedValues` copies all member buffers and members keep theirs (sync copies changed blocks).
- Evidence: 00-profile render-buffer table (retained/buffers ≈ 2.2×). Frame-time cost of unmerging not measured.
- Fix: `'auto'` also requires small members (e.g. mean vertex bytes below a budget); leave large meshes unmerged.
  Releasing member CPU arrays after merging is the alternative and is frozen (`gl/renderable/merged.ts`).
- Risk: draw-call count on WebGL/WebGPU for hundreds of surface members; measure before changing the default.

### Low

- P6: radii array is local to one unit's surface task, so the waste is cumulative zero-filled allocation (214 MB over
  a 9.5 s build), not resident memory; scales as units × model atoms (8 GB for 4000 units × 0.5M atoms). The structure
  variant (`:85`) uses serial indices and is not affected. Fix: compact x/y/z/r per unit with `indices =
  Interval(0, n)` (bit-identical, `molsurf.mjs`), or index radius by position in core.
- P8: precompute a `Uint8` membership mask over the bbox + 1 (one `Set`/LUT test per voxel), average from the mask
  with strides; same in `createSegmentTexture2d` (`volume/util.ts:367`) and `eachVolumeLoci` (`util.ts:128`). Do it
  with N1.
- P10: iterate the smaller side (`connectedIndices` with `OrderedSet.has`, or `e.indices` with `b.hasEdges(v)`).
- P9: style, reference position and shift direction are symmetric per undirected bond (shift direction only with a
  reference point; without one it uses `v1`). Radius is **not**: `uncertainty` size theme reads `aIndex` only
  (`theme/size/uncertainty.ts:47`), so radius stays per half-edge. Cache the symmetric parts per undirected bond.

## Bugs

| Id | Bug | Where | Severity | Verdict |
| --- | --- | --- | --- | --- |
| B1 | Volume and particle visuals lose state (transform, alpha, layers, clipping) on new data or recreate | [volume/representation.ts:79](../../packages/graphics/src/repr/volume/representation.ts#L79) | medium | confirmed (measured) |
| B2 | New units visuals never get the representation transform | [units-representation.ts:323](../../packages/graphics/src/repr/structure/units-representation.ts#L323) | medium | confirmed (measured) |
| N1 | Segment mask reads wrapped / out-of-grid voxels at grid faces | [segment.ts:154](../../packages/graphics/src/repr/volume/segment.ts#L154) | medium | confirmed (measured) |
| B3 | Label representation throws on coarse units | [label-text.ts:117](../../packages/graphics/src/repr/structure/visual/label-text.ts#L117) | medium | confirmed (measured) |
| N2 | Smoothed layers on merged surfaces use member 0's grid for all members | [visual.ts:303](../../packages/graphics/src/repr/visual.ts#L303), `gl/renderable/merged.ts:365` | medium | plausible |
| B5 | Impostor multi/aromatic bond offsets stale after size-theme change | [bond-intra-unit-cylinder.ts:339](../../packages/graphics/src/repr/structure/visual/bond-intra-unit-cylinder.ts#L339) | low | confirmed (measured) |
| B4 | Segment voxel loci use bbox-local ids as global cell indices | [segment.ts:111](../../packages/graphics/src/repr/volume/segment.ts#L111) | low | confirmed |
| B8 | Highlight removal mutates `tMarker` without cell update or status refresh | [units-representation.ts:140](../../packages/graphics/src/repr/structure/units-representation.ts#L140) | low | confirmed |
| B6 | Disk link style uses a stale midpoint in impostor and line builders | [link.ts:555](../../packages/graphics/src/repr/structure/visual/util/link.ts#L555) | low | confirmed (latent) |
| B7 | Complex visual / shape build the position iterator from old geometry | [complex-visual.ts:262](../../packages/graphics/src/repr/structure/complex-visual.ts#L262) | low | confirmed (latent) |
| N3 | `eachVolumeLoci` segment branch uses loop index as instance id; segment set accumulates across elements | [volume/util.ts:129](../../packages/graphics/src/repr/volume/util.ts#L129) | low | confirmed |

### B1: volume/particle state lost

`VolumeVisual.prepareUpdate` sets `createNew` for any non-identical grid or new key (`volume/visual.ts:184`);
`VolumeRepresentation.visual()` builds new visuals for new keys and `mustRecreate`; neither re-applies `_state`
(units/complex do via `setVisualState`). `ParticleRepresentation` (`particles/representation.ts:97`) has the same gap.
`volstate.mjs`: after `setState({visible: false, alphaFactor: 0.25})` and a new volume, the new object has visible
true, alpha 1, while `repr.state.visible` stays false. Real case: volume streaming update
(`plugin/core/src/behavior/dynamic/volume-streaming/transformers.ts:366`) re-sets only `visible`; the model transform
set at apply is dropped on the first streamed box, and layer transforms see an unchanged `repr.state` and do not
re-apply. Fix: track visuals created or whose `renderObject` changed in `createOrUpdate` and call
`setVisualState(visual, _state)` for them (with the transform applied unconditionally, see B2).

### B2: transform not applied to new visuals

`setVisualState(visual, group, _state)` for new visuals hits
`transform !== _state.transform || !Mat4.areEqual(...)`, false when `state === _state`. `transform.mjs`: 2 → 8 units
after `setState({transform: +100 Å x})`: objects 0–1 have tx 100, 2–7 have tx 0. Any component update adding groups,
or a `mustRecreate` (GPU/CPU switch), splits the representation in space. Fix: an `isNew` flag that bypasses the
guard (or `setTransform(_state.transform)` when not identity).

### N1: segment mask wraps at grid edges (new)

`getSegmentCells` expands the bbox by 2 without clamping to the grid; the centre read `o(x + minx, …)` is unclamped
and the upper "clamp" `Math.min(xn1 + maxx, x + 1 + minx)` is a no-op (always ≥ the argument). Reads at x = −2..−1 or
past the last column wrap into the neighbouring row (negative offsets give `undefined`). `createSegmentTexture2d`
(`volume/util.ts:396`) is identical. `verify/segedge.mjs`, 24³ grid, segment spanning the full x range,
y, z ∈ [8, 12): mesh x = [−2.00, 25.00] (grid is 0..23) and y = [6.50, 12.40] instead of [7.5, 11.5]; spurious
geometry outside the box. Common for membranes in tomograms. Fix: read through a clamped (or zero-outside) accessor
with the grid dimensions, in the P8 mask pass.

### B3: labels on coarse units

`createChainText` (`auth_asym_id`/`label_asym_id`), `createResidueText` and `createElementText` call atomic-only
accessors on every unit; `isApplicable` only needs `elementCount > 0`. `label.mjs` on 8zzh: all three levels throw
"Property only available for atomic models."; one coarse unit aborts a mixed IHM structure. Fix: coarse accessors at
chain level, skip (or label by coarse seq range) at residue/element level.

### N2: smoothed layers on merged surfaces (new, plausible)

Smoothing sets each member's `overpaintType` to `volumeInstance` and its own `tOverpaintGrid`/grid transform. The
merged object spreads non-own cells from the first member (`gl/renderable/merged.ts:365`), takes layout key and
texture granularity from the first member only, and is not re-checked with `canMergeValues` (which rejects
`volumeInstance`) until the next `createOrUpdate`. So every member samples member 0's grid: an overpaint on any chain
but the first of a big-bcif gaussian surface (merged, `smoothColors: auto`) should not show. Not rendered here (needs
GL). Repr-side fix without touching `gl`: skip layer smoothing for members of a merged object, or re-run
`mergeRenderObjectsByType` after layer changes. Root is in the frozen merged renderable.

### Low

- B5: `createLinkCylinderImpostors` bakes radius-scaled offsets (`absOffset`, `multipleOffset`, `aromaticOffset`,
  `link.ts:429/437/476`) into endpoints; update state recreates geometry on a size change only with
  `adjustCylinderLength`. `bondsize.mjs` on 7qpd: 12,744 of 38,088 `aStart` floats differ from a fresh build. Fix:
  also rebuild on a size-theme change when `multipleBonds !== 'off' || aromaticBonds || adjustCylinderLength` in all
  three impostor visuals: units intra (`:339`), structure intra (`:605` block, no size check at all) and inter
  (`bond-inter-unit-cylinder.ts:290` block, no size check at all). Mesh visuals rebuild anyway (no `uSize`). Costs a
  geometry rebuild per size change.
- B4: picking a voxel of a segment (granularity `'voxel'`, not the segmentation default `'object'`) reports a wrong
  cell and highlights a different voxel in isosurface/slice visuals. Fix with P2.
- B8: after reusing a visual for a changed structure, highlight removal writes the array only; when markers were not
  recreated, the GPU texture, a merged copy and a uniform `markerStatus` keep the old highlight. Use
  `visual.mark(EveryLoci, MarkerAction.RemoveHighlight)`.
- B6: the `Disk` branch reads `vm` set by the previous edge. Do not just compute the midpoint (that gives two
  off-centre disks per bond in the half-edge impostor); mirror the mesh: offset `va` and `vm` by `0.475·(vb − va)`.
  In-tree Disk users (rcsb clashes) use the mesh builder.
- B7: `createPositionIterator(geometry, …)` should take `newGeometry || geometry` (units-visual does); shape builds
  colors/sizes (`shape/representation.ts:200`) before refreshing `positionIt` (`:218`). Mesh position iterators
  capture `vertexBuffer` and `vertexCount` at creation, so vertex-granular themes read stale data when the builder
  returns a new object or grows the buffer.
- N3: `offset = j * cellCount` with `j` the loop index, not `OrderedSet.getAt(instances, j)`, and `set` is not reset
  per loci element. Wrong instance when a Segment loci names a non-zero instance on isosurface/slice visuals.

## Checked and dropped

None refuted. Downgraded: P6 high-memory claim (transient per unit, peak < 1 MB; low), P9 estimate (radius is
direction-dependent under the `uncertainty` theme; ~half the claimed saving). Corrected fix: B6 (midpoint fix alone
is wrong for half-edge impostors). P1 fix needs the merged-define constraint above.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/ac6cef1a-5954-43a5-ae83-c22fef35f065/scratchpad/perf-audit/graphics-repr/`
(bundle with `../build.mjs <in.ts> <out.mjs>`, esbuild, `molstar-src` condition):

| Script | Covers |
| --- | --- |
| `overpaint.ts` | P1 |
| `mark.ts`, `mark2.ts`, `markprof.ts`, `prof-mark/` | P3, P4, P5 |
| `eachel.ts` | P3, P4 |
| `segment.ts`, `segment2.ts`, `prof-seg/` | P2, P8 |
| `molsurf.ts` | P6 |
| `volstate.ts` | B1 |
| `transform.ts` | B2 |
| `label.ts` | B3 |
| `bondsize.ts` | B5 |
| `verify/segedge.ts` | N1 |
| `subtree.mjs`, `subtree2.mjs`, `prof-op/` | profile subtrees (P1, P9) |
