# Slice 07: graphics repr (v3)

Scope: [`packages/graphics/src/repr`](../../packages/graphics/src/repr), as in [v2](../perf-audit-2/07-graphics-repr.md):
units, complex, volume, particle and shape representations and visuals, `Visual.*` (mark, layers), loci iteration,
link builders, `merged.ts`. Mesh builders, marching cubes, `gl/` and Canvas3D are other slices (`gl/` frozen).

HEAD `16277d822` against v2 base `0944823d4` (A/B: the HEAD `tests/src/pipeline/bench.ts` bundled twice, once with the
base's `packages/{core,io,model,graphics}` sources). Measured 2026-10-25, Node 26, Apple M4, **load 25–35** (13 other
auditors): wall times swing 2–10× between runs and are given only as ranges; allocation (`total_allocated_bytes`),
retained MB, buffer bytes and counters are the evidence.

PRs since v2 that touched the slice (9): #2025 (07/P11), #2030 (labels: `withFontAtlasFonts`), #2032 (B2, B3), #2033
(B1, N1), #2034 (P6, B5), #2037 (periodic: unit cell origin), #2044 (N3), #2045 (P14), #2048 (13/P4, P8, P10, P16, P17,
P18, B6, B7, N1 again). #2042 did not touch the slice.

## Hot paths now

Repr-only stages (`--only 'interactions,repr:'`: required stages unmeasured, the interactions property in its own
stage). `buf` = the representation's render objects; `members` = arrays the merged objects' members keep (07/P7).

| Path | Fixture | Alloc HEAD / v2 base (MB) | Retained (MB) | Buffers + member copies (MB) | Where it goes |
| --- | --- | ---: | ---: | ---: | --- |
| `repr:cartoon` | 7qpd | 56–77 / 51–78 | 14–15 | 13.8 | unchanged (P19 plausible) |
| `repr:cartoon` | flat-32 | 135–220 / 215–293 | 101–104 | 35.1 (32.2 of it in members) | |
| `repr:ball-and-stick` | flat-32 | 196–244 / 400–455 | 382–387 | 165.2 + 168.1 | merged copy doubles (P7); half the alloc of v2 (not attributed) |
| `repr:spacefill` | flat-32 | 71–86 / 124–136 | 37 | 4.5 + 7.5 | |
| `repr:gaussian-surface` | flat-32 | 133–163 / 158–174 | 170 | 72.8 + 73.2 | merged copy doubles (P7) |
| `repr:molecular-surface` | flat-32 | 207–335 / 1513–1547 | 151 | 62.8 | **P6 fixed** (#2034): −80–85 % alloc |
| `repr:interactions` | flat-32 | 91–100 / 148–150 | 141 | 82.8 + 39.5 | **P11 fixed** (#2025): −35 % |
| `repr:interactions` | asm-192 | 128–133 / 272–278 | 266 | 265.0 | one non-instanced inter-unit visual over 254,528 contacts (P20), 6 vertices per cylinder (08/P7, frozen) |
| `repr:*` other | asm-192 | 4–35 / 6–49 | ≤ 11 | ≤ 8.9 | instanced, nothing merged |
| recolor (`chain-id`) | flat-192 spacefill, 1.5 s at load ~30 | 54 | – | – | `createGroupColor` 34 %, GC 25 %, **`mergeRenderObjectsByType` 113 ms (7.5 %)**: the merge check per member on every update (P17 regressed) |
| segment hover | 200³, 20 segments | – | – | 152.6 MB markers | grid-sized markers per segment (P2): 44 ms per hover+unhover pair |

Times (best of runs, load 25+): 7qpd cartoon 29–34 ms, b&s 32–39, spacefill 8, gaussian 356–379, molsurf 934–1079,
interactions 31–68; flat-32 cartoon 156–331, b&s 448–1392, gaussian 787–1659, molsurf 10.8–19.7 s (base 39–57 s),
interactions (repr only) 233–384 ms (base 550–1553).

## v2 → now

33 v2 items (22 perf, 11 bugs): **19 fixed**, 0 partly, **13 still open**, **1 regressed** (P17: fixed in #2048 but
ineffective), 0 dropped. New: 1 perf (P20), 0 bugs.

| v2 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| P2 | Segment group space = whole grid, per segment | still open | | re-measured: 200³, 20 segments of 6³: 152.6 MB `tMarker` for 4320 vertices (grid 30.5 MB); +44 ms per hover pair |
| P11 | Interactions impostors reserve 6× twice | fixed | #2025 | `create(cylinderCount(), max(1024, linkCount/4))` ([link.ts:391](../../packages/graphics/src/repr/structure/visual/util/link.ts#L391)); flat-32 repr alloc 148–150 → 91–100 MB, asm-192 272–278 → 128–133 MB (v2's patched estimate was 72–81) |
| P7 | `renderMerged: 'auto'` by count; merged copy | still open; **fixed in #2051** | #2051 | the members' arrays are the merged object's geometry and layer parts (`SegmentedArray`), no copy; `auto` from 30 units measured on both backends ([merged-render-objects.md](../plans/gl/merged-render-objects.md) steps 1, 3, 4); before: re-measured: member arrays 168.1 MB beside 165.2 MB merged (flat-32 b&s), 73.2 beside 72.8 (gaussian) |
| P12 | Serial id field per CPU isosurface | still open | | [isosurface.ts:91](../../packages/graphics/src/repr/volume/isosurface.ts#L91) unchanged |
| 13/P4 | First update after a build walks every cell | fixed | #2048 | `Visual.foreignResourcesCheck().created` ([visual.ts:129](../../packages/graphics/src/repr/visual.ts#L129)); reviewed, sound |
| P17 | Merge key re-derived per member per update | **regressed** (fix ineffective) | #2048 | the fingerprint watches `uLod` and `uClipObjectType`, which `BaseGeometry.updateValues` bumps on every update: see P17 below |
| P16 | `mustRecreate` with a layer leaves a surface unmerged | fixed | #2048 | layers of new visuals wait for `updateRenderObjectList` ([units-representation.ts:281](../../packages/graphics/src/repr/structure/units-representation.ts#L281)) |
| P13 | Segment/isosurface loci scan the grid on slice and direct volume | still open | | [volume/util.ts:127](../../packages/graphics/src/repr/volume/util.ts#L127); N3 part fixed |
| P14 | Clipping and wiggle clear every visual | fixed | #2045 | known-clear + `reachesAny` ([visual.ts:688](../../packages/graphics/src/repr/visual.ts#L688), `:744`) |
| P9 | Bond callbacks per half-edge | still open | | [link.ts:406](../../packages/graphics/src/repr/structure/visual/util/link.ts#L406) unchanged |
| 14/N4 | Cartoon reserves ~1.4× vertices | still open | | [polymer-trace-mesh.ts:86](../../packages/graphics/src/repr/structure/visual/polymer-trace-mesh.ts#L86) unchanged |
| P6 | Molsurf radii sized to the model per unit | fixed | #2034 | compact x/y/z/r per unit ([molecular-surface.ts:51](../../packages/graphics/src/repr/structure/visual/util/molecular-surface.ts#L51)); flat-32 alloc 1513–1547 → 207–335 MB |
| P15 | Bond impostors reserved at 2 per half-edge | still open | | [link.ts:393](../../packages/graphics/src/repr/structure/visual/util/link.ts#L393) |
| P18 | Size theme per atom only for the bounding radius | fixed | #2048 | evaluated only for a new sphere |
| P8 | Segment mask 7 `includes` per voxel | fixed | #2048 | `getSegmentMaskValue` with a box ([volume/util.ts:376](../../packages/graphics/src/repr/volume/util.ts#L376)) |
| 12/V1 | `createVolumeTexture3d` 7 `dataOffset` per voxel | still open; **fixed in #2051** | #2051 | read at `Tensor.strides` (`getStridedGrid`), output identical ([direct-volume-performance.md](../plans/gl/direct-volume-performance.md) step 7); before: [volume/util.ts:303](../../packages/graphics/src/repr/volume/util.ts#L303) unchanged |
| 12/V2 | `createVolumeTexture2d` | still open | | [volume/util.ts:200](../../packages/graphics/src/repr/volume/util.ts#L200) |
| 12/V4 | Volume dots `get` + `dataOffset` per voxel | still open | | [dot.ts:179](../../packages/graphics/src/repr/volume/dot.ts#L179) |
| 12/L1 | Slice `get`/`dataOffset` per pixel | still open | | [slice.ts:511](../../packages/graphics/src/repr/volume/slice.ts#L511) |
| P10 | `eachInterBond` scans connected atoms per loci atom | fixed | #2048 | iterates the smaller side ([bond.ts:438](../../packages/graphics/src/repr/structure/visual/util/bond.ts#L438)) |
| P19 | Scalar cartoon code deoptimized after GC (plausible) | still open | | 7qpd cartoon 29–34 ms, 56–77 MB: as v2 (34–41 ms, 48–82 MB) |
| 13/P2 | quality `auto` by root atoms | still open | | question for the maintainers, unchanged |
| B2 | New units visuals never get the repr transform | fixed | #2032 | `isNew` → `!Mat4.isIdentity(transform)` |
| B1 | Volume/particle visuals lose state | fixed | #2033 | `setVisualState` when the render object changes |
| N1 | Segment mask reads past grid faces | fixed | #2033, #2048 | zero outside the grid |
| B3 | Labels throw on coarse units | fixed | #2032 | coarse `asym_id`; residue/element levels skip coarse units |
| B5 | Impostor multi-bond offsets stale after size change | fixed | #2034 | `bondImpostorGeometryReadsSize` in all three impostor blocks |
| B4 | Segment voxel loci use bbox-local ids | still open | | [segment.ts:111](../../packages/graphics/src/repr/volume/segment.ts#L111) with ids `fillSerial` over the box (`:175`) |
| B8 | Highlight removal writes `tMarker` without a cell update | still open | | [units-representation.ts:169](../../packages/graphics/src/repr/structure/units-representation.ts#L169) |
| N3 | Segment branch: loop index as instance; set shared | fixed | #2044 | |
| B9 | `Visual.mark` takes overlapping intervals as full coverage | still open | | [visual.ts:235](../../packages/graphics/src/repr/visual.ts#L235) unchanged |
| B6 | Disk style stale midpoint | fixed | #2048 | |
| B7 | Position iterator from old geometry | fixed | #2048 | |

New bugs or costs from fixes since v2: **1** (P17's fingerprint never matches after an update, #2048: the fix is
inert, no harm beyond the fingerprint bookkeeping). Reviewed and sound: #2025 link count, #2032 `isNew` transform,
#2033 volume/particle state, #2034 compact molsurf positions (`id` stays ordinal-indexed, as `calcMolecularSurface`
reads it), #2044 per-element segment set, #2045 clipping/wiggle known-clear (same outcome as a full pass), #2048
foreign-resources `created` (conservative: a reused old geometry of another context is not marked clean), deferred
layers of new visuals (merged objects pick up members' later layer changes, as for any `setState`), P10 smaller-side
iteration, segment mask box (falls back to the stencil outside the box), `getCartoonParams` cache per structure (no
caller mutates the result; the `ctx` argument is unused by the params). #2030's `withFontAtlasFonts` is synchronous
when no web font is pending.

## Open performance findings (ranked)

| # | Id | Status | Title | File | Kind | Estimate | Effort | Frozen |
| ---: | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | P2 | still open | Segment visuals size markers and layers to the whole grid, per segment | [segment.ts:221](../../packages/graphics/src/repr/volume/segment.ts#L221) | memory+time | measured: −152.6 MB for 20 segments in 200³ (5× the grid), −44 ms per hover pair; ~2.7 GB for 100 in 300³ (est.) | S | no |
| 2 | P7 | **fixed in #2051**: members' arrays as parts, no copy; `auto` from 30 measured ([merged-render-objects.md](../plans/gl/merged-render-objects.md) steps 1, 3, 4) | `renderMerged: 'auto'` by count; members keep their arrays beside the merged copy | [merged.ts:107](../../packages/graphics/src/repr/merged.ts#L107) | memory | measured: +168 MB beside 165 MB (flat-32 b&s), +73 beside 73 (gaussian), +40 beside 83 (interactions); ~2× retained | S (policy) / M (release, frozen) | policy no / release yes |
| 3 | P20 | new | Inter-unit visuals (interactions, inter bonds) are one non-instanced complex visual on assemblies: every copy's contacts drawn as geometry | [interactions-inter-unit-cylinder.ts:326](../../packages/graphics/src/props/computed/representations/interactions-inter-unit-cylinder.ts#L326) | memory+time | measured: asm-192 interactions 265 MB buffers (other four reprs ≤ 9 MB each); inter contacts per copy constant (1,326 asm-192, 1,308 flat-32); est. 265 → ~2–3 MB + instances if they repeat per operator | L | no |
| 4 | P17 | **regressed**; **fixed in #2050** (`lodOf` by value, only what the merge check reads) | Merge fingerprint (#2048) misses on every update: `uLod`, `uClipObjectType` always re-versioned | [merged.ts:30](../../packages/graphics/src/repr/merged.ts#L30) | time | measured: 416/416 members re-checked per recolor or prop change (flat-32); flat-192 spacefill recolor: 113 ms of 1.5 s (7.5 %, load ~30) in `mergeRenderObjectsByType` | S | no |
| 5 | P12 | still open | Serial id field per CPU isosurface build | [isosurface.ts:91](../../packages/graphics/src/repr/volume/isosurface.ts#L91) | time+memory | v2 measured: 256 MB / 35 ms at 400³ | S | no |
| 6 | P9 | still open | Bond position, style, reference and shift per half-edge | [link.ts:406](../../packages/graphics/src/repr/structure/visual/util/link.ts#L406) | time | v2: ~200–250 ms of 4.1 s b&s on flat-192 | M | no |
| 7 | P13 | still open | Segment/isosurface loci scan the grid per mark on slice and direct volume | [volume/util.ts:127](../../packages/graphics/src/repr/volume/util.ts#L127) | time | v2: ~140 ms per hover move at 256³ | M | no |
| 8 | 14/N4 | still open | Cartoon reserves ~1.4× its vertices, then copies | [polymer-trace-mesh.ts:86](../../packages/graphics/src/repr/structure/visual/polymer-trace-mesh.ts#L86) | time+memory | v2: −0.3 s of `repr:auto` at N = 300 | S | no |
| 9 | P15 | still open | Bond impostors reserved at 2 cylinders per half-edge | [link.ts:393](../../packages/graphics/src/repr/structure/visual/util/link.ts#L393) | memory | v2: −37 % transient bond buffers | S | no |
| 10 | 12/V1, V2, V4, L1 | V1 **fixed in #2051** ([direct-volume-performance.md](../plans/gl/direct-volume-performance.md) step 7); V2, V4, L1 still open | Volume texture, dots, slice: `get`/`dataOffset` per voxel | [volume/util.ts:303](../../packages/graphics/src/repr/volume/util.ts#L303) | time | v2: 1.4–8× on the texture build | S | no |
| 11 | P19 | still open (plausible) | Cartoon curve code deoptimized after a full GC | [curve-segment.ts:164](../../packages/graphics/src/repr/structure/visual/util/polymer/curve-segment.ts#L164) | time | +5–8 ms per small build (v2) | M | no |
| 12 | 13/P2 | still open | quality `auto` by root atoms | [util.ts:169](../../packages/graphics/src/repr/util.ts#L169) | architecture | question | S | no |

### P2: segment group space = whole grid (still open, re-measured)

- `seg.ts` (200³ `Uint32` labels, 20 cube segments of edge 6, CPU mesh, `tryUseGpu: false`): 20 render objects, 4320
  vertices, **152.6 MB of `tMarker`** (20 × 8M groups) of 152.8 MB buffers; build 164 ms. Hover + unhover of one
  segment: **43.7 ms per pair**, from `eachVolumeLoci`'s grid-sized interval (`:110`) applied to the marked visual's
  8M-group marker array and its average.
- Fix as v2: `'object'` granularity (the segmentation default) → `groupCount = 1`; `'voxel'` → box-local groups mapped
  through `segmentation.bounds[key]` (fixes B4). Both CPU (`:221`) and texture mesh (`:384`, `:409`).

### P7: members keep their arrays (still open, re-measured)

`bench07` (HEAD bench + member byte count): flat-32 ball-and-stick 3 merged objects, 576 members holding **168.1 MB**
beside the merged objects' 165.2 MB, retained 382 MB; gaussian 73.2 beside 72.8 (retained 170); interactions 39.5
beside 82.8; spacefill 7.5 beside 4.5; cartoon partly merged, 32.2 of 35.1 MB in members. The merged copy is a CPU
duplicate of every member buffer; the release of member arrays after upload is in frozen `gl/renderable/merged.ts`.
Policy fix (not frozen, `geo/geometry/base.ts:127`): `'auto'` merges only small members (mean vertex bytes under a
budget). Needs a browser GPU draw-call check first (as v2).

### P20: inter-unit visuals are not instanced on assemblies (new)

- `repr:interactions` on asm-192 (192 operators × the 9 units of 1tup; every other representation instanced): 6
  render objects, **265.0 MB of buffers, 266 MB retained** (`aMapping`/`aStart`/`aEnd` 52.7 MB each, `elements` 35.1),
  against ≤ 9 MB for cartoon, ball-and-stick, spacefill and gaussian together; repr alloc 128–133 MB, 2.1–5.2 s at load.
  The inter-unit contacts: **254,528 = 1,326 per copy** on asm-192 and 41,864 = 1,308 per copy on flat-32: the count per copy
  does not grow with the copy count, consistent with contacts between units of one copy that repeat under each
  operator (not checked edge by edge).
- Cause: `InteractionsInterUnitImpostorVisual` (and `InterUnitBondCylinder*`, the `structure-*` visuals) is a
  `ComplexVisual`: one geometry over all inter-unit edges of the structure, instance count 1, no use of the symmetry
  groups the units visuals instance by.
- Fix (L): split inter-unit edges into those between units of one operator (geometry built once from the first
  copy's unit pair class, instanced by the operator's transforms, as `UnitsVisual` does) and the rest (cross-operator
  contacts, kept complex). Estimate: 265 MB → ~2–3 MB plus per-instance data on asm-192; asm-864 1.2 GB (v2) likewise.
  Picking ids need (operator, edge) mapping; shared with slice 09 (`props/computed/representations`). Independent of
  08/P7 (6 vertices per cylinder, frozen), which multiplies whatever is drawn.

### P17: the merge fingerprint never matches after an update (regressed, new cost of #2048's fix)

- `isMergeRelevantKey` (`merged.ts:27`) records `d*` cells, every key ending in `Type`, `uInstanceGranularity`, `uLod`,
  `lodLevels`. `BaseGeometry.updateValues` runs on every visual update that does not create a render object
  (`units-visual.ts:397`, `complex-visual.ts:295`) and does `ValueCell.update(values.uLod, …)`
  ([base.ts:265](../../packages/graphics/src/geo/geometry/base.ts#L265)) and `ValueCell.update(values.uClipObjectType,
  …)` (`:250`, matched by `endsWith('Type')` though `canMergeValues` never reads it) unconditionally. Every member's
  fingerprint is stale after any update, so `canMergeRenderObjects` runs as before #2048, plus the fingerprint
  rebuild (`Object.keys` and two arrays per member).
- Probe `fp.ts` (spacefill, `renderMerged: true`, `Array.prototype.sort` counted during the update; `getMergeKey`
  sorts once per member): 7qpd 7/7 members, flat-32 **416/416 sorts** for each of a recolor (`chain-id`), a second
  recolor and `alpha: 0.5`; uLod/uClipObjectType versions changed on 416/416 members. `canMergeRenderObjects` over the
  416 members: 13.6 ms (32.7 µs/member, load 27).
- flat-192 spacefill `theme-apply:spacefill:chain-id` CPU profile (`out/prof-sf192`, 1.51 s wall at load ~30):
  `mergeRenderObjectsByType` 113 ms total (7.5 %): `canMergeRenderObjects` 65 ms (`getMergeKey` 37), the fingerprint
  bookkeeping (`recordMergeFingerprint` 13 ms self, the `every` check 22 ms self) on top. v2 counted 4.3–9.9 ms of
  `getMergeKey` self per flat-192 recolor; the bookkeeping makes the path about 1.5× what it was without the fix.
- The unit test (`_test/merge-by-type.test.ts`) calls `mergeRenderObjectsByType` twice without updating the visuals in
  between, so it passes.
- Fix (S): fingerprint only what `canMergeValues` reads: `d*`, `colorType`, `sizeType`, the six layer `*Type` cells,
  `uInstanceGranularity`, `lodLevels`, and compare `uLod` by value (or `updateValues` uses `updateIfChanged` with
  `Vec4.exactEquals` for `uLod`); a test that recolors a merged representation and counts merge checks.

### Low

- P12, P9, P13, 14/N4, P15, 12/V*, P19, 13/P2: code unchanged since v2; v2's evidence and fixes stand.

## Bugs

| Id | Severity | Status | Title | File | Failure scenario |
| --- | --- | --- | --- | --- | --- |
| B4 | low | still open | Segment voxel loci use box-local ids as global cells | [segment.ts:111](../../packages/graphics/src/repr/volume/segment.ts#L111) | 'voxel' picking granularity on a segmentation: the wrong voxel picked and highlighted |
| B8 | low | still open | Highlight removal on reused visuals writes `tMarker` only | [units-representation.ts:169](../../packages/graphics/src/repr/structure/units-representation.ts#L169) | structure change while hovering: stale highlight on GPU and in the merged copy until the next mark |
| B9 | low | still open (plausible) | `Visual.mark` takes overlapping intervals as full coverage | [visual.ts:235](../../packages/graphics/src/repr/visual.ts#L235) | superposition atom entries (two atoms of one residue) on a polymer-granular visual: whole visual highlighted |

No new bugs found in the fixes since v2 (the P17 regression is a cost, not a correctness issue: merging is still
decided by the full check). Cosmetic: `getSegmentMaskValue` has two stacked doc comments, the first orphaned
([volume/util.ts:368](../../packages/graphics/src/repr/volume/util.ts#L368)).

## Already good, re-checked

- Interactions repr after #2025: one exactly filled chunk; flat-32 alloc −35 % against the v2 base in the same setup.
- Molecular surface after #2034: flat-32 alloc 1.5 GB → 0.2–0.3 GB, molsurf build 39–57 s → 11–20 s at the same load.
- Ball-and-stick alloc halved against the v2 base on flat-32 (400–455 → 196–244 MB) with identical buffers; not
  attributed to one PR (outside this slice's diff or spread over #2034/#2048).
- Render-buffer bytes identical HEAD vs base on every repr and fixture (7qpd, flat-32, asm-192): no fix changed what is
  drawn.

## Not audited / not measured

- Wall-time A/B: the machine ran at load 25–35; times are not comparable beyond ranges.
- GPU paths (texture mesh, compute isosurfaces): no GPU in Node.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/31ee827f-4715-4967-bec3-085f5013d11f/scratchpad/pa3/07/`
(`node build.mjs <in.ts> <out.mjs> [baseRoot]`; `base/` = `git archive 0944823d4` of core/io/model/graphics):

| Script | Covers |
| --- | --- |
| `runall.sh` → `out/log.txt` | repr stage A/B, 7qpd / flat-32 / asm-192 (`--only 'repr:'`, interactions property inside its stage) |
| `run2.sh`, `src/bench.ts` (member bytes) → `out/log2.txt` | repr-only A/B with member buffers (P7, P11) |
| `seg.ts` | P2 |
| `fp.ts` | P17 |
