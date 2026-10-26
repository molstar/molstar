# Slice 05: model-structure (v3)

Scope: `packages/model/src/model/structure` (`Structure`, `Unit`, element/Loci/Bundle, bonds and bond providers,
boundary, symmetry, query; this run also reviews `model/periodic.ts` and `Model.withAtomicConformation`, new in the
path since v2). HEAD `16277d822`, v2 base `0944823d4`, Node 26.12, Apple M4. Read-only.

Measurement conditions: the machine was overloaded throughout (load 23–31, swap 17.2 of 18.4 GB used, ~13 other
auditors). Wall times swung up to 5× between identical runs, so the A/B numbers below are **CPU time**
(`process.cpuUsage`, same process layout, base and HEAD interleaved twice) and **allocation**
(`total_allocated_bytes`); wall times are given only as ranges or ratios. A/B base: `git archive 0944823d4`
(`pa3/00/base`), bundled with the same bench through an alias plugin. Fixtures: `7qpd`, `flat-32` (186k atoms, 416
units), `asm-32` / `asm-192` (1tup assembly, 186k / 1.12M atoms), `flat-192` (single stages only), GRO
`split-peptide.gro` and the synthetic solvated system of `periodic-scaling.test.ts`. No asm-864 / flat-864 runs.

PRs since v2 touching the slice: #2016 (bond providers, new code), #2025 (05/P14), #2026, #2032 (05/B9, 05/P15
structure part), #2033, #2034 (05/B6, B10, B12), #2037 (periodic trajectories, new code), #2042 (loci/stats reused
maps), #2044 (05/B2), #2045 (05/P17, P19, P20 part), #2047 (05/B4, B5, B8, B11, P11, P16). Fixes in core that move this
slice's numbers: #2025's 02/P10 (Vec3 elements kind) and #2046.

## Hot paths now

| Path | Cost (fixture) | Where the time goes |
| --- | --- | --- |
| `Structure.interUnitBonds` | asm-192 364–369 ms CPU / 36–39 MB (base 539–577 ms / 95 MB); asm-32 113–129 ms (base 162–195); flat-32 87–98 ms (base 123–149) | B-grid query 42 % total, `getCandidates` 15 %, `findPairBonds` self 9 %; every symmetry-equivalent pair still recomputed (P13: 30,906 pairs in 241 classes on asm-192) |
| `Unit.bonds` (intra) | flat-32 272–336 ms CPU (base 348–406), 32–37 MB | per-atom 3.1 Å grid query 39 % self, `findBonds` 12 %; each pair found from both atoms (P21) |
| `lookup+boundary` | flat-32 60–95 ms CPU, 5–16 MB (base 120–144 ms, 41 MB); asm-192 37–41 ms, 6 MB (base 46–48 ms, 51 MB) | flat: 98-direction extrema per unit, `computeSphereExtrema` 26 % + `Vec3.dot` 13 % + `includePositionRadius` 6 % self (P2) |
| First move of a structure (P14 fix) | asm-192: small component's first `interUnitBonds` computes the root's: ~35× its own (P24) | `Structure.interUnitBonds` → `parent.interUnitBonds` |
| Interactive queries | 'Covalently Bonded Component' on 7qpd: 2.2 s loaded (v2 0.44 s quiet); 200 → 400 → 800 layers 52 → 265 → 1140 ms | fixed-point `includeConnected`, quadratic (P18) |
| Periodic frames | 16 frames × (frame + whole + box) kept by a module-level map: 191 MB at 250k atoms after every caller dropped them | `ModelPeriodicity.recent` (P23) |

## v2 → now

Counts: 27 v2 items (P2, P10–P22, B2–B14): **fixed 14, partly fixed 2, still open 11, regressed 0, dropped 0**.
New costs or bugs from fixes since v2: P24 (from #2025's P14 fix), B15 (from #2037). No regression in the A/B: every
stage's CPU time and allocation went down or stayed.

| v2 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| P2 | Unit boundaries at 98 EPOS directions up to 10k atoms | still open | | `unit.ts` and core `getBoundary` unchanged; extrema still 45 % self of the flat-32 stage. Allocation fell 41 → 5–16 MB through 02/P10 (#2025), not through P2 |
| P10 | One `Structure` per atom in singleton selections | still open | | [selection.ts:171](../../packages/model/src/model/structure/query/selection.ts#L171) unchanged |
| P11 | `ElementSetIntraBondCache` never evicts | fixed | #2047 | LRU of 64 sets, and 4 provider keys besides `default` (#2016) |
| P12 | `IntraUnitBondMapping`: 3 × Uint32 per directed bond instance | still open | | [mapping.ts:74](../../packages/model/src/model/structure/structure/mapping.ts#L74) unchanged |
| P13 | Inter-unit bonds per symmetry-equivalent unit pair | still open | | asm-192 prototype rerun: 642 → 157 ms (loaded), same 21,120 edges, 241 classes |
| P14 | `Structure.transform` and its components recompute inter-unit bonds | fixed | #2025 | asm-192: transformed 1.5–1.6 ms, polymer component 10–48 ms vs subset under the root 9–33 ms, same edges. Side cost: P24 |
| P15 | Structure and loci boundaries allocate transformed extrema | partly fixed | #2032 | structure: `includeTransformed` / `radiusTransformed` (asm-192 lb 51 → 6 MB); loci: [loci.ts:661](../../packages/model/src/model/structure/structure/element/loci.ts#L661) still `Sphere3D.transform`, 32–35 MB per call (below) |
| P16 | B-grid query at 4 Å for distance bonds | fixed | #2047 | `MaxPairingThreshold` = 2.954 Å, exact (same expressions as `getPairingThreshold`) |
| P17 | `hasCommonStructConnRecord` scans per pair | fixed | #2045 | partners per (struct_conn, elements) in a `WeakMap`; gone from the profile (was 9–13 % self) |
| P18 | Fixed-point `includeConnected` is quadratic | still open | | [modifiers.ts:358](../../packages/model/src/model/structure/query/queries/modifiers.ts#L358) unchanged; 2× layers → 4–5× time |
| P19 | Atom-group queries copy whole chains | fixed | #2045 | `addChainRanges` sets whole units; asm-192 polymer query 2 MB, 960 / 960 parent units reused |
| P20 | Whole-unit children rebuild every edge and the eager indices | partly fixed | #2045 | `wholeUnitBonds` shares the pair lists; the `InterUnitGraph` constructor still builds `edges`, `edgeKeyIndex`, `vertexKeyIndex` eagerly: 9.8 MB for 21k edges, 7.8 MB for the asm-192 polymer child |
| P21 | Intra pairs found twice by per-atom sphere queries | still open | | rerun `pairs-half.ts` flat-32: per-atom 349 ms vs half stencil 119 ms (loaded; v2 82 vs 29), same 563,456 pairs |
| P22 | `within` with `min-radius` is all-pairs | still open | | [structure-distance.ts:25](../../packages/model/src/model/structure/query/utils/structure-distance.ts#L25) unchanged |
| B2 | `entityIndices` uses atomic chain segments for coarse units | fixed | #2044 | coarse `chainElementSegments` |
| B3 | Bundle / Loci-expression grouping on a 32-bit hash | still open | | [bundle.ts:119](../../packages/model/src/model/structure/structure/element/bundle.ts#L119), [loci.ts:815](../../packages/model/src/model/structure/structure/element/loci.ts#L815) |
| B4 | `handleUnitChainsPartitioned` never advances `i` | fixed | #2047 | |
| B5 | `getStructureDeduplicate` dedups `this.currentUnit` | fixed | #2047 | |
| B6 | `pick` push element / pop structure; inverted `i % N` | fixed | #2034 | also `getSelectionProperties`, `merge` |
| B7 | `Bundle.toLoci` / `toStructure` create an Error, never throw | still open | | [bundle.ts:156](../../packages/model/src/model/structure/structure/element/bundle.ts#L156), :227 |
| B8 | `findMatesRadius`: world-space centre in a unit-local grid | fixed | #2047 | |
| B9 | `remapModel(m, parent)` swaps in the root's untransformed unit | fixed | #2032 | `remapx.ts` rerun: remapped component x 96.36 (= source) |
| B10 | `Bundle.toStructure` maps every unit through the first unit's elements | fixed | #2034 | `b10.ts` rerun: 20 atoms, copy 2 = 10..19 |
| B11 | `Loci.toExpression` merges models sharing an operator name | fixed | #2047 | `expr.ts` rerun: round trip exact |
| B12 | `includeSurroundings` pre-filters units by `radius` | fixed | #2034 | |
| B13 | Static component remap trusts the `(size, min, max, mid)` element hash | still open | | [structure.ts:1192](../../packages/model/src/model/structure/structure/structure.ts#L1192); guard in `structure-component.ts:200` unchanged |
| B14 | `within` adds the query atom's vdW radius, ignores `atom-radius` | still open | | `b14.ts` rerun: identical (126 vs 42 centre / 214 surface) |

Fixes since v2 reviewed hunk by hunk (bonds/subset/provider, structure getter, unit `bonds`, boundary, lookup3d,
symmetry, bundle, loci, stats, generators, filters, combinators, subset-builder, unit-merging, model/periodic):

- **#2025 (P14)**: correct (rigid-motion test from the first moved pair, key check for index pairs), but a child of a
  moved structure whose root has no bonds computes the **root's** bonds: **P24** below.
- **#2037**: `Model.withAtomicConformation` copies the frame's whole `_dynamicPropertyData` except the centre:
  **B15**; the made-frame cache is a module-level strong `Map`: **P23**.
- **#2016**: provider resolution, cache keys and `childBonds: 'root'` checked; `noDistance` keeps struct_conn, index
  pairs and dictionary bonds as documented; no bug found. `modelStructure()` has no in-repo caller (B15 makes it stale
  on periodic frames).
- #2045's `wholeUnitBonds` shares the parent's `UnitPairEdges` objects: safe (edges are by unit id and index, nothing
  mutates them). #2047's `MaxPairingThreshold` is exact. #2034's B10 fix adds a `SortedArray.areEqual` per copy only
  when the arrays differ (copies of a root share them). #2047's LRU: an index-pair model with > 64 units thrashes it on
  a structure rebuild, but the miss is a cheap index-pair lookup: not a finding.

## Open performance findings (ranked)

| Rank | Id | Status | Title | File | Kind | Estimate | Effort | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | P13 | still open | Inter-unit bonds recomputed for every symmetry-equivalent unit pair | [inter-compute.ts:406](../../packages/model/src/model/structure/structure/unit/bonds/inter-compute.ts#L406) | time | measured (prototype): asm-192 ÷4 (642 → 157 ms loaded; v2 quiet 172 → 45); asm-864 v2 774 → 131 ms; flat unchanged | M | no |
| 2 | P18 | still open | Fixed-point `includeConnected` is quadratic | [modifiers.ts:358](../../packages/model/src/model/structure/query/queries/modifiers.ts#L358) | time | measured: 7qpd 0.44 s quiet (2.2 s loaded) → a few ms estimated; quadratic in layers | M | no |
| 3 | P23 | new (from #2037); **fixed in #2050** (`WeakMap` by hierarchy, LRU per topology) | Periodic frames made whole are held by a module-level map | [periodic.ts:65](../../packages/model/src/model/structure/model/periodic.ts#L65) | memory | measured: 191 MB retained at 250k atoms after every reference is dropped; ~0.75 GB at 1M (estimated, linear) | S | no |
| 4 | P2 | still open | Unit boundaries use 98 directions up to 10k atoms | [unit.ts:435](../../packages/model/src/model/structure/structure/unit.ts#L435) | time | v2 estimate stands: `lookup+boundary` flat-192 ~÷3; extrema 45 % self of the flat-32 stage now | S | no |
| 5 | P21 | still open | Intra `findBonds` finds every pair twice through per-atom queries | [intra-compute.ts:387](../../packages/model/src/model/structure/structure/unit/bonds/intra-compute.ts#L387) | time | measured: pair enumeration ÷2.9 (349 → 119 ms loaded, flat-32); bonds-intra −25–40 % (v2 estimate) | M | no |
| 6 | P24 | new (from #2025) | First `interUnitBonds` of a small component of a moved structure computes the root's | [structure.ts:256](../../packages/model/src/model/structure/structure/structure.ts#L256) | time | measured: asm-192 non-polymer component 26–39 ms own vs 1.1–1.8 s via the root (loaded), once per root; asm-864 ~0.8 s quiet (v2 root cost) | S | no |
| 7 | P15 | partly fixed | Loci boundaries allocate transformed extrema per whole unit | [loci.ts:661](../../packages/model/src/model/structure/structure/element/loci.ts#L661) | memory | measured: asm-192 loci of 1,727 whole units 32–35 → 0.3–1.6 MB per call, same sphere; time equal with GC outside | S | no |
| 8 | P20 | partly fixed | `InterUnitGraph` builds its edge and vertex indices eagerly | [inter-unit-graph.ts:62](../../packages/core/src/math/graph/inter-unit-graph.ts#L62) | time+memory | measured: 9.8 MB / 5–12 ms per 21k-edge graph, paid by every child structure's bonds | S | no |
| 9 | P10 | still open | One `Structure` per atom in singleton selections | [selection.ts:171](../../packages/model/src/model/structure/query/selection.ts#L171) | time | v2: −19 ms per 186k singletons | M | no |
| 10 | P22 | still open | `within` with `min-radius` is an all-pairs scan | [structure-distance.ts:25](../../packages/model/src/model/structure/query/utils/structure-distance.ts#L25) | time | v2: O(\|a\|·\|b\|) → lookup; ties to B14 | M | no |
| 11 | P12 | still open | `IntraUnitBondMapping`: 3 × Uint32 per directed bond instance | [mapping.ts:74](../../packages/model/src/model/structure/structure/mapping.ts#L74) | memory | v2: −117 MB at asm-864 when used | M | no |

### P13: inter-unit bonds per symmetry-equivalent pair

- Unchanged mechanism (v2). `pairs.ts` asm-192: 1,728 units, 30,906 unit pairs (19,028 with a water or ion), 241
  classes (invariant, size, relative operator to 1e-3). `symbonds.ts` (crude class cache, two-unit structure per class):
  642 → 157 ms, edges 21,120 = 21,120. With P16/P17 done, the B-grid query is 42 % of the stage, so the class cache now
  removes most of what is left.
- Fix as in v2 (key by invariant ids, element identity and a rounded relative matrix, confirm with `Mat4.areEqual`,
  replay index pairs per class); a `Structure.eachUnitPairClass` would serve slice 06's interactions too.

### P18: fixed-point `includeConnected`

- Code unchanged. `conn.ts` 7qpd, seeded mid-unit: fixed point 2,185 ms loaded for 5,334 atoms; 50 / 100 / 200 / 400 /
  800 layers 3 / 9 / 52 / 265 / 1,140 ms. Each step rebuilds a subset structure over every atom reached. Fix as v2:
  frontier BFS with a visited mark, one subset at the end.

### P23 (new): periodic frame cache retains models

- `ModelPeriodicity.recent` is a module-level `Map<Model, Map<mode, Model>>` with 16 source frames, each with up to four
  made models (whole, box, box around a selection, no jump), every one with its own x/y/z. It is never cleared, so after
  the trajectory is removed from the state, its last 16 frames and their made copies stay until 16 other frames are
  made, possibly never in the session.
- `retain.ts` (synthetic 250k-atom solvated system, 16 frames, whole + box): heap + buffers 52 MB with the model,
  199 MB with the frames made, **199 MB after every reference is dropped and two GCs**.
- Fix: keep the LRU per topology, `WeakMap<AtomicHierarchy, LRU>` (frames of a topology share the hierarchy; the
  ephemeron lets made models point back at it), or clear it from the trajectory's dispose. Keeps the playback hit rate.

### P2, P21: unit boundaries and intra pairs

- Both unchanged in code. flat-32 profile, `lookup+boundary` only: `computeSphereExtrema` 25.7 %, `Vec3.dot` 12.9 %,
  `includePositionRadius` 5.7 %, `scaleAndSub` 4.8 % self. `bonds-intra` only: `GridLookup3D.query` 39.2 %,
  `findBonds` 11.7 %, `Result.add` 4.8 % self. Fixes as in v2 (14 directions for unit spheres with a padded structure
  sphere; a half-stencil pair walk in core).

### P24 (new, from #2025): the first move computes the root's inter-unit bonds

- `Structure.interUnitBonds` ([structure.ts:249](../../packages/model/src/model/structure/structure/structure.ts#L249)):
  when the parent has no bonds and the child is a moved image (`isMovedImage`), it computes `parent.interUnitBonds`
  (every unit of the root) and restricts them. The root is the untransformed structure (parents are flattened), so a
  ligand or water component of a moved structure pays for the whole root's bonds when the root never needed them (the
  default preset draws the polymer as cartoon, so only small components ask for inter-unit bonds).
- `xform2.ts` asm-192, fresh root without bonds, `Structure.transform`, non-polymer component (576 of 1,728 units):
  own computation 26–39 ms vs 1,090–1,830 ms through the root (loaded), 0 edges either way. Once per root: later
  moves reuse the root's bonds, which is the point of the fix.
- Fix: take the root route only when the child is a large part of the root (e.g. `child.elementCount * 4 >=
  parent.elementCount`, or the root already has bonds); otherwise compute the child's own, as before #2025.

### P15 (loci part), P20

- `locibnd.ts` asm-192: `Loci.getBoundingSphere` over 1,727 whole transformed units 15–26 ms, 32–35 MB per call;
  transforming each extremum into one scratch `Vec3` (the structure path's `includeTransformed`) 14–24 ms, 0.3–1.6 MB,
  radius 366.4816 both. Reached by focus and camera fit on large selections. Fix: reuse `includeTransformed` /
  `radiusTransformed` from `util/boundary.ts` in `getWholeUnitBoundarySphere`'s two callers.
- P20: `iub.ts` asm-192: `new InterUnitBonds(map)` 5–12 ms, 9.8 MB for 21,120 edges; `subsetInterUnitBonds` of the
  whole-unit polymer child 8–13 ms, 7.8 MB, nearly all of it the constructor. Fix: lazy `edges` / indices in core.

## Bugs

| Id | Severity | Status | Title | File | Failure scenario |
| --- | --- | --- | --- | --- | --- |
| B15 | medium | new (from #2037) | `Model.withAtomicConformation` copies coordinate-derived dynamic data | [model.ts:299](../../packages/model/src/model/structure/model/model.ts#L299) | periodic 'whole' shown with a tubular-helix cartoon, switched to 'box': helices of molecules moved into the box keep the whole frame's helix centres, 30 Å off (`stale.ts`) |
| B14 | medium | still open | `within` without `min-radius` adds the query atom's vdW, ignores `atom-radius` | [filters.ts:190](../../packages/model/src/model/structure/query/queries/filters.ts#L190) | PyMOL `around` / `within` over-select (7qpd, 3 Å: 126 atoms vs 42 by centre distance); `gap` misses the target radius (126 vs 214) |
| B7 | low | still open | `Bundle.toLoci` / `toStructure` create an Error, never throw | [bundle.ts:156](../../packages/model/src/model/structure/structure/element/bundle.ts#L156) | a bundle applied to another structure silently selects wrong atoms |
| B13 | low | still open (plausible) | Static component remap trusts the weak element hash | [structure.ts:1192](../../packages/model/src/model/structure/structure/structure.ts#L1192) | `[1,3,5,9]` vs `[1,4,5,9]` collide: a static component below a changing source stays stale |
| B3 | low | still open (plausible) | Bundle / Loci-expression grouping on a 32-bit hash without equality | [bundle.ts:119](../../packages/model/src/model/structure/structure/element/bundle.ts#L119) | rare merged groups in stored selections (P ≈ n²/2³³) |

### B15: stale dynamic data on periodic frames

- `withAtomicConformation` builds the made frame's `_dynamicPropertyData` with `Object.assign` from the source frame,
  deleting only the centre. Every made frame (whole, box, box around a selection, no jump) therefore starts with what was
  computed on the frame it came from, including coordinate-derived values: dynamic `CustomModelProperty` containers
  (`HelixOrientationProvider`, `type: 'dynamic'`), `BondProviderModelStructure` (a `Structure` of the old model, with
  bonds from the old coordinates), and the bond caches (harmless: only coordinate-free index-pair bonds are cached).
  `CustomModelProperty.attach` returns early when a value exists with the same params, so it is never recomputed.
- Failure (`stale.ts`, `split-peptide.gro`): frame 1 made whole, helix orientation attached (as the cartoon does on the
  model it shows), then `apply(frame, InBox)`: the peptide moves 30.00 Å into the box; the box model's helix
  orientation **is the whole model's object**, centre 0 at (35.09, 13.87, 9.00) vs (5.09, 13.87, 9.00) recomputed. In the
  app: the cartoon's `tubularHelices` reads `HelixOrientationProvider.get(unit.model)` (trace-iterator.ts:472), so
  after switching the periodic mode from 'whole' (the default) to 'box', or from 'off' to 'whole', tubular helices of
  moved or rejoined molecules are drawn at the previous mode's positions. A provider reading `modelStructure()` would
  read the previous frame's units and inter-unit bonds (latent: no in-repo caller).
- Fix: start the made frame's dynamic data empty and copy only what a reader set (trajectory info, periodic box,
  unwrapped conformation, model symmetry, index-pair bonds, coarse-grained), as `_trajectoryFromModelAndCoordinates`
  does for frames; or copy all but the keys of dynamic custom model properties and `BondProviderModelStructure`. Test:
  attach helix orientation to a whole frame, apply box, expect a new value.

### B14, B7, B13, B3

- B14: `b14.ts` rerun on 7qpd (30-atom target, 3 Å, per-atom singletons): lookup path 126 with or without
  `atom-radius`; centre distance 42; vdW surface distance 214 (= the `min-radius: 0` path). Code unchanged
  (`check(pos, maxRadius + c.r(e))`). Fix as v2 (drop `c.r(e)` for atomic units; with `elementRadius`, query at
  `maxRadius + r + maxTargetRadius` and test each hit exactly, which also fixes P22).
- B7, B13, B3: code unchanged since v2; scenarios and fixes as in v2.

## Already good, re-checked

- `subsetInterUnitBonds` with the rigid-motion test: a component of a moved structure restricts the root's bonds;
  identical edge sets on asm-192 (21,120 / 16,512).
- `hasCommonStructConnRecord` (P17), `MaxPairingThreshold` (P16), `addChainRanges` (P19), `structureUnion`,
  `eachUnitPair` itself (11 ms for 30,906 pairs on asm-192, no-op callback).
- Periodic per-frame costs (plan, whole, box) are linear and covered by `periodic-scaling.test.ts`; not re-measured.

## Not audited / not measured

- Quiet-machine wall times: the machine was overloaded and swapping; v2's quiet numbers remain the reference for
  absolute times, this run's ratios are CPU-time A/B.
- asm-864 / flat-864 (memory), flat-192 beyond one profiled run (13.7 s wall for bonds-intra, swap-bound).
- No-jump seeks in `periodic.ts` (each walked frame makes a whole model and enters the 16-frame cache); plausible
  cost, no fixture with many frames.

## Bench scripts

`/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/31ee827f-4715-4967-bec3-085f5013d11f/scratchpad/pa3/05/`.
`node build.mjs X.ts X.mjs [absolute base root]` (esbuild, `molstar-src`; with a base root, `@molstar/{core,io,model,graphics}`
resolve into it), run with `node --expose-gc --max-old-space-size=3000 X.mjs`. `bash ab.sh <fixture> <runs> <rounds>
[stages]` interleaves base and HEAD.

| Script | Measures |
| --- | --- |
| `ab.ts`, `ab.sh` | lookup+boundary, bonds-intra, bonds-inter on a fresh structure per run: CPU, wall, allocation, A/B |
| `prof.mjs` | self / total % of a `--cpu-prof` profile, source-mapped |
| `pairs.ts`, `symbonds.ts` | unit pairs and symmetry classes; class-cached inter-unit bonds (P13) |
| `xform.ts`, `xform2.ts` | transformed structure and component bonds (P14); first-move cost via the root (P24) |
| `iub.ts` | polymer component query, whole-unit subset, `InterUnitGraph` constructor (P19, P20) |
| `locibnd.ts` | loci bounding sphere over whole transformed units, current vs in place (P15) |
| `conn.ts` | fixed-point and layered `includeConnected` (P18) |
| `verify/pairs-half.ts` | per-atom queries vs half-stencil pair walk (P21) |
| `retain.ts` | periodic made-frame retention after dropping references (P23) |
| `stale.ts` | helix orientation of a box frame vs the whole frame it came from (B15) |
| `verify/b14.ts`, `verify/b10.ts`, `remapx.ts`, `expr.ts` | B14 (open), B10, B9, B11 (fixed) |
