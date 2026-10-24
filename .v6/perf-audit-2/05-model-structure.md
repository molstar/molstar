# Slice 05: model-structure (v2)

Scope: `packages/model/src/model/structure` (`Structure`, `Unit`, element/Loci/Bundle, bonds, boundary, symmetry,
query). Measured at `0944823d4` (= master), Node 25.6, Apple M4, on a shared, moderately loaded machine (±10–20 %).
Fixtures as in [00-profile.md](00-profile.md): `7qpd`, `1grm` (5 models), `flat-32` (= v1 big-bcif, 186k atoms, 416
units), `asm-32/192/864` (1tup assembly, 186k / 1.12M / 5.04M atoms, 9 symmetry groups), `flat-192/864` (unshared).
Read-only: no fix below is applied. MVS loading and animation is covered by
[.v6/mvs-optimizations/research](../mvs-optimizations/research/); nothing here touches frozen gl or Canvas3D code.

Since v1: #1994 (P2 structure/loci part, P3, P4, P5, P6), #2007 (P7), #2018 (P1, P9), #2022 (B1). #2000 is not an
ancestor of HEAD, but #1994 carries its P6 change. Two fixes brought new bugs: B9 (from #2022's B1 fix) and B10
(from #1994's P3 fix).

## Hot paths now

| Path | Cost (fixture) | Where the time goes |
| --- | --- | --- |
| `Structure.interUnitBonds` | 0.18 s asm-192, 0.78 s / 443 MB asm-864, 0.2 s flat-192 | B-grid query 45–52 %, `hasCommonStructConnRecord` 9–13 %; on assemblies every symmetry-equivalent pair redone (P13) |
| `Unit.bonds` (intra) | 0.82 s flat-192, 3.6 s flat-864 | per-atom grid query 55 %, each pair found from both atoms (P21) |
| `lookup+boundary` | 0.3 s flat-192, 171 ms asm-864 | flat: per-unit 98-direction `getBoundary` 83 % (P2); asm: `Sphere3D.transform` allocations, 79 % of the stage on asm-192 (P15) |
| Transformed structures and their components | flat-192 270 + 212 ms per component, asm-192 185 + 150 ms | inter-unit bonds recomputed after `Structure.transform` (P14) |
| Interactive queries | 'Covalently Bonded Component' 0.44 s on 7qpd | fixed-point `includeConnected`, quadratic (P18) |

Already good, re-checked: `structureUnion` (#2018), the carbohydrate BFS (#2018), `subsetIntraUnitBonds` (O(child)),
`Loci.size` and whole-unit shortcuts, `eachUnitPair` itself (25 ms for 155k pairs), symmetry building (12 ms at
asm-864), symmetry mates. v1's R3 (3k7f `lookup+boundary` ×2) does not reproduce warm: v1 base and HEAD both 1.45–1.51
+ 0.21 + 0.11–0.13 ms over 20 runs.

## v1 → now

Counts: 21 v1 items (P1–P12, B1–B8, and the 13-large-scale rows of this slice): fixed 8, partly fixed 2, still open 10,
regressed 0, dropped 1. New bugs from fixes: B9 (#2022), B10 (#1994).

| v1 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| P1 | `structureUnion` merges one structure at a time | fixed | #2018 | [structure-set.ts:28](../../packages/model/src/model/structure/query/utils/structure-set.ts#L28): each unit's sets merged once (union / concat / sort / marks), full-unit shortcut kept |
| P2 | 98-direction EPOS for structure, loci and unit boundaries | partly fixed | #1994 | structure and loci switch to 14 directions above 10k points; `Unit.boundary` still 98 up to 10k atoms (P2 below) |
| P3 | `Bundle.toStructure` / `Loci.toStructure` lose symmetry sharing | fixed | #1994 | child built once per group, `applyOperator` for copies; `Loci.toStructure` via `UnitEquivalenceBuilder`. Side effect: B10 |
| P4 | `findPairBonds` scans every atom of A per unit pair | fixed | #1994 | `getCandidates` takes A's atoms from A's grid ([inter-compute.ts:71](../../packages/model/src/model/structure/structure/unit/bonds/inter-compute.ts#L71)); asm-864 2.1 s / 3.6 GB → 0.78 s / 443 MB |
| P5 | Child structures recompute inter-unit bonds | fixed | #1994 | `subsetInterUnitBonds` when bond props and operators match; polymer component ÷11–÷39. Not done: (b) return the parent; misses: P14 (transformed), P20 (whole-unit children) |
| P6 | Child units recompute intra bonds | fixed | #1994 (content of #2000) | `subsetIntraUnitBonds`, O(child); the parent's getters are never called |
| P7 | Intra `findBonds` queries 4 Å, string lookups | fixed | #2007 | 3.1 Å query, dictionary pairs through the residue, Int8 element index, struct_conn residue set; flat-192 1490 → 821 ms |
| P8 | Index-pair bonds with operators: per-unit scans | dropped | | unreachable: no reader sets `operatorA/B` (unchanged) |
| P9 | Carbohydrate links: all-pairs loop, BFS per pair | fixed | #2018 | one depth-3 BFS per ring; elements are contiguous per unit |
| P10 | One `Structure` per atom in singleton selections; surroundings | still open | | [selection.ts:171](../../packages/model/src/model/structure/query/selection.ts#L171), [modifiers.ts:94](../../packages/model/src/model/structure/query/queries/modifiers.ts#L94) |
| P11 | `ElementSetIntraBondCache` never evicts | still open | | inflow smaller: restricted child bonds are not cached |
| P12 | `IntraUnitBondMapping`: 3 × Uint32 per directed bond instance | still open | | unchanged; ~117 MB at asm-864 when a structure-level bond visual asks |
| B1 | `remapModel` builds a private root per component | fixed | #2022 | `remapModel(m, parent)`, components share `a.root`. The whole-unit sharing introduced B9 |
| B2 | `entityIndices` uses atomic chain segments for coarse units | still open | | [structure.ts:561](../../packages/model/src/model/structure/structure/structure.ts#L561) |
| B3 | Bundle / Loci-expression grouping on a 32-bit hash | still open | | [bundle.ts:119](../../packages/model/src/model/structure/structure/element/bundle.ts#L119), [loci.ts:805](../../packages/model/src/model/structure/structure/element/loci.ts#L805) |
| B4 | `handleUnitChainsPartitioned`: `const i = 0` loop | still open | | [stats.ts:289](../../packages/model/src/model/structure/structure/element/stats.ts#L289) |
| B5 | `getStructureDeduplicate` dedups `this.currentUnit` | still open | | [subset-builder.ts:80](../../packages/model/src/model/structure/structure/util/subset-builder.ts#L80), no callers |
| B6 | `pick` push element / pop structure; inverted `i % N` | still open | | [filters.ts:24](../../packages/model/src/model/structure/query/queries/filters.ts#L24); `i % N` at :28, :90, :114, :128, combinators.ts:24 |
| B7 | `Bundle.toLoci` / `toStructure`: `new Error` never thrown | still open | | [bundle.ts:156](../../packages/model/src/model/structure/structure/element/bundle.ts#L156), :193 |
| B8 | `findMatesRadius`: world-space centre in a unit-local grid | still open | | [symmetry.ts:332](../../packages/model/src/model/structure/structure/symmetry.ts#L332), latent |
| 13/P2–P7 | 13-large-scale rows of this slice | partly fixed | #1994, #2007 | same ids as above; asm-864 bonds-inter 784 ms (was 2.1 s), polymer component 190 ms (was 2 s), flat-864 bonds-intra 3.6 s (was 7.3 s); P2's unit part open |

## Open performance findings (ranked)

| Rank | Id | vs v1 | Title | File | Kind | Impact | Estimate | Effort | Verdict | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | P13 | new | Inter-unit bonds recomputed for every symmetry-equivalent unit pair | [inter-compute.ts:403](../../packages/model/src/model/structure/structure/unit/bonds/inter-compute.ts#L403) | time+memory | high | asm-864 774 → ~130 ms, asm-192 172 → 45 ms (prototype); flat unchanged | M | confirmed | no |
| 2 | P14 | new; **fixed** #2025 | `Structure.transform` and its components recompute inter-unit bonds | [structure.ts:912](../../packages/model/src/model/structure/structure/structure.ts#L912) | time | high | per transform update: asm-192 185 ms + 150 ms per component → ~5 ms; flat-192 ~0.5 s → <5 ms; **measured** (`transform*` stages, loaded machine): asm-192 moved 320–730 → 0.6–1.8 ms, polymer component 310–540 → 17–38 ms; flat-32 57–99 → 0.3–1.4 ms and 52–88 → 0.4–0.6 ms; the root computes its bonds once if it has none | S | confirmed | no |
| 3 | P18 | new | Fixed-point `includeConnected` ('Covalently Bonded Component') is quadratic | [modifiers.ts:358](../../packages/model/src/model/structure/query/queries/modifiers.ts#L358) | time | high | 7qpd 438 ms → a few ms; long connected chains from tens of s to tens of ms | M | confirmed | no |
| 4 | P2 | partly fixed | Unit boundaries still use 98 EPOS directions up to 10k atoms | [unit.ts:418](../../packages/model/src/model/structure/structure/unit.ts#L418) | time | medium | `lookup+boundary` flat-192 ~300 → ~100 ms, flat-864 1.26 → ~0.4 s | S | confirmed | no |
| 5 | P21 | new | Intra `findBonds` finds every pair twice through per-atom sphere queries | [intra-compute.ts:396](../../packages/model/src/model/structure/structure/unit/bonds/intra-compute.ts#L396) | time | medium | bonds-intra −25–40 % (neighbour enumeration 82 → 29 ms on flat-32) | M | confirmed | no |
| 6 | P15 | new; partly fixed #2032 (structure boundary; loci open) | Structure and loci boundaries allocate transformed extrema per unit | [boundary.ts:42](../../packages/model/src/model/structure/structure/util/boundary.ts#L42) | time+memory | medium | asm-192 18 → 9.6 ms, asm-864 117 → 50 ms; −50 % allocation per call | S | confirmed | no |
| 7 | P17 | new; **fixed** #2045 | `hasCommonStructConnRecord` scans the smaller unit for every water/ion pair | [inter-compute.ts:449](../../packages/model/src/model/structure/structure/unit/bonds/inter-compute.ts#L449) | time | low | bonds-inter −9–13 % | S | confirmed | no |
| 8 | P16 | new; **fixed** #2047 | `findPairBonds` queries B's grid at 4 Å for distance bonds that end at 2.95 Å | [inter-compute.ts:226](../../packages/model/src/model/structure/structure/unit/bonds/inter-compute.ts#L226) | time | low | bonds-inter ≤ −15–20 % (upper bound) | S | plausible | no |
| 9 | P19 | new; **fixed** #2045 | Atom-group queries copy whole chains into JS arrays, then reuse the parent unit | [generators.ts:158](../../packages/model/src/model/structure/query/queries/generators.ts#L158) | time+memory | low | asm-864 polymer component query 25 ms / 148 MB → ~2 ms / ~0 | S | confirmed | no |
| 10 | P20 | new; partly fixed #2045 (a child of whole units shares the pair lists; the eager graph indices remain) | Whole-unit children rebuild every inter-unit edge and the eager graph indices | [subset.ts:170](../../packages/model/src/model/structure/structure/unit/bonds/subset.ts#L170) | time+memory | low | per component on asm-864 −20 ms / −60 MB | S | confirmed | no |
| 11 | P10 | still open | One `Structure` per atom in singleton selections; surroundings per query atom | [selection.ts:171](../../packages/model/src/model/structure/query/selection.ts#L171) | time | low | −19 ms per 186k singletons; up to 2× on compact surroundings | M | plausible | no |
| 12 | P22 | new | `within` with `min-radius` is an all-pairs scan | [structure-distance.ts:155](../../packages/model/src/model/structure/query/utils/structure-distance.ts#L155) | time | low | hand-written queries only: O(\|a\|·\|b\|) → lookup | M | confirmed | no |
| 13 | P11 | still open; **fixed** #2047 | `ElementSetIntraBondCache` never evicts | [element-set-intra-bond-cache.ts:23](../../packages/model/src/model/structure/structure/unit/bonds/element-set-intra-bond-cache.ts#L23) | memory | low | bounds growth on index-pair models | S | confirmed | no |
| 14 | P12 | still open | `IntraUnitBondMapping`: 3 × Uint32 per directed bond instance | [mapping.ts:74](../../packages/model/src/model/structure/structure/mapping.ts#L74) | memory | low | −117 MB at asm-864 when used | M | confirmed | no |

### P13: inter-unit bonds per symmetry-equivalent pair

- Mechanism: `findBonds` runs `findPairBonds` for every pair from `eachUnitPair`. Every decision in it (candidates
  from `inv(opA)·opB`, the image transform `inv(opB)·opA`, struct_conn distances, alt-loc and occupancy rules,
  thresholds) depends only on the two invariant units (model, elements) and their relative rigid operator. Copies of an
  assembly repeat each class once per copy.
- Evidence (`pairs.ts`, `symbonds.ts`, verifier rerun): asm-864 154,916 pairs in 241 classes (invariant, size, relative
  matrix to 1e-3); asm-192 30,906 in 241; flat-192 49,514 in 49,514. A crude cache (string key per pair, a two-unit
  `Structure` + `computeInterUnitBonds` per class, edges replayed) gives identical edge counts (21,120 / 95,040):
  asm-192 172 → 45 ms, 101 → 77 MB; asm-864 774 → 131 ms, 445 → 276 MB.
- Fix: in `findBonds`, key by (A.invariantId, A.elements identity, B.invariantId, B.elements identity, hash of the
  relative matrix rounded ~1e-4), confirm with `Mat4.areEqual`, replay cached `Int32Array` index pairs + shared props
  through `startUnitPair/add/finishUnitPair`. Only for rigid operators; keep `validUnitPair` per pair. A generic
  `Structure.eachUnitPairClass` would serve the interactions' `addStructureContacts` too (slice 06, 9.8 s at asm-864).
- Risk: copies differing below the tolerance share bonds (more consistent than today). Edge order per pair unchanged.
  Pin edge sets on 1tup assemblies and 1cbs cells.

### P14: transformed structures recompute inter-unit bonds

- Mechanism: `Structure.transform` gives every unit `SymmetryOperator.create(name, M, old)` and creates the structure
  with `parent: s` and no `interUnitBonds`. The getter ([structure.ts:240](../../packages/model/src/model/structure/structure/structure.ts#L240))
  calls `subsetInterUnitBonds`, which returns `undefined` on the first operator mismatch (subset.ts:153–158), so it
  recomputes. Components of the transformed structure get the untransformed root as parent (flattened at
  structure.ts:758) and recompute too. A common rigid motion keeps every relative operator: the bonds are identical.
- Reached by `TransformStructureConformation` ([hierarchy.ts:206](../../packages/plugin/core/src/state/transforms/structure/hierarchy.ts#L206)):
  each Prism structure-move release (`structure-move/gizmos.ts` commit, #2021), superposition, transform param edits.
  Each component with a bond visual pays its own recompute.
- Evidence (`xform.ts`, verifier rerun): asm-192 transformed 185 ms / 122 MB, polymer component 152 ms / 103 MB vs
  5.1 ms subset under the root; 7qpd 3.3 + 2.5 ms vs 0.06 ms; flat-192 270 + 212 ms vs 0.48 ms (auditor). Same edge
  counts.
- Fix: (1) `Structure.transform` passes `interUnitBonds: s.state.interUnitBonds` when present (unit ids are kept);
  (2) `subsetInterUnitBonds` accepts a child whose operators are `M·parentOp` for one rigid `M` (from the first unit,
  others checked with `Mat4.areEqual` 1e-6). (2) also covers the components.
- Risk: rigid only (`transform` rejects others); `dynamicBonds` unaffected (distances preserved).

### P18: fixed-point `includeConnected`

- Mechanism: with `fixedPoint`, `includeConnectedStep` repeats until the count stops growing. Each step
  (`expandConnected`) builds a new `StructureUniqueSubsetBuilder` and `Structure` over all atoms reached so far,
  binary-searches each in the input unit and walks all its bonds: O(diameter × size).
- Evidence (`conn.ts`, verifier rerun, 7qpd, seeded at the middle atom of the 3,365-atom unit): fixed point 438 ms for
  5,334 atoms; 50 / 100 / 200 / 400 / 800 layers: 1.4 / 3.5 / 12.6 / 50 / 214 ms (2× layers → 4× time). Reached from
  the selection menu: 'Covalently Bonded Component' ([manipulate.ts:97](../../packages/plugin/core/src/state/queries/structure/manipulate.ts#L97)),
  on the main thread.
- Fix: frontier BFS per unit: a visited mark over the input unit's elements, a queue of (unit, UnitIndex), expand only
  new atoms through `offset/b` and `interUnitBonds.getConnectedUnits(unit).getEdges(index)`, `atomicBond.test` per edge
  as now; build the subset once. `layerCount` mode: one BFS level per layer, `getWholeResidues` per layer when
  `wholeResidues`.
- Risk: per-layer whole-residue semantics; pin against today's query on 7qpd and glycans.

### P2: unit boundaries at 98 directions (rest of v1 P2)

- Mechanism: `Atomic.boundary` → core `getBoundary`, 98 directions up to 10k atoms; with radii, 98 dot products plus a
  `scaleAndAdd` per atom and direction.
- Evidence: flat-192 `lookup+boundary` 304–358 ms / 303 MB, `getBoundary` 83 % total, `Vec3.dot` 37 % self,
  `computeSphereExtrema` 18 %; flat-864 1.26 s. v1 `bnd.ts`: unit spheres 33 → 9 ms on flat-32 with 14 directions; unit
  spheres stay exact (the radius pass visits every atom), but structure spheres built from 14-direction unit extrema miss
  atoms by up to 0.9 Å (v1 `bnd2.ts`).
- Fix: 14 directions with refine for unit spheres, plus a pad of the structure sphere by the largest unit extremum
  error, or keep 98 only below ~1k atoms. v1 found the radius pass over unit centre + radius 5–17 % too large.
- Risk: camera fit and snapshots shift slightly.

### P21: intra pairs found twice

- Mechanism: each atom runs a full 3.1 Å query (27 cells) and drops half the hits (`bI <= aI`).
- Evidence (verifier, `verify/pairs-half.ts`, flat-32): per-atom queries alone 82 ms for 563,456 pairs; a half-stencil
  cell list (grid built per unit, each pair once) 29 ms, same pair count; full `computeIntraUnitBonds` 122 ms. Profile
  flat-192: `GridLookup3D.query` 55 % self.
- Fix: a core `GridLookup3D` pair walk with a half stencil filling a reused CSR of larger-index neighbours with squared
  distances; run today's per-atom logic over it.
- Risk: edge order changes (unit.md §5 "edge order follows the grid"); edge-set pins exist in intra-compute.test.ts.

### P15: transformed extrema allocation

- Mechanism: `Sphere3D.transform` maps `extrema` to 98 new `Vec3` (and an array) per non-identity unit, twice per
  `computeStructureBoundary`; `Loci.getBoundingSphere` likewise for whole units.
- Evidence (`bnd2.ts` interleaved, no forced GC, verifier rerun): asm-192 17.8–18.5 → 9.6–9.7 ms, same sphere
  (r 366.8638); asm-864 117 → 50 ms (auditor). `bnd.ts` with GC forced outside the timing: 10.7 vs 14.0 ms, 52 → 24.5 MB
  allocated: the gain is GC pressure, not arithmetic.
- Fix: transform each extremum into one scratch `Vec3` before `includePosition` / `radiusPosition` (fall back to the
  transformed centre + radius when a sphere has no extrema), or `BoundaryHelper.includeSphereTransformed` in core.
- Risk: none (same order and arithmetic).

### Low

- **P17**: cache per unit elements, in a `WeakMap` inside `computeInterUnitBonds`, the struct_conn partners outside
  the unit (most have none), and per invariant pair on assemblies. Profile self 13 % flat-192, 9 % asm-192.
- **P16** (plausible): only the distance loop needs the B-grid query; at `Math.min(maxRadius, 3.1)` for non-coarse
  models (the coarse path never reaches that loop). The auditor's upper bound (3.1 Å everywhere, which also cut
  struct_conn bonds: 21,120 → 18,816 edges) was flat-192 197 → 157 ms, asm-864 804 → 625 ms; the real fix gains less.
  Keep `maxRadius` for candidates and struct_conn. After P13 the assembly gain is small.
- **P19**: polymer/protein/nucleic components go through `atomGroupsSegmented` at chain level; `addElementRange`
  pushes every index into a `number[]`, then `_getStructure` reuses the parent unit. `iub.ts`: asm-864 25 ms / 148 MB.
  Fix: `setUnit(unit.id, unit.elements)` when the passing chain segments cover the unit.
- **P20**: `subsetInterUnitBonds` returns the parent bonds only for an all-units child; a polymer child goes edge by
  edge through `InterUnitGraph.Builder` (asm-864 23 ms / 63 MB); the `InterUnitGraph` constructor builds `edges`,
  `edgeKeyIndex`, `vertexKeyIndex` eagerly (14 ms / 42 MB for 95k edges). Fix: filter the parent's `UnitPairEdges`
  lists when all maps are null; lazy indices in core. (Lazy indices do not help valence-model, whose 23 % in
  `getEdgeIndices` is the lookups themselves.)
- **P10** (plausible): singleton fast paths into one subset builder; adaptive reverse check for `includeSurroundings`
  (v1: 42 → 20.5 ms on 1tup-2x2x2).
- **P22**: `withinMaxRadius` / `withinMinMaxRadius` loop every target element per selection atom. Only reached with an
  explicit `min-radius` (no transpiler emits it). Fix with the target lookup at `maxDist + max radii`; ties to B14.
- **P11**: LRU cap (~64) or skip sets under ~1k elements.
- **P12**: derive (group, instance, edge) from `unitGroupOffset` and per-group edge counts in the three repr consumers
  (`bond-intra-unit-cylinder.ts:442`, `-line.ts:237`, `util/bond.ts:227`; repr, not gl).

## Bugs

| Id | vs v1 | Severity | Title | File | Verdict |
| --- | --- | --- | --- | --- | --- |
| B9 | new (from #2022); **fixed** #2032 | medium | `remapModel(m, parent)` swaps in the root's untransformed unit for whole units | [structure.ts:418](../../packages/model/src/model/structure/structure/structure.ts#L418) | confirmed |
| B14 | new | medium | `within` without `min-radius` adds the query atom's vdW radius and ignores `atom-radius` | [filters.ts:189](../../packages/model/src/model/structure/query/queries/filters.ts#L189) | confirmed |
| B10 | new (from #1994); **fixed** #2034 | low | `Bundle.toStructure` maps every grouped unit through the first unit's elements | [bundle.ts:242](../../packages/model/src/model/structure/structure/element/bundle.ts#L242) | confirmed |
| B12 | new; **fixed** #2034 | low | `includeSurroundings` with `atom-radius` pre-filters units by `radius`, not the query radius | [lookup3d.ts:228](../../packages/model/src/model/structure/structure/util/lookup3d.ts#L228) | confirmed |
| B11 | new; **fixed** #2047 | low | `Loci.toExpression` merges models that share an operator name | [loci.ts:774](../../packages/model/src/model/structure/structure/element/loci.ts#L774) | confirmed |
| B2 | still open; **fixed** #2044 | low | `entityIndices` uses atomic chain segments for coarse units | [structure.ts:561](../../packages/model/src/model/structure/structure/structure.ts#L561) | confirmed |
| B6 | still open; **fixed** #2034 | low | `pick` pushes an element, pops a structure; inverted `i % N` checks | [filters.ts:24](../../packages/model/src/model/structure/query/queries/filters.ts#L24) | confirmed |
| B7 | still open | low | `Bundle.toLoci` / `toStructure` create an Error, never throw | [bundle.ts:156](../../packages/model/src/model/structure/structure/element/bundle.ts#L156) | confirmed |
| B13 | new | low | Static component remap trusts the `(size, min, max, mid)` element hash | [structure.ts:1111](../../packages/model/src/model/structure/structure/structure.ts#L1111) | plausible |
| B3 | still open | low | Bundle / Loci-expression grouping on a 32-bit hash without equality | [bundle.ts:119](../../packages/model/src/model/structure/structure/element/bundle.ts#L119) | plausible |
| B8 | still open; **fixed** #2047 | low | `findMatesRadius`: world-space centre against a unit-local grid | [symmetry.ts:332](../../packages/model/src/model/structure/structure/symmetry.ts#L332) | confirmed (latent) |
| B4 | still open; **fixed** #2047 | low | `handleUnitChainsPartitioned` never advances `i` | [stats.ts:289](../../packages/model/src/model/structure/structure/element/stats.ts#L289) | confirmed (latent) |
| B5 | still open; **fixed** #2047 | low | `getStructureDeduplicate` dedups `this.currentUnit` | [subset-builder.ts:80](../../packages/model/src/model/structure/structure/util/subset-builder.ts#L80) | confirmed (latent) |

- **B9**: `ofParent(u)` returns `parent.unitMap.get(u.id)` when the elements array is the same, without comparing
  operators; `updateStructureComponent` passes `a.root` (structure-component.ts:211). Below
  `TransformStructureConformation`, `a.root` is the untransformed remapped root (parent flattened), and the transformed
  units share its elements arrays, so whole units of the component get the root's untransformed units. Failure
  (`remapx.ts`, verifier rerun, 1grm): `transform(s0, +100 Å x)`, whole-unit component; frame 1 as the plugin does:
  source x 96.36, remapped component x −3.64 (it is `s1`'s unit), 96.36 without the parent argument. In the app: a moved
  or superposed structure (Prism move, superposition) with static components, during trajectory playback the
  components jump back by the transform. Fix: also require `Unit.areOperatorsEqual(p, u)` in `ofParent` (and say so in
  structure.md §3, which states the rule without it); test with a transform between structure and component over two
  frames.
- **B14**: since ea612c3ac (Oct 2025, also 5.x) atomic `conformation.r` is the vdW radius. `withinMaxRadiusLookup`
  tests `targetLookup.check(pos, maxRadius + c.r(e))` and never reads `elementRadius`; the target's radius is never
  added. So plain `within` (PyMOL `around`, `within … of`, `near_to`, `beyond`; MolScript without `min-radius`)
  selects by centre distance + the query atom's vdW, against the documented `atom-radius` default 0, and PyMOL `gap`
  (`atom-radius: vdw`) misses the target radius. Failure (`verify/b14.ts`, 7qpd, 30-atom target, 3 Å, per-atom
  singletons): lookup path 126 atoms with or without `atom-radius`; centre distance gives 42, vdW surface distance 214
  (= the `min-radius: 0` all-pairs path). Fix: `check(pos, maxRadius)` without `elementRadius` (keep `c.r` for coarse
  units if intended); with it, query `maxRadius + r + maxTargetRadius` and test each hit exactly (P22's path). Changes
  results of existing queries (to the documented ones).
- **B10**: `toStructure` converts the group's UnitIndex set through `ue = _units[0].elements` for all units, while
  `toLoci` uses each unit's own. Copies of one invariant unit in a child structure keep `invariantId` (`getChild`) but
  may hold different elements; since #1994, when the first child is the whole first unit, the others are pushed whole
  without the `unit.elements === ue` guard. Failure (`verify/b10.ts`, asm-32): child with copy 1 = 10 atoms, copy 2 =
  20 atoms, loci UnitIndex 0..9 in both (20 atoms): `toStructure` 30 atoms, copy 2 whole; `toLoci` 20. Before #1994 the
  else-branch already used `ue` (wrong atoms when the arrays differ). Bundles are mostly made on roots (copies share
  elements), so rare. Fix: per-unit indices when `unit.elements !== ue`; push whole only for equal elements.
- **B12**: `findIntoBuilderWithRadius` finds units with `unitLookup.find(x, y, z, radius)` but accepts atoms up to
  `pivotR + r + radius`; units between `radius` and `radius + pivotR` are skipped. `surr2.ts` on flat-32, radius 0.5,
  vdW: 1,180 of 186,496 query atoms miss 1,212 neighbours (small units: ions, waters, ligands). Only MolScript
  `include-surroundings` with `atom-radius`. Fix: `unitLookup.find(x, y, z, queryRadius)`.
- **B11**: `sourceIndexMap` is keyed by operator name only; in a multi-model structure (`ofTrajectory`) every model's
  units are `1_555`. Failure (`expr.ts`, 1grm 5 models): `{m0:u0:[0,1,2], m1:u2:[5,6,7]}` round-trips to
  `{m0:u0:[0,1,2,5,6,7]}`. No in-repo caller (public API). Fix: key by model id + operator name, and group `opData` by
  model too.
- **B2**: coarse units index `atomicHierarchy.chainAtomSegments` with coarse element indices; 8zzh gives `[]` instead of
  `{0, 1}` (v1 `ent.ts`); consumer mmCIF export `_entity` / `_entity_poly*`. Fix: the coarse `chainElementSegments`.
- **B6**: `pick` leaks a `ctx.element` entry and pops an unpushed `currentStructure`; `if (i % N)` runs the timeout
  check on 99 % of iterations (free while `timeoutMs` is 0). Fix: push/pop the structure, `i % N === 0`.
- **B7**: a bundle for another structure is applied by unit id and UnitIndex. Audit `selection.ts` bundle callers,
  then throw or return empty.
- **B13** (plausible): `areHierarchiesEqual` and `areEquivalent` compare `hashCode` (unit ids + `SortedArray.hashCode`
  = size, first, last, middle) and never the element sets. #2022's `remappable` guard
  (structure-component.ts:199–211) remaps a static component when that holds; `[1,3,5,9]` and `[1,4,5,9]` collide.
  A static component below a source whose elements change between frames (a nested component above it) could keep
  stale elements. The weakness predates #2022 (`areEquivalent` was the guard). Fix: `areUnitIdsAndIndicesEqual` in the
  guard (short-circuits on identical arrays).
- **B3** (plausible): two different (set, ranges) with equal `hash2(fnv, fnv)` merge; P ≈ n²/2³³. Fix: equality check
  on hash hits (also loci.ts:805).
- **B8, B4, B5**: latent (identity operators by default; `Partitioned` set nowhere; no callers). Fixes as in v1.

## Checked and dropped

- **P8** (v1): unreachable, unchanged.
- `Loci.toStructure` of a root's loci has no `Structure.parent` (`parent: loci.structure.parent`), so its inter-unit
  bonds would not restrict the root's. Not a finding: unit.md §5 records it, and the callers (selection manager,
  Prism superposition/selection, interactions extension) do not read its bonds.
- R3 (v1 profile, 3k7f `lookup+boundary`): does not reproduce warm.

## Bench scripts

`/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/1f4e2fbb-da6c-4948-aacb-385a4dab465c/scratchpad/perf-audit-2/model-structure/`.
Bundle with `node build.mjs X.ts X.mjs` (esbuild, conditions `molstar-src`, `import`, `default`), run with
`node --expose-gc --max-old-space-size=8000 X.mjs <fixture>`; `load.ts` holds the fixtures. `verify/` builds with
`node ../build.mjs`.

| Script | Measures |
| --- | --- |
| `pairs.ts` | unit pairs, symmetry classes, water/ion pairs, bonds-inter at 3.1 Å (P13, P16, P17) |
| `symbonds.ts` | class-cached inter-unit bonds vs current (P13) |
| `xform.ts` | inter-unit bonds of a transformed structure and its component (P14) |
| `bnd.ts`, `bnd2.ts` | structure boundary with in-place extrema transform, GC forced / interleaved (P15) |
| `conn.ts` | fixed-point and layered `includeConnected` (P18) |
| `iub.ts` | polymer component query + union, subset of inter-unit bonds, `InterUnitGraph` constructor (P19, P20) |
| `verify/pairs-half.ts` | per-atom 3.1 Å queries vs a half-stencil pair walk vs `computeIntraUnitBonds` (P21) |
| `remapx.ts` | remapped component of a transformed structure (B9) |
| `verify/b14.ts` | `within` paths vs centre and surface distance (B14) |
| `verify/b10.ts` | `Bundle.toStructure` on a child with unequal copies (B10) |
| `surr.ts`, `surr2.ts` | surroundings reverse check; atom-radius unit pre-filter (P10, B12) |
| `expr.ts` | `Loci.toExpression` round trip on a multi-model structure (B11) |
| `mates.ts`, `r3.ts` | symmetry mates; 3k7f `lookup+boundary` (R3) |
