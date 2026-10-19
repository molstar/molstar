# Slice 05: model-structure

`packages/model/src/model/structure` (`Structure`, `Unit`, element/Loci/Bundle, bonds, boundary, symmetry, query).
Measured 2026-10-19 on `v6-opts-261019` (`90d717f49`), Node 25.6, Apple M4. Fixtures are those of
[00-profile.md](00-profile.md): `7qpd` (11.5k atoms, 7 units), `big-bcif` (186k atoms, 416 units, nothing shared),
`1tup-2x2x2` (186k atoms, 288 units, 9 symmetry groups), `8zzh` (coarse only). This is a read-only audit: the fixes
below are recommendations and none of them is applied. Every finding was checked again by a verifier, and the
numbers below are the verifier's reruns.

## Hot paths

Building the structure is cheap. `Structure.ofModel` and the symmetry and assembly builders take 9 ms for 186k atoms,
because units share element sets and props across operators. The cost is in the lazy caches and in the interactive
query and Loci layer.

| Path                                                     | Cost (fixture)                         | Where the time goes                                                                                                 |
| -------------------------------------------------------- | -------------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| `Unit.boundary` + `Structure.boundary` + `lookup3d`      | 74 ms (big-bcif), 16 ms (1tup-2x2x2)   | EPOS 98-direction dot products: unit spheres 33–52 ms, structure sphere 11–26 ms, grids 6 ms                        |
| `Unit.bonds` (intra)                                     | 242–283 ms (big-bcif)                  | grid query at 4 Å for every atom (67 % total), string-keyed element lookups per neighbour                            |
| `Structure.interUnitBonds`                               | 54–59 ms, 90 MB allocated (1tup-2x2x2) | 1.65M A-atom iterations over 3,956 unit pairs, of which 89 % fail the B-sphere test after a full transform          |
| `structureUnion` (Sequence or singletons → `Structure`)  | quadratic in unit size                 | `SortedArray.union` copy for each input structure                                                                   |
| Child structures (components, selections, focus, bundles) | full recompute                         | intra-unit bonds, inter-unit bonds, boundary and lookup are not derived from the parent; bundle children lose symmetry sharing |

Already good, no action needed: schema selection (element ranges with scratch buffers, identity shortcut);
`StructureSubsetBuilder` prop sharing through `UnitEquivalenceBuilder`; the `Loci.size` and whole-unit boundary
shortcuts; geometry sharing through symmetry groups.

## Performance findings

Confirmed and plausible findings, ranked.

| Id  | Title                                                         | File                                                                                                         | Kind        | Impact | Estimate                                                                         | Effort | Verdict   | Frozen |
| --- | ------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------ | ----------- | ------ | -------------------------------------------------------------------------------- | ------ | --------- | ------ |
| P1  | `structureUnion` merges one structure at a time               | [structure-set.ts:24](../../packages/model/src/model/structure/query/utils/structure-set.ts#L24)             | time+memory | high   | 100k singletons in one unit: 8.7 s → 18 ms; residue union 36 → 2.7 ms             | S      | confirmed | no     |
| P6  | Child units recompute intra bonds instead of subsetting       | [unit.ts:305](../../packages/model/src/model/structure/structure/unit.ts#L305)                               | time        | medium | half-unit children, big-bcif: 90 → 1.6 ms (55×)                                   | M      | confirmed | no     |
| P2  | 98-direction EPOS for structure, loci and unit boundaries     | [util/boundary.ts:15](../../packages/model/src/model/structure/structure/util/boundary.ts#L15)               | time        | medium | structure sphere 11–14 → 3–4 ms; unit spheres 33 → 9 ms; lookup+boundary 74 → ~30 ms | S      | confirmed | no     |
| P3  | `Bundle.toStructure` / `Loci.toStructure` lose symmetry sharing | [bundle.ts:236](../../packages/model/src/model/structure/structure/element/bundle.ts#L236)                   | time+memory | medium | 1tup-2x2x2 half selection: boundary 42 → 10 ms, per-copy lookups/bonds avoided    | S      | confirmed | no     |
| P4  | `findPairBonds` scans every atom of A for each unit pair      | [inter-compute.ts:105](../../packages/model/src/model/structure/structure/unit/bonds/inter-compute.ts#L105)  | time+memory | medium | inter bonds 55 → ~25–30 ms; 90 → <10 MB allocated                                | M      | confirmed | no     |
| P5  | Child structures recompute inter-unit bonds                   | [structure.ts:240](../../packages/model/src/model/structure/structure/structure.ts#L240)                     | time        | medium | ~50 ms → <1 ms per large child when the root has bonds                           | M      | confirmed | no     |
| P7  | Intra `findBonds` queries 4 Å and uses string lookups         | [intra-compute.ts:241](../../packages/model/src/model/structure/structure/unit/bonds/intra-compute.ts#L241)  | time        | medium | bonds-intra −15–25 % (grid part 150 → 115 ms at 3.1 Å)                           | M      | plausible | no     |
| P10 | One `Structure` per atom in singleton selections; surroundings | [selection.ts:164](../../packages/model/src/model/structure/query/selection.ts#L164)                         | time        | low    | forEach 19 ms per 186k atoms; surroundings up to 2× faster                       | M      | plausible | no     |
| P9  | Carbohydrate links: all-pairs loop with a BFS per pair        | [carbohydrates/compute.ts:264](../../packages/model/src/model/structure/structure/carbohydrates/compute.ts#L264) | time        | low    | O(E²) → O(E·k); matters only for polysaccharide-heavy entries                    | M      | plausible | no     |
| P11 | `ElementSetIntraBondCache` never evicts                       | [element-set-intra-bond-cache.ts:13](../../packages/model/src/model/structure/structure/unit/bonds/element-set-intra-bond-cache.ts#L13) | memory      | low    | bounds growth in long sessions on index-pair models                              | S      | confirmed | no     |
| P12 | `IntraUnitBondMapping`: 3 × Uint32 per directed bond instance | [mapping.ts:74](../../packages/model/src/model/structure/structure/mapping.ts#L74)                           | memory      | low    | −12 B per directed bond instance (4 MB big-bcif, ~24 MB for 60-mers)              | M      | confirmed | no     |

### P1 `structureUnion` is quadratic per unit

- **Mechanism.** For every input structure and unit it calls `SortedArray.union(accumulated, u.elements)`, which
  allocates and copies the whole accumulated set. The `fullUnits` shortcut fires only once a unit is complete.
- **Callers.** `StructureSelection.unionStructure`, so every Sequence selection that becomes a structure: components
  from expressions, `toLociWithSourceUnits`, the targets of `within`/`intersectBy`/`exceptBy`/`isConnectedTo`. Also
  `LinearBuilder`/`UniqueStructuresBuilder.getSelection` when all results are singletons: `pick`, `within`,
  `withSameAtomProperties`, and the PyMOL/VMD `within`/`around` transpilations.
- **Evidence** (`uni.ts`, `uni2.ts`, `ms.ts`):

  | Case                                   | Current | concat + sort + dedupe per unit |
  | -------------------------------------- | ------: | ------------------------------: |
  | residues() union, big-bcif (32k)       |   36 ms |                          2.7 ms |
  | 186k singletons, 416 units (max 1535)  |  328 ms |                               – |
  | 10k / 25k / 50k / 100k singletons, one 186k-atom unit | 104 / 524 / 2,063 / 8,735 ms | 1.9 / 3.5 / 8.3 / 17.9 ms |
  | `within(all, all, 1 Å)`, big-bcif      |  793 ms |                   (mostly the union) |

  Single-unit systems such as MD/PDB/GRO files without chain IDs are where this hurts: seconds per query.
- **Fix.** Collect the element arrays per unit id and merge each unit once: concat, `sort`, `SortedArray.deduplicate`,
  or a mark array over the source unit's elements. Keep the `fullUnits` and single-structure shortcuts. The output
  was verified identical. `getSelection` for all-singletons can feed elements straight into one subset builder.
- **Risk.** The output must stay a sorted, deduplicated set per unit in the same unit order. The query tests cover
  this.

### P6 Child units recompute intra-unit bonds

- **Mechanism.** `Atomic.getChild` creates empty `AtomicProperties`, so the `bonds` getter reruns
  `computeIntraUnitBonds` (grid plus a 4 Å query per atom). Distance-based bonds of a child are cached nowhere.
- **Restriction is exact.** Every decision in `findBonds` and `findIndexPairBonds` is per pair: struct_conn, the
  component dictionary, the hydrogen fallback, alt-locs, the thresholds, coarse trace. `canRemap` can be inherited:
  a subset of watery, dictionary-based or sequenced bonds keeps the property.
- **Applies when** the parent unit's bonds exist. Full-unit components share the root's `Unit` object, so focus,
  surroundings and selection children of those units usually find the parent bonds computed.
- **Evidence** (`sub.ts`, big-bcif, every second atom of 320 units): recompute 90 ms, subset 1.6 ms, identical edge
  counts (19,584).
- **Fix.**
  1. In `getChild`, pass a reference to the parent props (or a lazy thunk).
  2. In the `bonds` getter, if `parent.props.bonds` exists, map the parent UnitIndex to the child UnitIndex by walking
     both sorted arrays.
  3. Keep the edges whose two ends are present, carrying flags, order and key; inherit `canRemap`.
  4. Keep `ElementSetIntraBondCache` for the case where the parent bonds are absent.
- **Risk.** Do not retain the parent props forever: hold a reference to `bonds` only, or clear it after use. Add a
  test that subset bonds equal recomputed bonds, with struct_conn and alt-loc fixtures.

### P2 Boundary: 98 directions where 14 plus refine is enough

- **Mechanism.** Three call sites use the 98-direction helper:
  - `computeStructureBoundary` whenever there are 500 units or fewer: `includeSphere` over 196 extrema per unit × 98
    directions;
  - `Unit.boundary` for units of up to 10k atoms (`core/math/geometry/boundary.ts:22`);
  - `Loci.getBoundingSphere` ([loci.ts:631](../../packages/model/src/model/structure/structure/element/loci.ts#L631)),
    per element.

  `refineExtrema` projects the 14-direction extrema onto the 98 directions anyway.
- **Verifier check: enclosure.** The structure radius pass runs over the unit spheres' extrema, not over the atoms,
  so the extrema quality decides whether atoms end up outside the sphere (`bnd2.ts`):

  | Unit helper / structure helper | big-bcif time | r      | Worst atom outside                 |
  | ------------------------------ | ------------: | -----: | ---------------------------------- |
  | 98 / 98 (current)              | 11.4 ms       | 202.05 | 0.02 Å (7qpd 0.01)                 |
  | 98 / 14                        | 3.7 ms        | 197.76 | 0.02 Å (7qpd 0.01): same as now    |
  | 14 / 98                        | 11.1 ms       | 202.00 | 0.12 Å (7qpd 0.88 Å, 4 atoms)      |
  | 14 / 14                        | 5.7 ms        | 197.66 | 0.13 Å (7qpd 0.79 Å)               |

  Unit spheres stay exact with 14 directions, because their radius pass visits every atom. Unit boundaries:
  33–52 ms → 9–19 ms (big-bcif). A radius pass on the unit spheres' centre plus radius encloses all atoms but is
  5–17 % larger, so it is not worth it.
- **Fix.**
  - Safe: always use '14' in `computeStructureBoundary` and `Loci.getBoundingSphere`. The loci radius pass visits
    every element, so it stays exact.
  - Conditional: lower the `Unit.boundary` threshold too, and accept structure spheres that miss atoms by up to about
    0.9 Å, or add a small pad.
  - Optional: deduplicate `Sphere3D` extrema (at most 28 distinct points with 14 directions).
- **Estimate.** lookup+boundary on big-bcif 74 ms → ~55 ms (safe part only) or ~30 ms (both parts). The auditor's
  "high" is downgraded to medium: this costs tens of ms once per structure.
- **Risk.** Camera fit and visual-test snapshots may shift slightly.

### P3 Bundle/Loci → Structure lose symmetry sharing

- **Mechanism.** `Bundle.toStructure` calls `unit.getChild(indices)` for every unit of a grouped-unit list, and each
  call creates fresh props. `Loci.toStructure` ([loci.ts:175](../../packages/model/src/model/structure/structure/element/loci.ts#L175))
  does the same. `StructureSubsetBuilder` avoids this with `UnitEquivalenceBuilder` + `pivot.applyOperator`.
- **Callers.** Bundle components (`structure-component.ts:161/239`); `StructureSelectionFromBundle`
  (`selection.ts:414/434`), which the focus representation uses; the MolScript `bundle` generator (`internal.ts:127`);
  Prism selection, component and superposition; plugin UI superposition.
- **Evidence** (`ms.ts`, every second atom of 1tup-2x2x2): 195 distinct props against 9 for the subset builder;
  structure boundary 42 vs 10 ms; all-unit intra bonds 80 vs 4.7 ms.
- **Verifier correction.** Symmetry groups still form (9), so instanced visuals build geometry only for
  `group.units[0]`. The 16× bond figure applies only to code that touches every unit's props: structure boundary,
  per-unit `lookup3d` in inter-unit bonds and contacts, picking. Downgraded to medium.
- **Fix.** Per grouped list: `first = _units[0].getChild(indices)`; if `first === _units[0]`, push the original units;
  otherwise push `first.applyOperator(u.id, u.conformation.operator, true)` for the rest. `Loci.toStructure`: go
  through `structure.subsetBuilder` or a `UnitEquivalenceBuilder`.
- **Risk.** `applyOperator` keeps the pivot's `chainGroupId`, which is the subset builder's behaviour today. Unit ids
  and operators stay unchanged.

### P4 `findPairBonds` full scan per unit pair

- **Mechanism.** Every atom of A is transformed into B's frame and then tested against B's sphere (:109). Only 11 % of
  them pass (185,450 of 1,655,418 on 1tup-2x2x2). The function deoptimizes, and heap sampling attributes 83 of 90 MB
  per call to `findPairBonds` itself (~52 B per iteration) for 3,520 edges.
- **Evidence** (`pairs.ts`): `computeInterUnitBonds` 54–59 ms, 90–92 MB. Fetching the same 185,450 atoms with
  `A.lookup3d.find(inv(opA)·opB·cB, rB + maxRadius)` takes 5.9 ms, against a 35.5 ms scan.
- **Fix.**
  1. Compute B's centre in A's frame and query `unitA.lookup3d` into a reused result.
  2. Iterate only the returned indices. The struct_conn and index-pair branches stay per candidate atom; they already
     sit after the sphere test, so the candidate set is identical.
  3. Use inlined matrix arithmetic held in locals.

  `interactions.ts:464` can reuse the same pattern.
- **Risk.** The edge insertion order changes. Compare edge sets on assemblies.

### P5 Child structures recompute inter-unit bonds

- **Mechanism.** The `interUnitBonds` getter reuses the parent's bonds only when the parent has 0 edges.
  `StructureSubsetBuilder._getStructure` returns a new `Structure` even when every unit was reused.
- **Evidence** (`ms.ts`): polymer child after the root's bonds were computed: big-bcif 48 vs 48 ms; 1tup-2x2x2 56 ms
  root vs 62 ms child.
- **Applies when** the root has bonds. Valence model, interactions, `includeConnected` (focus) and `isConnectedTo`
  all compute them on the root.
- **Fix.**
  - (a) If the parent bonds exist and `dynamicBonds`, `interBondsValidUnit` and `interBondsValidUnitPair` are the
    parent's, filter the parent's edges: keep pairs whose units are both present, remap partial units through
    `SortedArray.indexOf`, drop missing atoms. All criteria are pairwise; the water/ion policy is per unit.
  - (b) Return `this.parent` when every unit was reused and the counts match, as `substructureFromRows` does.
- **Risk.** (b) changes `child !== parent`. `structure-component` already handles `component === a`; check the
  hierarchy manager.

### P7 Intra `findBonds` radius and per-neighbour lookups (plausible)

- **Radius.** The largest distance threshold reachable from the tables is ~2.95 Å: element thresholds ≤ 2.88 give
  (2.88 + 2.88)/1.95, and pair thresholds are ≤ 2.87. The 4 Å query exists only so that component-dictionary bonds
  are found at any distance.
- **Evidence** (`q.ts`): unit queries alone take 150–156 ms at 4 Å (12.3 neighbours/atom), 115–118 ms at 3.1 Å,
  96 ms at 2.5 Å.
- **Fix.**
  - Use a 3.1 Å radius (keep 4 Å for coarse models).
  - Resolve dictionary pairs through the residue atom range instead of the grid.
  - Use `derived.atom.atomicNumber` instead of `getElementIdx(type_symbol.value())`.
  - Keep a per-residue struct_conn flag before the cantor `Set` lookup.
- **Risk.** Dictionary bonds longer than 3.1 Å (bad geometry) change unless the dictionary path is made independent
  of the radius. Not measured end to end.

### Low

- **P10** (plausible). `StructureSelection.forEach` on Singletons builds a `Unit` and a `Structure` per atom: 19 ms
  per 186k atoms, then a P1 union. `within(all, 8.4k-atom target, 5 Å)` takes 63 ms. A reverse check for
  `includeSurroundings` halves the time on compact queries (42 → 20.5 ms, 1tup-2x2x2) and ties on big-bcif (23 vs
  26 ms), so choose it adaptively. Add singleton fast paths into one subset builder.
- **P9** (plausible, not measured). An E²/2 loop over all ring elements, with `areVertexSetsConnected(…, 3)` for
  every same-unit pair in different residues. Cross-unit pairs are rejected cheaply. Only polysaccharide or
  glycan-heavy units are affected. Fix: a depth-3 BFS per ring through an atom → ring-element index.
- **P11** (confirmed). `set()` only pushes, and entries live as long as the model. They are cacheable for every
  index-pair model, including mmCIF with `molstar_bond_site`. Each focus or selection child adds its element set
  and graph. Fix: an LRU cap, or skip caching small sets. P6 makes most of this cache unnecessary.
- **P12** (confirmed). `unitIndex`, `unitEdgeIndex` and `unitGroupIndex` are derivable from `unitGroupOffset` and the
  per-group `edgeCount`. Consumers: `bond-intra-unit-line.ts:237`, `bond-intra-unit-cylinder.ts:442`,
  `util/bond.ts:227` (repr, not `gl/`). The structure-level bond geometry is several times larger than this mapping.

## Bugs

| Id  | Title                                                                       | File                                                                                                       | Severity | Verdict   |
| --- | --------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------- | -------- | --------- |
| B1  | `remapModel` builds a private root per component                            | [structure.ts:410](../../packages/model/src/model/structure/structure/structure.ts#L410)                   | medium   | confirmed |
| B2  | `entityIndices` uses atomic chain segments for coarse units                 | [structure.ts:527](../../packages/model/src/model/structure/structure/structure.ts#L527)                   | low      | confirmed |
| B4  | `handleUnitChainsPartitioned`: `const i = 0` infinite loop                  | [stats.ts:289](../../packages/model/src/model/structure/structure/element/stats.ts#L289)                   | low      | confirmed |
| B6  | `pick`: push element / pop structure; `i % N` timeout checks inverted        | [filters.ts:24](../../packages/model/src/model/structure/query/queries/filters.ts#L24)                     | low      | confirmed |
| B7  | `Bundle.toLoci` / `toStructure`: `new Error` never thrown                   | [bundle.ts:155](../../packages/model/src/model/structure/structure/element/bundle.ts#L155)                 | low      | confirmed |
| B3  | Bundle/Loci-expression grouping on a 32-bit hash without an equality check  | [bundle.ts:118](../../packages/model/src/model/structure/structure/element/bundle.ts#L118)                 | low      | plausible |
| B5  | `getStructureDeduplicate` deduplicates `this.currentUnit`                   | [subset-builder.ts:80](../../packages/model/src/model/structure/structure/util/subset-builder.ts#L80)      | low      | confirmed |
| B8  | `findMatesRadius` checks a world-space centre in a unit's local grid        | [symmetry.ts:332](../../packages/model/src/model/structure/structure/symmetry.ts#L332)                     | low      | confirmed |

- **B1.** Remapping a component also remaps its parent, `parent?.remapModel(m)`, so every component gets its own copy
  of the root.
  - *Trigger.* A trajectory frame change: `hierarchy.ts:154` remaps the root and `structure-component.ts:199` remaps
    each static component.
  - *Effect.* `comp.root !== root`, and two remaps of the same component give different roots (`remap.ts`). Root-type
    custom props (secondary structure/DSSP, ASA; `custom-structure-property.ts:72` uses `data.root`) are recomputed
    on each private root, for K components on every frame.
  - *Fix.* `remapModel(m, parent?)` and pass `a.root` (or `a`) from `updateStructureComponent`.
- **B2.** On 8zzh (0 atoms, 4 coarse units), `entityIndices` is `[]` but should be `{0, 1}` (`ent.ts`). In mixed
  models the segments are arbitrary. The only consumers are the mmCIF export `_entity` / `_entity_poly*` categories.
  Fix: use the coarse `chainElementSegments`, as `Stats` does.
- **B4.** The loop never advances `i`. A loci with two or more `Partitioned` units in one chain-operator group, where
  the first element is not in `firstCI`, hangs the main thread (labels, hover). `Unit.Trait.Partitioned` is set
  nowhere in the repo, so this is latent.
- **B6.** `pick` leaks one `ctx.element` stack entry per call and pops a `currentStructure` it never pushed, which
  clobbers an outer one. `getSelectionProperties` leaves `currentStructure` set as well. `if (i % 100)` runs on 99 %
  of the iterations (`filters.ts:28,90,114,128`, `combinators.ts:24`); this costs nothing while `timeoutMs === 0`.
- **B7.** A bundle whose hash does not match the structure is applied anyway, so UnitIndex sets land on another
  structure's units. That gives wrong atoms, or `undefined` when an index is out of range. Audit
  `selection.ts:414/434` before making it throw.
- **B3.** Two different (set, ranges) pairs with the same `hash2(fnv, fnv)` merge, and the second gets the first one's
  atoms. P ≈ n²/2³³ for n distinct per-unit selections, about 1 % at 10k. Fix: `Map<number, Entry[]>` plus
  `SortedArray.areEqual`. The same applies to `loci.ts:780`.
- **B5.** Every unit gets the last committed unit's elements. There are no callers. Fix or remove.
- **B8.** `closeUnit.lookup3d.check` gets a world-space centre and ignores the unit's operator. This is wrong only for
  input units with non-identity operators. The default input is `Structure.ofModel`, so it is latent.

## Checked and dropped

- **P8** (index-pair bonds with operators: per-unit scan of `bySameOperator`, `atomsToUnits` Set per atom). The
  complexity claim is right, but the path is unreachable in the repo. `hasOperators` needs `operatorA`/`operatorB` in
  the pair data, and no format or extension sets them: `molstar_bond_site`, topology and mol2 pass `indexA`/`indexB`
  only. `canGetFromIndexPairBonds` and `canAddFromIndexPairBonds` are therefore always false. Revisit if a producer
  of operator-mapped pairs appears.
- Already measured as a non-win: unit grids over compacted per-unit Float32 coordinates with Interval indices were
  40 % slower (polymorphic query).

## Bench scripts

`/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/ac6cef1a-5954-43a5-ae83-c22fef35f065/scratchpad/perf-audit/model-structure/`.
Bundle with `node ../build.mjs X.ts X.mjs` and run with `node --max-old-space-size=8000 X.mjs`. Fixtures are in
`../fixtures/`.

| Script     | Measures                                                                                                     |
| ---------- | ------------------------------------------------------------------------------------------------------------ |
| `ms.ts`    | residue union, surroundings / reverse check, forEach, `within`, Bundle vs subset builder, child inter/intra bonds (`big` / `sym` / `7qpd`) |
| `uni.ts`   | singleton union and `within(all, all, 1)`; `merge` arg = 7qpd as one unit                                    |
| `uni2.ts`  | (verifier) singleton union inside one 186k-atom unit vs concat + sort                                         |
| `bnd.ts`   | unit / structure boundary, 98 vs 14 directions                                                               |
| `bnd2.ts`  | (verifier) structure-sphere enclosure of every atom for each helper combination                               |
| `pairs.ts` | unit pairs, in-sphere atoms, grid fetch, `computeInterUnitBonds` time / allocation, heap sample `hp/alloc.json` |
| `sub.ts`   | (verifier) child intra bonds: subset parent graph vs recompute                                                |
| `q.ts`     | per-atom unit grid queries at 4 / 3.1 / 2.5 Å                                                                 |
| `remap.ts` | `remapModel` root identity and property sharing (B1)                                                          |
| `ent.ts`   | `entityIndices` on 8zzh (B2)                                                                                  |
