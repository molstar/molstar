# Slice 06: model props and script

Scope: [`packages/model/src/props`](../../packages/model/src/props) (computed properties: ASA, interactions, valence
model, secondary structure / DSSP, volume streamlines, `CustomStructureProperty`) and
[`packages/model/src/script`](../../packages/model/src/script) (the MolScript compiler and runtime table). Query
generators (`model/structure/query`) and the theme callers in `graphics/src/props` sit in other slices. Fixtures and
method: [00-profile.md](00-profile.md). Nothing here touches Canvas3D or `gl`; no finding is in a frozen area.
Measured 2026-10-19, Node 25.6, Apple M4.

## Hot paths at 186k atoms

| Stage | 1tup-2x2x2 | big-bcif | Where the time goes |
| --- | ---: | ---: | --- |
| `asa` | 1.92 s / 8.3 GB alloc | 1.83 s / 7.5 GB | `computeRange` ~55 % (per-atom arrays, sort); `structure.lookup3d.find` ~40 % (750 ms on 1tup-2x2x2) |
| `interactions` | 0.49 s | 0.89 s | inter-unit loop: feature scan + transform 150 / 72 ms, grid queries ~200 ms; `getBoundary` of feature grids ~150 ms on big-bcif; intra contacts 34 % of big-bcif (416 groups), 4 % of 1tup-2x2x2 (9 groups) |
| `valence-model` | 48 ms / 92 MB | 48 ms / 93 MB | per unit, not per group; a 4-tuple per atom |
| DSSP (`type=dssp`) | 36–53 ms | | fine |
| secondary structure from the model | 2 ms | | fine |
| MolScript | | | compile 2 µs; compiled atom-groups query 2.4–5.7 ms vs 1.4–3.3 ms native (~1.7×, acceptable) |

Already efficient: features and intra contacts per symmetry group (`interactions.ts:383/391`), the DSSP trace lookup,
the per-pair bounding-sphere test, the valence cache per unit and props.

## Performance findings (confirmed and plausible, ranked)

| Id | Title | File | Kind | Impact | Estimate | Effort | Verdict | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| P1 | ASA `computeRange` allocates per atom (array-of-arrays neighbours, sort, destructuring) | [area.ts:64](../../packages/model/src/props/computed/accessible-surface-area/shrake-rupley/area.ts#L64) | time+memory | high | ASA 1.6–2.4×; alloc 8 GB → < 10 MB at 186k atoms | S | confirmed | no |
| N1 | Inter-unit loop transforms and sphere-tests every feature of A for every close pair | [contacts.ts:231](../../packages/model/src/props/computed/interactions/contacts.ts#L231) | time | high | with P6: inter loop 352 → 138 ms (1tup-2x2x2, −44 % of interactions), 270 → 132 ms (big-bcif) | S | confirmed | no |
| P7 | ASA neighbours through the full structure lookup (per-unit inverse transforms, H, water) | [area.ts:58](../../packages/model/src/props/computed/accessible-surface-area/shrake-rupley/area.ts#L58) | time | medium | lookup 750 → 345 ms on 1tup-2x2x2; 29 → 24 ms on 7qpd | M | confirmed | no |
| P2 | Interactions feature cache keyed by all props; never evicted | [interactions.ts:356](../../packages/model/src/props/computed/interactions/interactions.ts#L356) | time+memory | medium | +13 MB retained per param change on big-bcif (unbounded); ~15 % of each recompute | S | confirmed | no |
| P6 | Feature grids use the default 32 per cell, too coarse for 3.5–6 Å queries | [features.ts:116](../../packages/model/src/props/computed/interactions/features.ts#L116) | time | medium | grid queries 1.4–1.6× (−80 ms on 186k atoms) | S | confirmed | no |
| P5 | Feature and subset grids pay for the 98-direction `getBoundary`; subsets rebuilt per tester per run | [features.ts:143](../../packages/model/src/props/computed/interactions/features.ts#L143) | time | medium | −110–140 ms (13–15 %) of big-bcif interactions | S | confirmed | no |
| P4 | `validPair` runs before any type-compatibility check | [contacts.ts:244](../../packages/model/src/props/computed/interactions/contacts.ts#L244) | time | medium | 8–10 % of interactions (40–60 ms at 186k) | S | confirmed | no |
| P3 | Inter-unit contacts recomputed for every symmetry-equivalent unit pair | [interactions.ts:464](../../packages/model/src/props/computed/interactions/interactions.ts#L464) | time | medium | 1tup-2x2x2 has 441 pair classes for 4282 pairs; gain on what is left after N1/P6, up to ~60× on icosahedral assemblies | L | plausible | no |
| P8 | Valence model per unit (not per group), tuple and closure per atom | [valence-model.ts:371](../../packages/model/src/props/computed/chemistry/valence-model.ts#L371) | time+memory | low | 48 → ~20 ms, 92 → < 10 MB | M | confirmed | no |
| P9 | Intra refinement repeated per symmetry member; string keys | [refine.ts:58](../../packages/model/src/props/computed/interactions/refine.ts#L58) | time | low | refine 9–20 → ~3 ms on 1tup-2x2x2 | S | confirmed | no |

### P1: ASA `computeRange` allocation

- Mechanism: per atom `neighbors = []` of 5-element arrays (`:64`, marked TODO), comparator sort (`:88`),
  destructuring `[sX,sY,sZ]` / `[, sqRadius, nX, nY, nZ]` inside the 92-points × neighbours loop, `unitIndexMap.get`
  per lookup hit, `bElementIndex` computed twice.
- Evidence (`asa-bench`, rerun): 7qpd 124 → 51 ms, 988 → 1 MB allocated; 1tup-cell 203 → 126 ms, 712 → 5 MB; area
  difference 0.
- Fix: module-level growable `Float64Array`s for (sqRadius, nX, nY, nZ); sphere points flattened once in
  `initialize()`; drop the sort and test the last occluder first (Eisenhaber); cache the unit offset while `lUnits[iI]`
  is the same unit; `bI === aI` for the self test. Reference: `computeAreaFast` in `asa-bench.ts`.
- Risk: none seen (the sort only affects early exit). Add a test comparing areas with today's output on a fixture.

### N1: Inter-unit feature scan

- Mechanism: `addStructureContacts` loops over all `featuresA.count` features of unit A, transforms each into B's
  space and tests it against B's sphere. On 1tup-2x2x2 that is 4.25M transforms for 0.44M survivors (10 %).
- Evidence (`verify breakdown`): scan only 150 ms (1tup-2x2x2) / 72 ms (big-bcif). Prefilter through A's own feature
  grid — query `featuresA.lookup3d` at B's sphere centre mapped into A space (`inv(imageTransform)`) with radius
  `rB + maxDistance` — gives the same 438,794 / 416,010 candidates in 24 / 21 ms. Full loop (scan + find): default
  352 / 270 ms; prefilter 225 / 213 ms; prefilter + 4 per cell (P6) 138 / 132 ms.
- Fix: per pair, `Mat4.invert` once, `featuresA.lookup3d.find(c', rB + maxD)`, iterate its indices (keep the sphere
  test), then the existing per-feature code. Pass a module-level `Result` to the outer `find`
  (`GridLookup3D.find(x, y, z, r, result)`): members of one symmetry group share one `Features` object, so
  `featuresA.lookup3d === featuresB.lookup3d` for same-group pairs and the inner query would overwrite the outer result.
- Risk: same candidates; they arrive in grid-cell order, not index order, which only changes inter edge numbering
  (check that no consumer depends on it, or sort the few candidate indices).

### P7: ASA neighbour lookup

- Mechanism: `structure.lookup3d.find` per atom: unit-boundary query, inverse transform per close unit, a
  32-per-cell grid over all elements (H, water, non-polymer included and then skipped), radius fixed at 6.54 Å.
- Evidence (`asa-lookup`): one world-space `GridLookup3D` over the contributing atoms (Float32 xyz, 4 Å cells,
  `getFastBoundary`): 1tup-2x2x2 750 ms → 18 ms build + 327 ms queries (9.19M → 8.77M hits); 7qpd 29 → 24 ms. Gain is
  mostly the per-unit overhead in many-unit structures; models with explicit H gain more (not measured).
- Fix: build that grid in `computeArea`, map grid indices back to serial indices; query radius
  `probe + r1 + probe + maxAssignedVdw` (also fixes B6). Removes the `unitIndexMap` lookup.
- Risk: +12 B per contributing atom during the computation. Areas must equal today's (test).

### P2: Feature cache key

- Mechanism: `cacheKey = JSON.stringify(p)` covers every tester, contact and bridge param; `findUnitFeatures` stores
  `features-${cacheKey}` on `unit.transientCache`. Features depend only on the enabled feature providers (no provider
  reads params; checked all `Features.Provider(...)`), so any distance/angle change recomputes all features, and old
  entries stay for the unit's lifetime.
- Evidence (`probe featcache`, rerun on 7qpd): 10 changes of `hydrogen-bonds.distanceMax` → 10 `features-` entries on
  unit 0, +0.5 MB each; big-bcif +13 MB each (auditor). Recompute 840 ms on big-bcif of which features ~79 ms and their
  `getBoundary` ~66 ms.
- Fix: key from the selected feature providers' types (plus the topology identity, see B2); one slot per unit (fixed
  symbol key holding `{ key, features }`).
- Risk: none if providers stay param-free; document that in the spec.

### P6: Feature grid cell size

- Mechanism: `GridLookup3D(position, boundary)` with the default ~32 points per cell → ~8 Å cells at feature density;
  a 6 Å query scans ~27 cells.
- Evidence (`probe cells`, `verify breakdown`): 1tup-2x2x2 query loop 352 → 256–274 ms at 4 per cell (2 per cell and
  4 Å cells the same); with N1 223 → 138 ms. Build time unchanged (3–4 ms).
- Fix: pass `4` (or `Vec3.create(4, 4, 4)`) in `Features.create` and `createSubset`; leave the core default.
- Risk: slightly more memory for sparse subsets. Results unchanged.

### P5: Feature `getBoundary`

- Mechanism: under 10k points `getBoundary` uses the 98-direction helper (49 dots per point, plus a distance pass)
  only to get a box and a sphere for the grid and the inter prefilter. `_addUnitContacts` calls
  `features.subset(tester.requiredFeatures)` for each tester per group per computation (no memo), so 8 subset index
  arrays and grids are rebuilt even when features come from the cache.
- Evidence (`verify boundary`): the 416 full feature sets of big-bcif: `getBoundary` 60–65 ms vs `getFastBoundary`
  3 ms. Auditor profile: subset grids another ~90 ms.
- Fix: `getFastBoundary` for feature and subset grids; memoize subsets on the `Features` object by a type-set key.
- Risk: a looser sphere (from the box) lets more features through the N1/inter prefilter; they fail in the grid query.

### P4: Type check before `validPair`

- Mechanism: for each hit at the max tester distance (6 Å, cation-pi), `validPair` reads altLoc strings, compares
  residues and calls `connectedTo` before the tester loop rejects incompatible types.
- Evidence (`probe inter`, rerun): 2,225,450 hits, 124,038 (5.6 %) with a tester whose required features hold both
  types within its distance. Auditor profile: `validPair` 47 ms of 512 (1tup-2x2x2), 60 of 950 (big-bcif).
- Fix: before the loop a 14×14 `FeatureType × FeatureType → tester bitmask` and max distance² per mask; skip a hit
  when the mask is 0 or the distance exceeds it, then `validPair`, then only the testers in the mask, in the same order.
- Risk: none if every `getType` returns falsy for types outside its `requiredFeatures` (true for all eight; keep the
  order so the first accepting tester still wins).

### P3: Symmetry-equivalent unit pairs

- Mechanism: `addStructureContacts` runs for each of the 4282 close pairs of 1tup-2x2x2, which fall into 441
  (invariantA, invariantB, relative transform) classes.
- Not exact as proposed: `getType` is not pair-local (H-bond geometry walks `eachBondedAtom(structure, …)`, which
  includes inter-unit bonds that differ at block edges), and line of sight queries `structure.lookup3d`, which sees
  third units. A safe variant caches per class the `(i, j, type)` candidates only for units without inter-unit bonds and
  reruns `validPair` + line of sight per concrete pair.
- Estimate: after N1/P6 the inter loop on 1tup-2x2x2 is ~140 ms plus testers and line of sight; P3 removes most of
  the candidate search there. Worth it only for large assemblies; needs a spec decision.
- Risk: high (semantics); test equality against today's output on 1tup-2x2x2.

### Low

- **P8** valence model: `calcValenceModel` caches per unit; `applyOperator` copies get no cache, so 288 units are
  computed, not 9 groups. Per atom a `[charge, implH, totH, geom]` tuple, closures in `bondToElementCount` /
  `explicitValence`, 5–7 `interUnitBonds.getEdgeIndices` lookups (each miss allocates `[]`). Fix: write into the
  output arrays, inline without closures, skip inter-unit lookups for units without inter-unit edges. Per-group sharing
  needs copy-on-write near inter-unit bonds (block edges); do it only with a test on 1tup-2x2x2.
- **P9** refine: 288 `unitsContacts` keys share 9 objects; template-string keys and `[dist, edge]` tuples per
  hydrophobic edge; `elementKeyIndex` uses `${index}|${unitId}` strings and returns a new `[]` on misses. Measured 9–20
  ms on 1tup-2x2x2. Fix: `cantorPairing` numeric keys and a shared empty array. Refining each shared object once
  changes the result of the weak-H-bond refiner (see N2): decide N2 first.

## Bugs

| Id | Title | File | Severity | Verdict |
| --- | --- | --- | --- | --- |
| B1 | ASA `getValue` binary-searches an unsorted, duplicated array | [shrake-rupley.ts:133](../../packages/model/src/props/computed/accessible-surface-area/shrake-rupley.ts#L133) | high | confirmed |
| B2 | Interaction features cached on shared `Unit`s but computed from the local substructure's bonds | [interactions.ts:406](../../packages/model/src/props/computed/interactions/interactions.ts#L406) | medium | confirmed |
| B3 | ASA throws on any structure with coarse units | [radii.ts:39](../../packages/model/src/props/computed/accessible-surface-area/shrake-rupley/radii.ts#L39) | medium | confirmed |
| N2 | Weak-H-bond refiner filters shared intra contacts from one member's inter contacts | [refine.ts:189](../../packages/model/src/props/computed/interactions/refine.ts#L189) | low | confirmed (code) |
| B4 | `calcPlaneAngle` skips the second neighbour when it is `UnitIndex` 0 | [geometry.ts:147](../../packages/model/src/props/computed/chemistry/geometry.ts#L147) | low | confirmed |
| B6 | ASA lookup radius ignores radii appended to `VdWLookup` | [shrake-rupley.ts:89](../../packages/model/src/props/computed/accessible-surface-area/shrake-rupley.ts#L89) | low | confirmed |
| B7 | Water bridges ignore `validUnitPair` (cross-model bridges) | [water-bridges.ts:211](../../packages/model/src/props/computed/interactions/water-bridges.ts#L211) | low | confirmed |
| B5 | `getArray` drops the first object-form argument (ring fingerprint) | [table.ts:632](../../packages/model/src/script/runtime/query/table.ts#L632) | low | confirmed |
| B8 | DSSP H-bond graph stores energies in input order, not slot order | [backbone-hbonds.ts:153](../../packages/model/src/props/computed/secondary-structure/dssp/backbone-hbonds.ts#L153) | low (latent) | confirmed |
| B10 | Streamline `prev` is module state, not reset per line | [basic.ts:42](../../packages/model/src/props/volume/streamlines/basic.ts#L42) | low | confirmed |
| B9 | `CustomStructureProperty.attach` has no in-flight dedup | [custom-structure-property.ts:76](../../packages/model/src/props/common/custom-structure-property.ts#L76) | low | plausible |

- **B1** `SortedArray.indexOf(ofSortedArray(root.serialMapping.elementIndices), location.element)`: `elementIndices`
  concatenates every unit's elements, so it is unsorted with several models or operators and repeats each element per
  copy; `location.unit` is ignored. Failure (`probe asa`, rerun): 1tup-cell 9,910 of 23,312 elements get another
  element's value (max error 144 Å²); 1tup-2x2x2 86,312 of 186,496; 7qpd 0. Affects the ASA color theme
  (`graphics/src/props/computed/themes/accessible-surface-area.ts:50`) and `is-buried` / `is-accessible` on every
  assembly or symmetry structure. Fix: map to the root unit (`root.unitMap.get(location.unit.id)`), then
  `root.serialMapping.getSerialIndex(rootUnit, location.element)`; better, the theme caches the unit offset. Also
  ~2.8× faster theme eval (12.4 → 4.4 ms). Test on 1tup with a symmetry range and on a component substructure.
- **B2** `findUnitFeatures(structure, unit, …)` runs providers against the local structure (`eachBondedAtom`,
  `bondCount` see only `structure.interUnitBonds` of the component), stores the result in `unit.transientCache`, and
  `Unit.getChild` returns the same `Unit` when a component holds the whole unit. Interactions is a `local` property, so
  this is the normal component flow. Failure (`probe featorder`, rerun on 7qpd): computing on single-unit components
  first, then on the root, gives units 1 and 3 one feature more (2888 vs 2887, 2505 vs 2504) than a root-only load; for
  a covalent ligand the hydrophobic flag and donor/acceptor counts of the bonded atom change. Fix: compute features with
  the root topology as `calcUnitValenceModel` does (`valence-model.ts:310-321`), or put the topology identity in the key.
  Related, by design: group members share `group.units[0]`'s features although their inter-unit bonds can differ.
- **B3** `assignRadiusForHeavyAtoms` calls `residue.key` and `label_atom_id` on coarse units ("Property only available
  for atomic models."); `isApplicable` is `() => true`. Failure: 8zzh (IHM coarse) and mixed atomic + coarse entries
  throw in `asa` and the ASA theme (00-profile). Fix: skip non-atomic units (`atomRadiusType` 0,
  `serialResidueIndex` -1, advance `m`); `isApplicable: s => s.units.some(Unit.isAtomic)`.
- **N2** intra contacts of a symmetry group are one shared object (`interactions.ts:392-395`), refined once per member.
  `weakHydrogenBondsRefiner.handleIntraContact` filters a weak H-bond when the acceptor has an H-bond, checking the
  member's own inter contacts, and writes the flag into the shared array. Failure: in a symmetry range, an acceptor
  that makes an inter-unit H-bond in one copy (a crystal contact) hides that weak H-bond in every copy, including edge
  copies where the contact does not exist. Fix: refine intra contacts per group from intra data only, or give members
  their own flag arrays when an inter-dependent refiner fires.
- **B4** `if (ni === 1 && unitX1 && indexX1)`: index 0 is falsy, so an acceptor/donor whose only heavy neighbour is
  the first atom of its unit (e.g. acetate O–C1) gets `undefined` and skips the out-of-plane check in H-bond geometry.
  Fix: `!== undefined`.
- **B6** `maxLookupRadius = 2*probe + 2*1.87` is fixed before `assignRadiusForHeavyAtoms` appends larger radii (Se 1.9
  in MSE; metals with `nonPolymer`). Failure: pairs with `d` in (6.54, 2*probe + r1 + r2) are missed, ASA slightly
  overestimated near Se and ions. Fix: radius after radii assignment, or per-atom radius with P7.
- **B7** `findWaterBridgeContacts` accepts hits from any unit; inter contacts apply `Structure.validUnitPair`.
  Failure: an ensemble shown as all models, or an MD frame set with waters: a water of model 1 bridges atoms of model 2.
  The `unitPairTest` option is ignored too. Fix: skip hits with `!validUnitPair(structure, unitW, hitUnit)` and require
  it between donor and acceptor.
- **B5** object-form branch starts at `i = 1` without pushing `keys[0]` (cf. `math.sub`). Failure: a
  `ring-fingerprint` expression with object args (JSON / named args) loses its first element. Fix: start at 0.
- **B8** `createGraph({ energies })` passes the input-order `number[]`; `_energies` (slot order) is filled and thrown
  away. No consumer reads `edgeProps.energies` today. Fix: `createGraph({ energies: _energies })`.
- **B10** `prev` is never reset, so a new line's first point is dropped when it lies within tolerance of the previous
  line's last point; a `Vec3` is allocated per step even when not pushed. Fix: reset `prev` to NaN per direction call,
  allocate after the distance test.
- **B9** `attach` returns early only when a value exists; two concurrent attaches both run `obtain` (≈1 s interactions
  / ≈2 s ASA on big-bcif) and the last `set` wins, possibly with the other caller's props. Plausible only if the plugin
  or Prism ensure custom properties concurrently (state transforms are normally sequential). Fix: a pending-promise
  map per structure and descriptor, cleared on rejection.

## Checked and dropped

- None of the auditor's findings were refuted. Corrections: P2 impact high → medium (13 MB per change, not per slider
  frame; sliders commit on release); P5 "49 dots in each of two passes" — the second pass is one distance per point;
  the auditor's inter "query loop" numbers include the feature scan of N1, which is a separate and larger cost than the
  grid queries' share suggested.
- MolScript compile/runtime: not a hotspot (closure layering ~1.7× native) — no finding.
- DSSP: 36–53 ms on 186k atoms with the trace lookup already efficient — no finding.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/ac6cef1a-5954-43a5-ae83-c22fef35f065/scratchpad/perf-audit/model-props-script/`,
bundled with `../build.mjs <in.ts> <out.mjs>` (esbuild, `molstar-src` condition).

| Script | Modes / what |
| --- | --- |
| `probe.ts` | `asa` (B1 correctness, getValue timing), `inter` (pair classes, hits, refine), `inter2`, `cells` (P6), `featcache` (P2, needs `--expose-gc`), `featorder` (B2), `ms` (MolScript vs native), `dssp`; `node probe.mjs <mode> <fixture>` |
| `asa-bench.ts` | P1: original vs `computeAreaFast`; `node asa-bench.mjs 7qpd 1tup-cell` |
| `v/verify.ts` | `prefilter` / `breakdown` (N1 + P6 combinations), `boundary` (P5), `refine` (P9) |
| `v/asa-lookup.ts` | P7: structure lookup vs world-space grid of contributing atoms |

Fixtures `7qpd`, `1tup-cell`, `1tup-2x2x2`, `big-bcif` as in [00-profile.md](00-profile.md); big ones need
`--max-old-space-size=8000`.
