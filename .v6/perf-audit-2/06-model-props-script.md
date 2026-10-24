# Slice 06 (v2): model props and script

Scope: [`packages/model/src/props`](../../packages/model/src/props) (ASA, interactions, valence model, DSSP,
streamlines, `CustomStructureProperty`) and [`packages/model/src/script`](../../packages/model/src/script) (MolScript
compiler and runtime table). Re-audit of [v1](../perf-audit/06-model-props-script.md) at `0944823d4`. Fixtures and
method: [00-profile.md](00-profile.md) (flat-N = N 1tup copies as separate entries, one group per unit; asm-N = the same
as an assembly, 9 groups). Single-threaded on a shared machine (load 2–4): ±10–15 %. No finding touches Canvas3D or
`gl`. MVS is out of scope (`.v6/mvs-optimizations/research/`).

Changed since v1: #1986 (ASA `computeRange` without allocation), #1989 (ASA `getValue` on assemblies), #1996 (inter-unit
candidates from A's grid), #2004 (H-bond hydrogen across units, `filterInter`, salt-bridge refiner direction), #2017
(custom-model-property asset owner). Per-symmetry ASA (13/P6) not started.

## Hot paths now

| Stage | flat-32 / asm-32 | flat-192 / asm-192 | Where the time goes |
| --- | ---: | ---: | --- |
| `interactions` | 0.74 / 0.35 s | 4.8 / 2.2 s | inter-unit loop 42 % / 86 % at 32, 50 % / 83 % at 192 (grid `query` 32 % self); intra 33 % on flat (one group per unit); `getBoundary` of feature + subset grids 19 % of flat-192 (932 ms); line of sight 11 %; `validPair` 7 %; features 7 % |
| `asa` | 0.89 / 0.95 s | 5.6 / 5.6 s | asm-192 profile (5.9 s): `computeRange` self 45 %, lookup `query` + `add` + `_find` 48 % |
| `valence-model` | ~0.05 s (v1; scratch runs include lazy bonds) | 0.31 / 0.32 s | per unit, not per group (1728 computations on asm-192); `getEdgeIndices` 23–25 %; 484 MB alloc |
| MolScript, DSSP, helix orientation, SIFTS | | | not hot; nothing new |

- Inter loop on flat-32 (scratch replay, features cached): 309–312 ms for 2,267,946 grid hits → 41,864 edges.
- Intra on flat-32: 245 ms, of which halogen 46 ms and metal 51–53 ms find 0 contacts.
- Inter contacts retain 508 B per edge (20.3 of 24.0 MB per computation on flat-32).
- Allocation flat-32 (sampled): interactions 539 MB (H-bond geometry 26 %, `getBoundary` 22 %); valence 93 MB; ASA
  14 MB (v1 7.5 GB).
- Already good: the N1 candidate prefilter (exact, ordered), the ASA kernel, B1's root mapping, the valence model's
  root topology.

## v1 → now

23 v1 items: **3 fixed**, 0 partly fixed, **20 still open**, 0 regressed, 0 dropped.

| v1 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| P1 | ASA `computeRange` allocates per atom | fixed | #1986 | flat `Float64Array` neighbours (`area.ts:55-72`), flattened sphere points, last-occluder test, cached unit offset; ASA alloc 14 MB sampled (v1 7.5 GB) |
| N1 | Inter-unit loop transforms every feature of A | fixed | #1996 | `getCandidates` (`contacts.ts:219`) queries A's grid around B's centre; own `Result`, sorted; 1 % self in profiles |
| B1 | ASA `getValue` binary-searches an unsorted array | fixed | #1989 | `shrake-rupley.ts:134-144` maps through `root.unitIndexMap`; test added |
| P7 | ASA neighbours through the structure lookup | still open | | `area.ts:98`; world grid re-measured −23 % (flat-32), −28 % (asm-32), areas identical |
| P2 | Feature cache keyed by all props, never evicted | still open | | `interactions.ts:356`, `:407` |
| P6 | Feature grids at 32 per cell | still open | | ~9 % of the inter loop after P4 (v1: 1.4–1.6×) |
| P5 | `getBoundary` for feature/subset grids; subsets not memoized | still open | | 19 % of flat-192 interactions |
| P4 | `validPair` before any type check | still open | | now the largest exact win (−48 % inter loop) |
| P3 | Inter contacts per symmetry-equivalent pair | still open | | asm-192: 33,814 pairs in 267 classes; v1's safe variant would cover nothing (1536 / 1728 units have inter-unit bonds) |
| P8 | Valence model per unit; tuple and closures per atom | still open | | `valence-model.ts:371-381`, `:334`, `:83` |
| P9 | Intra refinement per member; string keys | still open | | `refine.ts:59-84`, `:123-124`; `common.ts:95` |
| B2 | Features cached on shared Units, computed from the local structure | still open | | `interactions.ts:406-421` unchanged |
| B3 | ASA throws on coarse units | still open | | `radii.ts:39`; `isApplicable` `() => true` |
| N2 | Weak-H-bond refiner writes shared intra flags from one member | still open | | `refine.ts:190-193` |
| B4 | `calcPlaneAngle` skips neighbour at UnitIndex 0 | still open | | `geometry.ts:147` |
| B6 | ASA lookup radius ignores appended radii | still open | | `shrake-rupley.ts:90` |
| B7 | Water bridges ignore `validUnitPair` | still open | | `water-bridges.ts:211-227` |
| B5 | `getArray` drops the first object-form argument | still open | | `table.ts:632` |
| B8 | DSSP H-bond energies in input order | still open | | `backbone-hbonds.ts:153` |
| B10 | Streamline `prev` is module state | still open | | `basic.ts:42` |
| B9 | `CustomStructureProperty.attach` without in-flight dedup | still open | | `custom-structure-property.ts:76` |
| 13/P6 | ASA ignores symmetry | still open | | interface share re-measured: 85 % of contributing atoms on 1tup (18.5 % on 7qpd) |
| 12/S1 | Streamline trace: Map per line, Vec3 per step | still open | | `basic.ts:64`, `:90` |

## Open performance findings (ranked)

| Id | vs v1 | Title | File | Kind | Impact | Estimate | Effort | Verdict | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| P4 | still open; **fixed** #2025 | Inter queries at the 6 Å max distance for every feature type; `validPair` before any type check | [contacts.ts:283](../../packages/model/src/props/computed/interactions/contacts.ts#L283) | time | high | inter loop −48 % (flat-32 309 → 165 ms, asm-32 293 → 155 ms); ~−1.2 s flat-192, ~−0.9 s asm-192; exact, same edge order; **measured** with P10 (06 alone, A/B): whole `interactions` stage flat-32 707–732 → 409–544 ms, asm-192 2.29–2.65 → 1.10–1.54 s; edges and order equal (1tup flat-32 / asm-32, defaults and all providers) | S | confirmed | no |
| P10 | new; partly fixed #2025 | Feature types with no possible partner are computed, gridded and queried (halogen acceptors, metal partners) | [interactions.ts:373](../../packages/model/src/props/computed/interactions/interactions.ts#L373) | time+memory | high | with P4: inter loop −57 % (310 → 133 ms flat-32); intra −40 % (−100 ms flat-32, ~−13 % of flat-192); loop part done (`ContactTester.acceptsTypes`: partner-type skips inter and intra); the features are still computed and gridded (provider pruning open) | M | confirmed | no |
| P5 | still open; **fixed** #2045 | Feature and subset grids pay the 98-direction `getBoundary`; subsets rebuilt per computation | [features.ts:143](../../packages/model/src/props/computed/interactions/features.ts#L143) | time+memory | high | grids 180–189 → 50–56 ms on flat-32; ~−10–15 % of flat interactions (flat-192 up to ~0.7 s); −100 MB alloc flat-32 | S | confirmed | no |
| P3 | still open | Inter contacts recomputed for every symmetry-equivalent unit pair | [interactions.ts:464](../../packages/model/src/props/computed/interactions/interactions.ts#L464) | time | high | after P4+P10 the search is 43 % of the loop (52 of 122 ms, asm-32); per-class reuse removes most of it: ~−40 % of the remaining inter loop on assemblies (~−2 s at asm-864) | L | plausible | no |
| P7 | still open | ASA neighbours through the structure lookup, not one world grid | [area.ts:98](../../packages/model/src/props/computed/accessible-surface-area/shrake-rupley/area.ts#L98) | time | medium | asa −23 % (flat-32 888–944 → 687–703 ms) to −28 % (asm-32); ~−1.3–1.6 s at 1.1M atoms | M | confirmed | no |
| P13 | new | Inter contacts retained as JS objects, three Maps and string keys: ~510 B per edge | [common.ts:91](../../packages/model/src/props/computed/interactions/common.ts#L91) | memory | medium | retained ÷10 (~0.5 GB at flat-864); ~5–8 % of interactions time on asm | M | confirmed | no |
| P2 | still open; **fixed** #2044 | Feature cache keyed by all params, never evicted | [interactions.ts:407](../../packages/model/src/props/computed/interactions/interactions.ts#L407) | time+memory | medium | −7 % per recompute after a param change; stops +13 MB growth per change at 186k atoms | S | confirmed | no |
| P11 | new | Line of sight queries the full structure lookup per accepted contact | [contacts.ts:103](../../packages/model/src/props/computed/interactions/contacts.ts#L103) | time | low | ~−5 % of interactions (half of its 11 %) | M | plausible | no |
| P8 | still open | Valence model per unit; tuple, closures and inter-unit lookups per atom | [valence-model.ts:371](../../packages/model/src/props/computed/chemistry/valence-model.ts#L371) | time+memory | low | asm-192 0.32 → ~0.05 s with group sharing; alloc ~÷10 | M | confirmed | no |
| P12 | new; **fixed** #2047 | H-bond geometry allocates arrays, closures and objects per candidate | [hydrogen-bonds.ts:241](../../packages/model/src/props/computed/interactions/hydrogen-bonds.ts#L241) | memory | low | −140 MB alloc on flat-32 (26 % of the stage); a few % time | S | confirmed | no |
| 13/P6 | still open | ASA computed for every atom of every assembly copy | [area.ts:87](../../packages/model/src/props/computed/accessible-surface-area/shrake-rupley/area.ts#L87) | time | low | ≤ (1 − interface share): ~15 % on 1tup assemblies, more with large chains | L | plausible | no |
| P6 | still open | Feature grids at the default 32 points per cell | [features.ts:116](../../packages/model/src/props/computed/interactions/features.ts#L116) | time | low | ~−9 % of the inter loop after P4 | S | confirmed | no |
| P9 | still open | Intra refinement per symmetry member; string keys, tuples | [refine.ts:59](../../packages/model/src/props/computed/interactions/refine.ts#L59) | time | low | −3–5 % of interactions on assemblies | S | confirmed | no |
| 12/S1 | still open | Streamline trace: Map per line direction, Vec3 per step | [basic.ts:64](../../packages/model/src/props/volume/streamlines/basic.ts#L64) | time+memory | low | 10–20 % of a ~40 ms trace | S | confirmed | no |

### P4: per-type query radius and type-pair mask

- Mechanism: per candidate i of A, `lookup3d.find` on B at `maxDistance` of all testers (cation-pi 6.0 Å), then
  `validPair` (altLoc strings, residue compare, `connectedTo`) and every tester's `getType` per hit. Default testers:
  cation-pi 6.0, pi-stacking 5.5, halogen 4.0, H-bond 3.5 / 4.1 (S), metal 3.0 Å; most features belong to the short ones.
- Evidence (`inter.ts`, re-run with an edge-by-edge check in `v/inter-eq.ts`): flat-32 2,267,946 hits → 392,066,
  `validPair` 2,267,946 → 129,414, loop 309–312 → 163–171 ms; asm-32 293–317 → 153–168 ms; 7qpd 11,340 → 1,568 hits.
  Edges identical in order and type on 7qpd, flat-32, asm-32 (`unitA,unitB,indexA,indexB,type` sequence equal).
- Fix: once per computation `radiusOf[type]` (max `maxDistance` of testers requiring the type) and a 16×16
  `mask[typeA][typeB]` of testers (max distance² per mask). Per candidate: query B at `radiusOf[typesA[i]]` (skip at 0);
  per hit skip on mask 0 or distance ≥ mask max; then `validPair`; then only the masked testers, in order (first
  accepting wins).
- Risk: exact while every `getType` returns undefined outside its type pairs (true for all 8 today: each tests the
  types first). Put it in the `ContactTester` contract; keep an equality test (the replay does it).

### P10: features with no partner

- Mechanism: `featureProviders` keeps every provider whose types meet `requiredFeatures` (`:368-373`). Defaults
  include halogen bonds and metal coordination, so every O/N/S bonded to C/N/P/S becomes a `HalogenAcceptor` (with an
  `eachBondedAtom` walk) and every backbone O, side-chain O/N/S and non-standard O/N a `DativeBondPartner` and/or
  `IonicTypePartner`, whether or not the structure has a halogen donor or a metal. Testers declare only a type set, so
  partner-partner and acceptor-acceptor hits reach `validPair` in the inter loop and the intra loop.
- Evidence: flat-32 has 275,904 features, 62 % halogen acceptors or metal partners, with 0 halogen donors and 96 metals
  (all in their own units). Intra per tester (`intra2.ts`, re-run): halogen 46–47 ms and metal 51–53 ms, 0 contacts, of
  ~245 ms. Inter with exact type pairs plus a skip when B has no partner type (`v/inter-eq.ts`): 309–312 → 132–133 ms
  (flat-32), 293 → 122–125 ms (asm-32); 293,638 hits; edges identical in order.
- Fix: each tester declares its accepted pairs (H-bond Donor–Acceptor; weak WeakDonor–Acceptor; halogen
  HalDonor–HalAcceptor; metal TM–Dative, TM–TM, IonicMetal–IonicPartner; pi Ring–Ring; cation-pi Ring–Positive; ionic
  Negative–Positive; hydrophobic H–H) and P4's mask uses them. Per `Features` a present-types bitset: inter skips a
  candidate whose partner bitset misses B's; intra skips a tester in a unit lacking one side of every pair. Optionally
  run the anchor providers (halogen donors, metals) first and drop `HalogenAcceptorProvider` / `MetalBindingProvider`
  when no unit has the anchor (not measured); then the provider set must enter the feature cache key (P2).
- Risk: the bitset part is exact (measured). Pruning providers makes feature sets depend on the structure's content;
  consumers read features only through contacts and bridges, but a test of edges with and without pruning is needed.

### P5: feature grid boundary and subset memo

- Mechanism: `Features.create` (`:116`) and `createSubset` (`:143`) build grids with `getBoundary`: ≤ 10k points the
  98-direction helper (49 dots per point + a distance pass), 10k–250k the 14-direction one, only to get a box and a
  sphere. `_addUnitContacts` calls `features.subset(...)` per tester per group per computation (`contacts.ts:166`); no
  memo, so recomputes with cached features rebuild 5 subset index arrays and grids.
- Evidence: flat-192 profile `getBoundary` 932 ms = 19 % (subset grids 557, full grids 427 ms); `Vec3.dot` 386 ms and
  `computeSphereExtrema` 181 ms self. flat-32 (`intra.ts`, re-run): 416 full + 5 × 416 subset grids 180–189 ms with
  `getBoundary`, 49–56 ms with `getFastBoundary`, 32–34 ms for the index arrays alone. Alloc: `getBoundary` 119 of
  539 MB. Counter-effect: on box-derived spheres the inter loop is 191–233 vs 164–171 ms (looser B-sphere prefilter and
  `getCandidates` radius), eating part of the gain.
- Fix: `Box3D` in one pass plus a centroid/max-distance sphere (two passes, `CentroidHelper`), not
  `Sphere3D.fromBox3D`. Memoize `subset()` on the `Features` object by a sorted type-set key.
- Risk: none for results (same points per grid). The memo keeps one grid per tester type set per unit for the
  features' lifetime (cached per unit already, P2).

### P3: symmetry-equivalent unit pairs

- Mechanism: `addStructureContacts` per close pair; pairs fall into (invariantA, invariantB, `inv(opB)·opA`) classes,
  and the candidate geometry (i, j, distSq after the masks) is identical within a class.
- Evidence (`classes.ts`): asm-32 4,282 pairs / 265 classes; asm-192 33,814 / 267 (127×); 256 / 288 and 1536 / 1728
  units have inter-unit bonds (Zn and struct_conn), so a reuse limited to units without inter-unit bonds (v1) covers
  nothing. Search-only replay after P4+P10 (`v/inter-eq2.ts`, asm-32): 52 of 122 ms (43 %) is candidate search and
  grid queries; `validPair` + `getType` + line of sight the rest.
- Fix: per class, compute once the (i, j, distSq) list passing P4/P10; per concrete pair replay `validPair`, the masked
  testers and line of sight on it. Full result reuse only with a spec decision (H-bond geometry walks inter-unit bonds,
  line of sight sees third units).
- Risk: candidate reuse is exact (rigid relative operator, same features per group); needs an equality test on asm-32
  and asm-192 and a key tolerant to float noise in the operator (rounded matrix, as `classes.ts`).

### P7: ASA world grid

- Mechanism: per contributing atom `structure.lookup3d.find` (unit-boundary query, inverse transform per close unit,
  per-unit 32-per-cell grid incl. H, water, non-polymer, skipped afterwards), radius fixed 6.54 Å; neighbour
  coordinates through `conformation.x/y/z` (a matrix transform each on assembly copies).
- Evidence (`asa.ts`, re-run; radii assignment is 6–10 ms, not a bias): flat-32 869–944 → 674–703 ms; asm-32 932–952 →
  670–678 ms; 7qpd 45–49 → 43–53 ms; max area difference 0. asm-192 profile: lookup 48 % of 5.9 s.
- Fix: in `computeArea`, one world `GridLookup3D` over atoms with `atomRadiusType != 0` (Float32 xyz, serial index,
  vdw), 4 Å cells; query `probe + vdw1 + probe + max(VdWLookup)` after radii assignment (fixes B6); read neighbours from
  the grid arrays. Update `properties.md §2` ("neighbours prefiltered with `structure.lookup3d`").
- Risk: +20 B per contributing atom during the computation; B6 changes areas slightly near Se / metal ions (intended).

### P13: inter-contact storage

- Mechanism: `InteractionsInterContacts` extends core `InterUnitGraph` (`inter-unit-graph.ts:61-95`): per-pair maps of
  `{indexB, props}`, a spread `Edge` per edge, `edgeKeyIndex` (Map of Maps), `vertexKeyIndex` (Map → `number[]`), plus
  this class's `elementKeyIndex` (`Map<string, number[]>`, `` `${index}|${unitId}` `` per member, `:93-121`).
- Evidence (`retained.ts`, full GC, re-run): flat-32 24.0 MB retained per computation, inter contacts 20.3 MB for
  41,864 edges = 508 B/edge (asm-32 515). 13 measured 1,143 MB retained for interactions at flat-864. Build on asm-192:
  `InteractionsInterContacts` 67 + `_InterUnitGraph` 41 + `getContacts` 73 ms (~8 %).
- Fix: typed columns (unitA, unitB, indexA, indexB Int32; type, flag Int8) sorted by (unitA, indexA) with CSR offsets
  per unit pair; numeric element keys (per-unit offset arrays or `cantorPairing`). Keep an `edges[i].props`-shaped view
  or migrate the readers (refine.ts, the interactions representation, the extension).
- Risk: API change for `contacts.edges[i].props`; refiners mutate `flag` (writable typed array).

### P2: feature cache key

- Mechanism: `cacheKey = JSON.stringify(p)` (`:356`) spans every tester, contact and bridge param; features depend
  only on the provider set; each key leaves a `features-…` entry on `unit.transientCache` for the unit's lifetime.
- Evidence: code unchanged; v1: +0.5 MB per change on 7qpd, +13 MB on 186k atoms. `findUnitFeatures` 351 ms (7 %) of
  flat-192 interactions, redone per param change.
- Fix: key by the sorted provider types (plus P10's pruning set and B2's topology identity), one slot per unit under a
  fixed symbol holding `{ key, features }`. Document "providers are param-free" in `properties.md §2`.

### Low

- **P11** `checkLineOfSight` runs `structure.lookup3d.find` (3 Å × `distFactor`) at the midpoint of every accepted
  contact, intra and inter, then `typeSymbol` + `VdwRadius` + altLoc per hit, H skipped. 534 ms (11 %) of flat-192, 269
  ms (12 %) of asm-192 (`_find` 490 / 244 ms). Fix: a lazily built world grid of heavy atoms (Float32 xyz, vdw, unit
  and element back-refs) shared with P7's builder. Gain not measured: plausible ~half.
- **P8** `calcValenceModel` caches on each unit; `applyOperator` copies have their own cache, so asm-192 computes 1728
  units for 9 groups (0.32 s / 484 MB). Per atom a 4-tuple (`:334`), `forEach` closure (`:83`), closures in
  `bondToElementCount`, 5–7 `getEdgeIndices` whose misses allocate `[]`. Fix: write into the output arrays, inline the
  loops, skip inter-unit lookups when `getConnectedUnits(unit.id)` is empty, shared empty array in core (02/P9); group
  sharing with copy-on-write for atoms with inter-unit bonds (test on asm-32).
- **P12** `checkGeometry` per candidate: `calcAngles` two arrays + closure (twice), `calcPlaneAngle` `[Vec3(), Vec3()]`
  + two closures, `closestHydrogen` an object + closure (#2004), the `donor` object (`:252`), `.some` closures,
  `[don, acc]` arrays (`:294`, `:325`). 141 MB of 539 MB sampled on flat-32; 233 ms (5 %) CPU on flat-192. Not reduced
  by P4/P10 (real donor-acceptor pairs). Fix: module scratch, `Bond.ElementBondIterator` loops, out-param objects,
  if/else instead of destructuring.
- **13/P6** `computeArea` runs over every serial atom of the root. Interface share (`asa.ts`, re-run): 84.9 % on
  flat-32, 18.5 % on 7qpd, so per-group reuse saves ≤ 15 % on 1tup assemblies. After P7: compute the representative
  against its own unit, flag atoms with a neighbour in another unit, copy the rest to members.
- **P6** `GridLookup3D(position, boundary)` without a cell size: ~8 Å cells. Prebuilt grids on flat-32: 4 per cell
  173–179 vs 191–233 ms. Pass 4 together with P5.
- **P9** the intra refinement iterates every unit id of `unitsContacts` though members share one object; template-string
  keys and `[dist, edge]` tuples (`:123-124`), two `Features.Info` per inter edge (`:47-50`), string keys and a new `[]`
  per miss in `getContactIndicesForElement` (`common.ts:95`). asm-192 `refineInteractions` 84 ms (4 %). Refining each
  object once changes N2's outcome: decide N2 first.
- **12/S1** `traceOneDirection` creates a `Map` per call (`:64`) and a `Vec3` per step before the tolerance test
  (`:90`). 12 estimated ≤ 10–20 % of a 38 ms trace. Fix: a reusable typed visit grid; allocate after the test (with B10).

## Bugs

| Id | vs v1 | Severity | Title | File | Verdict |
| --- | --- | --- | --- | --- | --- |
| B12 | new; **fixed** #2032 | medium | Weak-H-bond refiner throws when intra contacts are skipped (interactions extension) | [refine.ts:165](../../packages/model/src/props/computed/interactions/refine.ts#L165) | confirmed |
| B2 | still open | medium | Interaction features cached on shared Units but computed from the local substructure's bonds | [interactions.ts:406](../../packages/model/src/props/computed/interactions/interactions.ts#L406) | confirmed |
| B3 | still open; **fixed** #2032 | medium | ASA throws on any structure with coarse units | [radii.ts:39](../../packages/model/src/props/computed/accessible-surface-area/shrake-rupley/radii.ts#L39) | confirmed |
| B11 | new; **fixed** #2033 | low | Every non-backbone N of a non-polymer residue becomes a `NegativeCharge` feature | [charged.ts:156](../../packages/model/src/props/computed/interactions/charged.ts#L156) | confirmed |
| N2 | still open | low | Refiners write per-member decisions into a group's shared intra flags | [refine.ts:190](../../packages/model/src/props/computed/interactions/refine.ts#L190) | confirmed (code) |
| B4 | still open; **fixed** #2034 | low | `calcPlaneAngle` skips the second-neighbour search when the neighbour is UnitIndex 0 | [geometry.ts:147](../../packages/model/src/props/computed/chemistry/geometry.ts#L147) | confirmed |
| B6 | still open; **fixed** #2044 | low | ASA lookup radius fixed before radii > 1.87 Å are appended | [shrake-rupley.ts:90](../../packages/model/src/props/computed/accessible-surface-area/shrake-rupley.ts#L90) | confirmed |
| B7 | still open; **fixed** #2034 | low | Water bridges ignore `validUnitPair` and `unitPairTest` | [water-bridges.ts:211](../../packages/model/src/props/computed/interactions/water-bridges.ts#L211) | confirmed |
| B5 | still open; **fixed** #2047 | low (latent) | `getArray` drops the first object-form MolScript argument | [table.ts:632](../../packages/model/src/script/runtime/query/table.ts#L632) | confirmed |
| B8 | still open; **fixed** #2047 | low (latent) | DSSP backbone H-bond graph stores energies in input order | [backbone-hbonds.ts:153](../../packages/model/src/props/computed/secondary-structure/dssp/backbone-hbonds.ts#L153) | confirmed |
| B10 | still open | low | Streamline tracer's `prev` is module state, never reset | [basic.ts:42](../../packages/model/src/props/volume/streamlines/basic.ts#L42) | confirmed |
| B9 | still open; **fixed** #2045 | low | `CustomStructureProperty.attach` has no in-flight dedup | [custom-structure-property.ts:76](../../packages/model/src/props/common/custom-structure-property.ts#L76) | plausible |

- **B12** With `skipIntraContacts`, `unitsContacts` stays empty (`interactions.ts:389`), but `refineInteractions` still
  runs `weakHydrogenBondsRefiner.handleInterContact`, whose `hasHydrogenBond` destructures
  `interactions.unitsContacts.get(acc.unit.id)` (undefined). Failure (`b12.ts`, re-run on 7qpd): weak H-bonds on +
  `{ skipIntraContacts: true }` → "TypeError: Cannot destructure property 'edgeProps' of
  'interactions.unitsContacts.get(...)' as it is undefined." (the full computation finds 120 inter edges). Only caller:
  `packages/extensions/interactions/src/compute.ts:52`, which forwards the user's interaction options, so the extension crashes
  whenever weak H-bonds are enabled. The other refiners read `unitsContacts` only from `handleIntraContact` (loop over
  its keys) and are safe. Fix: skip the intra check when there is no entry; test `computeInteractions` with
  `skipIntraContacts` and every provider on.
- **B2** `findUnitFeatures` runs providers against the local structure (`eachBondedAtom`, `bondCount` see that
  structure's `interUnitBonds`) and stores the result on `unit.transientCache`, shared when a component holds the whole
  unit (`Unit.getChild` returns the same object); the extension also copies units with `propagateTransientCache: true`
  (`compute.ts:37`). Failure (v1 probe, code unchanged): 7qpd units 1 and 3 have one feature more or less depending on
  whether a component or the root computed first. Fix: providers on the root topology as `calcUnitValenceModel` does
  (`valence-model.ts:310-321`), or the topology identity in the key (with P2).
- **B3** `assignRadiusForHeavyAtoms` calls `residue.key` / `label_atom_id` on every unit; coarse units throw "Property
  only available for atomic models."; `isApplicable` is `() => true` (`accessible-surface-area.ts:62`). Failure: 8zzh
  and mixed IHM entries: `asa`, the ASA theme and `is-buried` / `is-accessible` throw (00-profile). Fix: non-atomic
  units get `atomRadiusType` 0 and `serialResidueIndex` -1, advance `m`; `isApplicable: s => s.units.some(Unit.isAtomic)`.
- **B11** `addUnitNegativeCharges`' non-polymer branch (`:155-160`) does `startState`, `pushMember` for a non-backbone
  N, `finishState(NegativeCharge)` per atom: a copy of the positive-charge code. Failure (`negn.ts`, re-run): 7qpd has 4
  single-N negative features, all NAG N2 (amide); with ionic on they make "ionic" contacts with Lys/Arg/His, and since
  #2004 `saltBridgeRefiner` then hides H-bonds between the same groups. Ionic is off by default. Fix: delete the block;
  keep the acid groups and `charge[j] < 0` atoms; test on a NAG or amine ligand.
- **N2** members of a symmetry group share one intra-contacts object (`interactions.ts:391-395`), refined once per
  member. `weakHydrogenBondsRefiner.handleIntraContact` decides from the member's own inter contacts (`:174`) and
  writes the shared flag; `waterBridgeRefiner` (`:438-443`) does the same with legs keyed by the member's unit id.
  Failure: in a symmetry range, a weak H-bond hidden by one copy's crystal-contact H-bond disappears in every copy.
  Fix: refine intra flags from intra data only, or per-member flag arrays when an inter-dependent refiner fires.
- **B4** `if (ni === 1 && unitX1 && indexX1)`: index 0 is falsy. Failure: an acceptor/donor whose only heavy neighbour
  is atom 0 of its unit (e.g. acetate O–C, C first in an ACT unit) gets `undefined` and skips the out-of-plane test.
  Fix: `indexX1 !== undefined`.
- **B6** `maxLookupRadius = 2·probe + 2·VdWLookup[2]` (6.54 Å) is set in `initialize`; `handleNonStandardCase` then
  appends larger radii (Se 1.9; metals with `nonPolymer`). Failure: pairs at d in (6.54, 2·probe + r1 + r2) are missed;
  ASA overestimated near Se / ions. Fix: radius after radii assignment, or per atom with P7.
- **B7** `findWaterBridgeContacts` accepts every `structure.lookup3d` hit. Failure: an ensemble shown as all models: a
  water of model 1 bridges atoms of model 2; the extension's contact groups (`unitPairTest`) are ignored for bridges.
  Fix: skip hits with `!validUnitPair(structure, unitW, hitUnit)`, require it between donor and acceptor, pass
  `unitPairTest` through.
- **B5** object-args branch loops from `i = 1` without pushing `keys[0]` (the `math.sub` code at `:121-123` uses
  `keys[0]` as its initial value; `getArray` does not). Failure: `ring-fingerprint` with named / JSON args loses its
  first element (the builder emits arrays, so latent). Fix: start at 0.
- **B8** `createGraph({ energies })` passes the input-order array; slot-ordered `_energies` is discarded. No reader of
  `edgeProps.energies` today. Fix: `{ energies: _energies }`.
- **B10** `prev` (module `Vec3`, starts at 0,0,0) persists across lines and directions. Failure: a new line whose first
  midpoint lies within `stepSize / 1.1` of the previous line's last written point loses that point. Fix: set `prev` to
  NaN at the start of `traceOneDirection`.
- **B9** `attach` returns early only when a value exists. Two concurrent attaches both run `obtain` (interactions
  1–5 s, ASA 1–6 s) and the last `set` wins, possibly with the other caller's props. No concurrent caller found (state
  transforms are sequential). Fix: pending-promise map per (structure, descriptor, props key), cleared on settle.

## Checked and dropped

Nothing refuted. Corrections to the auditor's notes:

- ASA "5.7 s on flat-192" with "lookup 42 %": the profile is asm-192 (5.9 s; lookup 48 %); 13's table gives
  5.6 s on both flat-192 and asm-192.
- P5 "below 250k points the 98-direction helper": ≤ 10k uses 98 directions, 10k–250k 14, > 250k `getFastBoundary`;
  feature sets (~650 per unit) take the 98-direction path. "#1996 also builds A-side grids": only for A units (the
  smaller of the pair, `structure.ts:1274`) that are never a B and fail the all-candidates shortcut: minor.
- P4 "same edge order": now verified edge by edge (was only counts and type histograms).
- P3 impact: the search share after P4+P10 is measured (43 %); full reuse remains a spec decision.
- P11 impact medium → low (an M-effort fix for ~5 %, unmeasured).
- v1 P3's "safe variant" (units without inter-unit bonds): dropped, covers 11 % of asm-192 units.
- `VdWLookup` growth across ASA runs: not a leak (`indexOf` before append, bounded by distinct radii).

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/1f4e2fbb-da6c-4948-aacb-385a4dab465c/scratchpad/perf-audit-2/model-props-script/`,
built with `node build.mjs <in.ts> <out.mjs>` (esbuild, `molstar-src` condition); big fixtures need
`--max-old-space-size=8000`, `retained` needs `--expose-gc`. Fixtures in `load.ts` (7qpd, 1tup-2x2x2, flat-32/192,
asm-32/192 from `tests/build/pipeline/fixtures/`).

| Script | What |
| --- | --- |
| `inter.ts` | inter-loop variants: per-type radius + mask (P4), exact pairs + present types (P10), 4 per cell and fast boundary (P5, P6); counts, timings |
| `v/inter-eq.ts`, `v/inter-eq2.ts` | edge-by-edge equality of P4 and P4+P10 vs `addStructureContacts`; search-only share (P3) |
| `intra.ts`, `intra2.ts` | grid boundary cost (P5); intra time per tester (P10) |
| `asa.ts`, `v/asa-radii.ts` | world grid vs structure lookup with area diff (P7); interface share (13/P6); radii-assignment cost |
| `classes.ts` | unit pairs, symmetry classes, units with inter-unit bonds (P3) |
| `retained.ts` | retained memory per interactions computation (P13) |
| `negn.ts`, `b12.ts` | B11, B12 reproductions |
| `alloc-flat-32/`, `pf192/`, `pa192/` (+ `.txt`) | allocation profile flat-32; CPU profiles flat-192 and asm-192 |
