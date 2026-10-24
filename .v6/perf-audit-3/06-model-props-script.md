# Slice 06: model props and script (v3)

Scope: [`packages/model/src/props`](../../packages/model/src/props) (interactions, ASA, valence model, DSSP,
streamlines, `CustomStructureProperty`) and [`packages/model/src/script`](../../packages/model/src/script). HEAD
`16277d822`, v2 base `0944823d4`. v2 report: [06-model-props-script.md](../perf-audit-2/06-model-props-script.md).

Measurement: Apple M4, load 3 at the start, **23–31 for most of the run** (other auditors). Wall times are 2–5× v2's and
swing ±50 %: base and HEAD were built from the same bench (`git archive 0944823d4 packages` in scratch, esbuild with the
`molstar-src` condition) and run back to back, but only allocation (`total_allocated_bytes`), retained MB, counters,
hashes and profile shares are reported as reliable. Fixtures: 7qpd, flat-32 / asm-32 (`tests/build/pipeline/fixtures/`),
1tup-2x2x2.

PRs since v2 touching the slice: #2025 (06/P4 + P10 loop part), #2032 (B3, B12), #2033 (B11), #2034 (B4, B7), #2044 (P2,
B6), #2045 (P5, B9), #2047 (P12, B5, B8). Outside the slice but visible here: #2046 (core) and #2047 (model) cut DSSP
and valence allocation.

## Hot paths now

| Stage | Fixture | Base → now | Where the time goes now (CPU profile shares) |
| --- | --- | --- | --- |
| `interactions` (cold features) | flat-32 | 2.5–3.1 → 0.83–1.0 s; alloc 477–484 → 100–102 MB | inter loop 38 %, intra loop 33 %, features 17 %, refine 6 %, inter storage 5 %; **line of sight 21 %** (inter 740 + intra 728 ms of 7.0 s profiled), H-bond geometry 12 % |
| `interactions` (cold features) | asm-32 | 2.4–3.7 → 0.52–0.71 s; alloc 237–257 → 61–81 MB | inter loop 72 %: line of sight 24 %, feature-grid search ~26 %, H-bond geometry 9 %; inter storage 9.5 %; refine 8 % |
| `asa` | asm-32 / flat-32 | unchanged code; alloc 11–13 MB both | `structure.lookup3d` `_find` 59 % of `computeRange` (asm-32), kernel self 36 % |
| `valence-model` | flat-32 / asm-32 | alloc 79 → 52 MB (outside the slice) | still per unit: asm-32 (9 groups) allocates the same 52 MB as flat-32 (416 groups) |
| DSSP (`computeUnitDSSP` per group) | flat-32 | alloc 83–97 → 14–16 MB (outside the slice) | output hash identical |
| MolScript, streamlines, helix orientation | | | not hot; no change but B5 |

- Contacts are **identical to the v2 base**: FNV hash over every inter edge (unitA, unitB, indexA, indexB, type, flag, in
  order) and every intra edge, defaults and all providers + water bridges, equal on 7qpd, flat-32 and asm-32. The only
  difference: 7qpd with every provider on (ionic included) has 4 features fewer (18,876 → 18,872, the NAG N2 of B11's
  fix) and so a different inter hash with the same 200 edges and the same type histogram (intended).
- Per computation retained: 23.7 MB flat-32, of which inter contacts 20.3 MB = 508 B/edge (unchanged, P13). Feature
  cache retained 13.9 MB flat-32 (base 15.0) and 0.6 MB asm-32: the subset memo of #2045 costs nothing net.

## v2 → now

26 v2 items (14 perf, 12 bugs): **12 fixed**, **2 partly fixed**, **12 still open**, 0 regressed, 0 dropped.

| v2 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| P4 | Inter queries at 6 Å for every feature; `validPair` first | fixed | #2025 | `ContactTester.acceptsTypes`, per-type query radius and per-pair max ([contacts.ts:391](../../packages/model/src/props/computed/interactions/contacts.ts#L391)); edges equal to base in order |
| P10 | Feature types with no partner computed, gridded, queried | partly fixed | #2025 | loop part done (inter type-pair mask, intra skip of unpartnered testers and features); providers not pruned: halogen acceptors 110 + metal binding 183 ms = 4 % of flat-32 interactions for 0 contacts |
| P5 | 98-direction `getBoundary` for feature and subset grids; subsets per computation | fixed | #2045 | exact box + centroid sphere, `getFastBoundary` for subsets, subset memo per type set; `getFeaturesBoundary` 1.1 % of flat-32 |
| P3 | Inter contacts per symmetry-equivalent unit pair | still open | | asm-32: search ~26 %, line of sight 24 %, geometry 9 % of interactions; see P3 below |
| P7 | ASA neighbours through the structure lookup | still open | | `_find` 59 % of `computeRange` on asm-32 (v2: 48 % asm-192); world-grid A/B too noisy to re-time (v2: −23–28 %) |
| P13 | Inter contacts retained as objects, Maps, string keys | still open | | 508 B/edge flat-32, 515 asm-32, unchanged; build (`InteractionsInterContacts` + `getContacts`) 5 % flat, 9.5 % asm |
| P2 | Feature cache keyed by all params | fixed | #2044 | key = provider indices ([interactions.ts:374](../../packages/model/src/props/computed/interactions/interactions.ts#L374)); bounded by provider sets |
| P11 | Line of sight through the structure lookup | still open, **now medium** | | 21 % (flat-32) / 24 % (asm-32) of interactions after P4/P5/P10 (v2: 11–12 %) |
| P8 | Valence model per unit, tuple and closures per atom | still open | | code unchanged; asm-32 allocates as much as flat-32 |
| P12 | H-bond geometry allocates per candidate | fixed | #2047 | module scratch, out-params, loops; part of the −380 MB alloc on flat-32 |
| 13/P6 | ASA for every atom of every copy | still open | | interface share 84.9 % flat-32 / asm-32, 18.5 % 7qpd (unchanged) |
| P6 | Feature grids at 32 points per cell | still open | | not re-timed (noise); **changes the inter edge order** (hash differs at 4 and 8 per cell, same as v2's set), so it needs an order decision or a sort |
| P9 | Intra refinement per member; string keys, tuples | still open | | `refineInteractions` 6 % flat / 8 % asm; `Features.Info` per inter edge 46–59 ms self |
| 12/S1 | Streamline Map per line, Vec3 per step | still open | | code unchanged |
| B12 | Weak-H-bond refiner throws with `skipIntraContacts` | fixed | #2032 | intra check skipped when there is no entry |
| B2 | Features cached on shared Units, computed from the local structure | still open | | `findUnitFeatures` unchanged ([interactions.ts:407](../../packages/model/src/props/computed/interactions/interactions.ts#L407)) |
| B3 | ASA throws on coarse units | fixed | #2032 | coarse elements radius 0 / residue -1; `isApplicable` needs an atomic unit; `getNormalizedValue` short-circuits on -1 |
| B11 | Every non-backbone ligand N a `NegativeCharge` | fixed | #2033 | 7qpd: −4 features with ionic on, no other change |
| N2 | Refiners write per-member decisions into shared intra flags | still open | | [refine.ts:196](../../packages/model/src/props/computed/interactions/refine.ts#L196), [:445](../../packages/model/src/props/computed/interactions/refine.ts#L445) |
| B4 | `calcPlaneAngle` skips neighbour index 0 | fixed | #2034 | `indexX1 !== undefined` |
| B6 | ASA lookup radius before appended radii | fixed | #2044 | but costs every atom: P14 |
| B7 | Water bridges ignore `validUnitPair` / `unitPairTest` | partly fixed | #2034 | `validUnitPair` for legs and donor–acceptor; `unitPairTest` still not passed to bridges |
| B5 | `getArray` drops the first object-form argument | fixed | #2047 | from key 0 |
| B8 | DSSP H-bond energies in input order | fixed | #2047 | `_energies`; DSSP output hash equal to base |
| B10 | Streamline `prev` is module state | still open | | [basic.ts:42](../../packages/model/src/props/volume/streamlines/basic.ts#L42) |
| B9 | `CustomStructureProperty.attach` without in-flight dedup | fixed | #2045 | pending map per structure; see review notes |

### New bugs or costs from fixes since v2

- **#2044 (B6) → P14**: the ASA lookup radius is now `2·probe + 2·max assigned radius` for every atom; one Zn with
  `nonPolymer` on raises it 6.54 → 7.0 Å and the hits per atom +19 % on 1tup-2x2x2 (52.8 → 63.0), with the area sum
  unchanged (976,177.5 both). See P14.
- **#2045 (B9)**, low, not a regression: an attach with props different from the one in flight starts its own
  computation and replaces the pending entry; if the older one finishes last, its `set` wins (stale props). Same as
  before the fix. A waiter also does not observe its own abort while it awaits the other requester's promise.
- **#2025 (P4/P10)**: exact in order (hashes above). The `acceptsTypes` contract is documented on `ContactTester`; a
  tester without it is tested on every pair, as before. `getTesterTable` is rebuilt per computation (new testers array):
  16 ms on flat-32, negligible.
- **#2045 (P5)**: centroid sphere ≥ every point, box exact; subset memo retains nothing net (13.9 vs 15.0 MB).
- **#2047 (P12)**: module scratch in `geometry.ts` / `hydrogen-bonds.ts` is not reentrant, but no caller nests (water
  bridges call `checkGeometry` sequentially); the scratch drops its unit references after each call.
- No new cost from #2032, #2033, #2034 (B7 adds one `validUnitPair` per water-lookup hit; bridges are off by default).

## Open performance findings (ranked)

| Id | Status | Title | File | Kind | Estimate | Effort | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- |
| P3 | still open | Inter contacts recomputed per symmetry-equivalent unit pair | [interactions.ts:464](../../packages/model/src/props/computed/interactions/interactions.ts#L464) | time | measured shares, asm-32: candidate reuse removes the search, ~26 % of interactions; reusing line of sight and geometry too (a spec decision) ~60 % | L | no |
| P11 | still open (low → medium) | Line of sight queries the structure lookup at 3 Å per accepted contact | [contacts.ts:100](../../packages/model/src/props/computed/interactions/contacts.ts#L100) | time | measured share 21–24 % of interactions; midpoint replay: structure lookup 330–574 ms → heavy-atom world grid 37–178 ms per 3 replays, build 136–284 ms (net only when shared with P7 or reused) | M | no |
| P7 | still open | ASA neighbours through the structure lookup | [area.ts:98](../../packages/model/src/props/computed/accessible-surface-area/shrake-rupley/area.ts#L98) | time | lookup 59 % of `computeRange` (asm-32); v2 measured −23–28 % with a world grid, areas identical | M | no |
| P13 | still open | Inter contacts as objects, three Maps, string keys: 508 B/edge | [common.ts:92](../../packages/model/src/props/computed/interactions/common.ts#L92) | memory+time | retained ÷10 (measured 20.3 MB / 41,864 edges flat-32); build 5–9.5 % of interactions | M | no |
| P9 | still open | Intra refinement per member; `Info` and string keys per inter edge | [refine.ts:59](../../packages/model/src/props/computed/interactions/refine.ts#L59) | time | `refineInteractions` 6–8 % of interactions (measured share) | S | no |
| P10 | partly fixed | Halogen-acceptor and metal-binding features built with no donor / metal anywhere | [interactions.ts:372](../../packages/model/src/props/computed/interactions/interactions.ts#L372) | time+memory | 4 % of flat-32 interactions (measured share), 62 % of the features (v2) | M | no |
| P8 | still open | Valence model per unit, not per group | [valence-model.ts:366](../../packages/model/src/props/computed/chemistry/valence-model.ts#L366) | time+memory | asm-32 52 MB alloc = flat-32's; ÷(units/groups) with group sharing (estimated) | M | no |
| P14 | **new** (from #2044) | ASA lookup radius `2·probe + 2·max radius` for every atom | [shrake-rupley.ts:102](../../packages/model/src/props/computed/accessible-surface-area/shrake-rupley.ts#L102) | time | `nonPolymer` + Zn: +19 % hits per atom (measured); Ca / K: ×1.46 / ×2.04 (volume ratio, estimated); per-atom bound `cutoff1 + maxR` | S | no |
| P6 | still open | Feature grids at 32 points per cell | [features.ts:118](../../packages/model/src/props/computed/interactions/features.ts#L118) | time | v2: −9 % of the inter loop; not exact in edge order (measured: hash differs) | S | no |
| 13/P6 | still open | ASA for every atom of every copy | [area.ts:87](../../packages/model/src/props/computed/accessible-surface-area/shrake-rupley/area.ts#L87) | time | ≤ 15 % on 1tup assemblies (interface share 84.9 %) | L | no |
| 12/S1 | still open | Streamline trace: Map per line, Vec3 per step | [basic.ts:64](../../packages/model/src/props/volume/streamlines/basic.ts#L64) | time+memory | 10–20 % of a ~40 ms trace (v2 estimate) | S | no |

### P3: symmetry-equivalent pairs (the largest left)

- After #2025 the inter loop on asm-32 (profile, 5 computations, 3.9 s in `computeInteractions`) is 2.78 s: grid
  queries on feature grids ~0.98 s + `getCandidates` 0.22 s + loop self 0.21 s (the search), line of sight 0.92 s,
  H-bond geometry 0.35 s, `validPair` 0.08 s. v2's classes: 4,282 pairs in 265 classes on asm-32, 33,814 in 267 on
  asm-192.
- Fix as v2: per class the (i, j, distSq) candidate list after the type mask, replayed per pair with `validPair`, the
  testers and line of sight. Exact (rigid operator); needs an equality test like this audit's hash on asm-32 / asm-192.
  Reusing line of sight and geometry per class too (~60 %) needs a spec decision: both see third units and inter-unit
  bonds.

### P11: line of sight

- Counts (`los.ts`): flat-32 41,300 calls, 268,364 hits at 3 Å, of which 96,808 (36 %) lie within the largest heavy-atom
  vdW radius present (2.1 Å, Zn); asm-32 21,584 calls, 112,627 → 47,704 hits. Only atoms closer than `vdw ·
  distFactor` can block, so querying at `maxVdwPresent · distFactor` is exact (hashes equal in the probe) but **no
  time win measured**: the cost is the unit-sphere query, the inverse transform and the cell scan of 32-per-cell unit
  grids, not the hits.
- Midpoint replay (same process, 3 passes): structure lookup at 3 Å 574 (flat-32) / 331 ms (asm-32); one world
  `GridLookup3D` of heavy atoms at 2–4 Å cells: 150–424 / 37–136 ms depending on cell and radius; building it 136–284 ms
  under load. Net positive only if the grid is built once per structure and shared with ASA (P7) or kept across
  recomputes (param changes, the interactions extension's groups).
- Fix: a lazily built heavy-atom world grid on the structure (Float32 xyz, vdW, unit and element back-refs, ~2–3 Å
  cells) used by line of sight and ASA; query at `maxVdwPresent · distFactor`.

### P14: ASA lookup radius since #2044 (new)

- `getMaxLookupRadius` ([shrake-rupley.ts:102](../../packages/model/src/props/computed/accessible-surface-area/shrake-rupley.ts#L102))
  takes the largest radius assigned and applies `2·probe + 2·maxR` to every atom. Probe `asarad.ts` (structure lookup
  wrapped): 1tup-2x2x2 `nonPolymer: true`: radius 6.54 → 7.0 Å, hits per call 52.8 → 63.0 (+19 %), area sum 976,177.5
  before and after (the missed pairs B6 described contribute nothing here); 7qpd unchanged (appended radii 1.7, 1.55,
  1.52 < 1.87).
- The needed bound per atom is `cutoff1 + maxR = 2·probe + vdw1 + maxR` (≤ 6.77 Å for C next to Zn, −11 % hits vs
  7.0). Fix: per-atom radius in `computeRange`, or fold into P7 (the world grid query is per atom anyway).

### Low

- **P13**: unchanged from v2 (typed columns with CSR per unit pair, numeric element keys).
- **P9**: `Features.Info` is built twice per inter edge in `refineInteractions` (46–59 ms self); string keys and
  `[dist, edge]` tuples ([refine.ts:105](../../packages/model/src/props/computed/interactions/refine.ts#L105)). Decide
  N2 first.
- **P10**: run the anchor providers (halogen donors, metals) first; skip `HalogenAcceptorProvider` /
  `MetalBindingProvider` when no unit has an anchor; the provider set is already the cache key (#2044).
- **P8**, **13/P6**, **12/S1**: unchanged from v2.
- **P6**: 4 points per cell reorders hits within a query, so inter edges change order (counts and types unchanged);
  either accept an order change in the spec or sort edges per pair.

## Bugs

| Id | Severity | Status | Title | File | Failure scenario |
| --- | --- | --- | --- | --- | --- |
| B2 | medium | still open | Features cached on shared Units but computed from the local substructure's bonds | [interactions.ts:407](../../packages/model/src/props/computed/interactions/interactions.ts#L407) | 7qpd units 1 and 3 have one feature more or less depending on whether a component or the root computed first (v1/v2 probe; code unchanged) |
| B13 | low | **new** | Ligand acid-group negative charges push bonded O of another unit with the own unit's index | [charged.ts:171](../../packages/model/src/props/computed/interactions/charged.ts#L171) | phosphate ligand covalently bonded to a Ser OG in another chain, ionic on: the phosphate feature gets the ligand's own atom at the OG's index (O5, not bonded to P), centre shifted |
| N2 | low | still open | Refiners write per-member decisions into a group's shared intra flags | [refine.ts:196](../../packages/model/src/props/computed/interactions/refine.ts#L196) | in a symmetry range, a weak H-bond hidden by one copy's crystal contact disappears in every copy |
| B7 | low | partly fixed | Water bridges ignore `unitPairTest` | [water-bridges.ts:224](../../packages/model/src/props/computed/interactions/water-bridges.ts#L224) | the interactions extension's contact groups do not restrict bridges |
| B10 | low | still open | Streamline `prev` is module state, never reset | [basic.ts:42](../../packages/model/src/props/volume/streamlines/basic.ts#L42) | a line whose first midpoint lies within `stepSize / 1.1` of the previous line's last point loses that point |

- **B13** (new). `addUnitNegativeCharges`' non-polymer branch (`:171`) and the nucleotide branch (`:147`) call
  `eachBondedAtom(structure, unit, j, (_, k) => …)`, which also visits inter-unit bonds, and then use `k` with the
  current unit: `typeSymbol(unit, k)`, `x[elements[k]]`, `pushMember(…, k)`. Probe (`b13.ts` on `b13.cif`: SER chain A,
  a 6-atom phosphate ester ligand in chain B, P 1.6 Å from OG; ionic on): P bonds to unit 0 index 5 (OG) and unit 1
  O1–O3; the `NegativeCharge` Phosphate feature has members `5:O5 2:O2 3:O3 1:O1`, so O5 (bonded to C4, 3.6 Å from P)
  replaces OG and the centre moves. When the index is out of the unit's range, `typeSymbol` reads past it and the member
  is skipped by luck. Applies equally to sulfates, sulfonic acids and carboxylates, and to nucleotide phosphates bonded
  across chains. Fix: skip `unitX !== unit` (the bridging O is not a charge carrier), or use `eachIntraBondedAtom` as
  the positive-charge code (`:93`, `:102`) does; test on a covalent phosphate adduct.
- **B2, N2, B10**: code unchanged since v2; scenarios as v2.
- **B7**: `unitPairTest` stays undecided for bridges (bugs.md note of #2034).

## Already good, re-checked

- #2025's contact loops: exact against the base on three fixtures and four provider configurations (hashes above).
- B3: coarse units skipped in radii, `isApplicable` requires an atomic unit, `getValue` returns -1 for them.
- B8: DSSP output hash equal to the base on 7qpd, flat-32, asm-32 (no reader of `energies` in input order).
- Slice tests: `pnpm test:unit packages/model/src/props packages/model/src/script`: 16 files, 51 tests pass.

## Not audited / not measured

- Absolute times: load 23–31 during most runs; ASA and P6 A/Bs within noise; P7 not re-timed.
- B2 not re-probed (code unchanged).
- The source of the DSSP and valence allocation drop (outside the slice: #2046 / #2047).

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/31ee827f-4715-4967-bec3-085f5013d11f/scratchpad/pa3/06/`:
`node build.mjs src/<x>.ts out/<x>.mjs now|base` (esbuild, `molstar-src`; `base/` is `git archive 0944823d4 packages`
with `node_modules/@molstar/*` links).

| Script | What |
| --- | --- |
| `src/ab.ts` | base vs HEAD: interactions (cold/warm, alloc, edge hashes: defaults, all providers without / with ionic), ASA, DSSP, valence |
| `src/prof.ts`, `sum.mjs`, `callers.mjs`, `incl2.mjs` | CPU profiles (`prof/flat32`, `prof/asm32`, `prof/asa32`) and their summaries |
| `src/los.ts` | line-of-sight calls, hits at 3 Å vs max vdW, clamp A/B with hashes, midpoint replay against a world grid (P11) |
| `src/p6.ts` | feature grids at 32 / 8 / 4 per cell, warm, with hashes (P6) |
| `src/asarad.ts` | ASA lookup radius and hits per atom with `nonPolymer` off/on, base vs HEAD (P14) |
| `src/retained.ts`, `src/featret.ts` | retained per computation (P13); feature-cache retained (subset memo) |
| `src/b13.ts`, `b13.cif` | B13 reproduction |
