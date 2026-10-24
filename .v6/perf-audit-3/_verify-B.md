# Verify B: slices 04, 05, 06 (adversarial pass)

HEAD `16277d822`. Read-only. Load fell from ~30 to 5–8 during this pass, so some reruns are quieter than the
auditors' numbers; ratios are from one process. New probes: `pa3/vB/poly.ts` (polymer.ts test), `pa3/vB/ranges.ts`
(`areBackboneConnected` allocation). Every other number comes from rerunning the auditor's own probe.

| Id | Verdict | Severity / estimate | Reason |
| --- | --- | --- | --- |
| 04/P20 | CONFIRMED | as stated: +12 B/atom/frame (−1.2 GB at 100k × 1000) | `lmp.mjs` rerun: 12.7 vs 24.3 B/atom/frame. The probe's +12 B is exactly the three Float32 `ux/uy/uz`. Only `model.ts:143` and `periodic.ts` read them, and lazily |
| 04/P21 | CONFIRMED | as stated (rerun 1k 18, 10k 199, 30k 850 ms) | linear scan from the start per frame. A DCD with `NSAVC > 0` always sets `step`. Reached through `hierarchy.ts:65/77` |
| 04/B15 | CONFIRMED | low | `seqsort.mjs` rerun: `ALA1 GLY2 ALA1 HOH0`. `aI` is a position in `indices`, and `valueKind` reads the row at that position |
| 05/B15 | CONFIRMED, severity CORRECTED | low (was medium) | `stale.mjs` rerun: the box model gets the same helix-orientation object, centre 35.09 vs 5.09. It shows only with `tubularHelices` (off by default) plus a periodic-mode switch |
| 05/P23 | CONFIRMED | as stated (100k: 83 MB kept after drop; linear) | module-level strong `Map`, never cleared. `retain.mjs 100000`: 24 → 83 MB, and it stays at 83 MB after every reference is dropped |
| 05/P24 | CONFIRMED, estimate CORRECTED | ~0.25 s once per root on asm-192 (12–60×), not 1.1–1.8 s | `xform2.mjs` at load 5: root route 239–307 ms vs own 4–21 ms. asm-32: 38–51 vs 4–10 ms. Edges 0 = 0 |
| 06/P14 | CONFIRMED, scope narrowed | low; opt-in only | it applies only with a radius above 1.87: `nonPolymer: true` (default false, no in-repo caller sets it) with metals, or polymer Se (+~3 % volume). The per-atom bound `2p + vdw1 + maxR` is correct |
| 06/B13 | CONFIRMED | low | `b13.mjs` rerun: members `5:O5 2:O2 3:O3 1:O1`. The `(_, k)` callback drops the unit. Both `eachIntraBondedAtom` and skipping `unitX !== unit` are correct fixes |
| 05/P13 | PLAUSIBLE, estimate CORRECTED | asm-192 ÷2.8–3.8 (was ÷4) | `symbonds.mjs` reruns: 137 → 48 ms and 329 → 87 ms, edges 21,120 = 21,120. The crude prototype allocates more (66 vs 37 MB/run) |
| 05/P18 | CONFIRMED | as stated (quadratic) | `conn.mjs` rerun: fixed point 0.97 s, layers 50→800: 3, 8, 29, 108, 449 ms (×4 per doubling). Reached from `manipulate.ts:102` |
| 05/P21 | PLAUSIBLE | as stated (bonds-intra −25–40 %, upper bound ~−43 %) | `pairs-half.mjs` rerun: 204 vs 73 ms, same 563,456 pairs. `findBonds` 301 ms drops half its hits (`bI <= aI`). Per-atom-A state (threshold, structConnAdded) has to move into the pair loop. Edge order may change |
| 06/P3 | PLAUSIBLE | ~26 % (search only) stands | profile rerun: `findInterUnitContacts` 82 % of `computeInteractions` on asm-32. Grid queries plus candidates add up to about 31 % |
| 06/P11 | share CONFIRMED, gain CORRECTED | low: ≤ ~5 % of interactions per computation, net ~0 unless the grid is shared | `los.mjs asm-32` at load 5: interactions 271–323 ms. The structure-lookup replay is 60 ms per 3 passes vs grid 12–28 ms plus a 20 ms build. The 2.1 Å clamp is exact but gains no time |
| 01 note → 05/P25 | CONFIRMED | low: −40–60 % of `getAtomicPolymerElements` (flat-32 11–12 → 4.5–7 ms; flat-192 ~6×) | `move()` yields only non-empty segments, so the test is always true (0 false on 416 units). Output identical with the test removed. New id **05/P25** (free after P24) |
| 04 note: `areBackboneConnected` 4 MB | REFUTED | not a finding | it is 9 samples of a 512 KB sampling profile. Measured: the whole `getAtomicRanges` 0.9–2.5 MB / 6–10 ms; the backbone lookups alone 0.3–1.3 MB for 32k pairs |

## Evidence notes

- **04/P20**: `lammpsFrameCoordinates` builds `ux/uy/uz` for every frame inside the conversion loop of
  `coordinatesFromLammpsTrajectory`. The fix can stay inside the reader: a getter on the frame object, or the image
  columns plus the box. A thunk-typed `Frame.unwrapped` would change the public `Coordinates.Frame` type
  (`coordinates.ts:35`), and that change would go into `breaking-v6-changes.md`. The parsed file keeps `frame.image`
  anyway.
- **04/P21**: an advancing pointer gives the same output even when entries are not sorted. Frame steps
  (`ISTART + i·NSAVC`) only increase, and the first `k` with `entries[k].step > s` cannot move backwards as `s` grows.
  The `break` already depends on that order.
- **04/B15**: models after the first also read model 1's rows by position (`indices` start at `start`). For NMR models
  with identical layouts the result does not change, so the visible case is out-of-order residues after a moved bucket.
- **05/B15**: in `apply`, 'box' is made from the whole model, so it inherits the whole model's dynamic data
  (`HelixOrientationProvider`, type `'dynamic'`). 'whole' inherits the source frame's, so off → whole is stale too. Other
  model-level dynamic providers (rcsb validation) do not depend on coordinates.
- **05/P23**: the proposed `WeakMap<AtomicHierarchy, LRU>` is sound. Made models spread `...model`, so they share the
  key, and the ephemeron lets them point back at it.
- **05/P24**: the probe's own compute uses the getter's `ignoreWater` and `ignoreIon` defaults. It is reached only after
  a `Structure.transform` (superposition, MVS transforms), and the cost is paid once per root.
- **06/P14**: base set `2p + 2·VdWLookup[2]` (6.54 Å) before the radii were appended. HEAD's correct-but-wider bound
  matters only when a radius above 1.87 is assigned.
- **06/P11**: in the auditor's loaded run the replay looked like a 2–4× win. At low load, building the world grid
  (20 ms) costs more than it saves for a single computation (~13–16 ms/pass). So the rank "medium" holds only if the grid
  is shared with ASA (P7) or kept across recomputes.
- **05/P25** (new, from slice 01): `polymer.ts:25`
  `OrderedSet.areIntersecting(Interval.ofRange(elements[start], elements[end - 1]), elements)` tests a residue segment's
  own elements against the set. The fix is to drop the test. Effort S, output identical (`vB/poly.ts`: `SortedArray`s
  equal on flat-32). It also saves ~0.8 MB/run of allocation on flat-32.
- **04 note**: `findAtomByNames` reads a dictionary-backed `label_atom_id` (no string allocation). The remaining bytes
  come from iterators, range arrays and an unoptimized first run. This is not 4 MB and does not need an id.
