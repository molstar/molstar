# Plan: Kekulé assignment of aromatic sets inside Hückel

Status: proposed, not started. The [benchmark](benchmark.md) and its baseline (phase 1) come before any change to the
algorithm (phase 2).

## 1. The gap

`perceiveAromaticRings` (`huckel.ts`) types each atom of a ring set (π contribution 0, 1 or 2) and checks with
`feasibleRingKekule` that the 1e atoms have a perfect matching over the ring bonds. Only part of the result reaches
`assignBondOrders` (`kekule.ts`, step 8):

| Passed on                                                          | Lost                                                       |
| ------------------------------------------------------------------ | ---------------------------------------------------------- |
| `Aromatic` flag on ring bonds (read by `isAromaticBond`)           | which atoms are 1e, i.e. must take exactly one ring double |
| 0e/2e atoms: ring bonds as assigned singles, C=O, pyrrole `hCount` | the perfect matching that `feasibleRingKekule` found       |

Step 8 recomputes demand from geometry (`isSp2OrAmbiguousDoubleDemand`) and maximizes a score instead of requiring a
perfect matching. Consequences:

- **G1 coverage**: the only hard constraint is `mustCoverDouble` (Trigonal C, degree 3; its TODO names this gap). A
  pyridine-type N, a CH ring carbon (degree 2 without H), an ambiguous atom can stay without a double.
- **G2 later passes override the typing**: an uncovered 1e atom is open to step 8b (`assignCarbonylsByDistance`: an enol
  `*C(O)*` resolved to 1e becomes a ketone) and step 9a (`isRingBond`: a double into a neighbouring non-aromatic ring).
  The π count no longer matches what Hückel checked; the bonds keep `Aromatic`.
- **G3 per-ring fallback checks rings independently**: when a system fails as a whole, each ring is checked for its own
  perfect matching as if the fusion atoms were free in both. Two rings can each pass while their union has none.
- **G4 stale `open`**: `computeOpenValence` runs once per residue in `perceiveAromaticRings`. Only one write reaches a
  sibling ring: step 8d charging a bridgehead N (typed `*N(*)*` within one ring, the other ring's bond counting as
  exocyclic). The sibling then types it 1e but `feasibleRingKekule` reads its old `open = 0`.
- **G5 no Trigonal assignment**: the `huckel.ts` header says aromatic atoms become sp2; the code never sets `geometry`.
  Ambiguous aromatic atoms can be turned Tetrahedral by `preventConjugationNextToPyramidalNitrogen` and are skipped by
  `isSp2DoubleDemand` (step 9d).
- **G6 Kekulé form depends on coordinate noise**: in step 8, `edgeWeight` is 2 for every aromatic edge plus a length
  bonus (0.1/0.25 Å buckets). Aromatic bond lengths carry no Kekulé signal, so the written form follows the residue's
  coordinate error; two copies of a ligand in different models can differ.

Verified, not a gap: a failed `perceiveAromaticSet` writes nothing to `State` (all writes follow the last
`return false`). The one earlier write, step 8d's `charge = 1`, always leaves `sum % 4 === 2`, so the set then succeeds.
This holds by arithmetic only; phase 2 makes it structural.

Out of scope: sub-systems between "whole system" and "one ring" (a bicyclic aromatic fused to a saturated ring is judged
ring by ring); systems above 60 demand atoms (fullerenes and graphene-like fragments).

## 2. Phase 1: baseline

The [benchmark](benchmark.md) run on the current code, its baseline committed. Then the cases that pin this change,
found in the benchmark's output where possible:

| Case                                                            | Gap         | Where to find it                               |
| --------------------------------------------------------------- | ----------- | ---------------------------------------------- |
| aromatic set violating the 4n+2 / 1e diagnostic                 | G1, G2      | benchmark diagnostics (`benchmark.md §5.2`)    |
| aromatic ring fused to a saturated ring (per-ring fallback)     | fallback    | indane, tetrahydrocarbazole components         |
| bridgehead N in a fused system                                  | G4, step 8d | indolizine, imidazo[1,2-a]pyridine components  |
| two rings each feasible alone, union without a perfect matching | G3          | to be found: tier A failures in fused systems  |
| ambiguous ring atom next to a pyramidal N                       | G5          | to be found                                    |
| Kekulé form under coordinate noise                              | G6          | a one-off check, outside the benchmark (below) |

G6 one-off: perceive components again with their coordinates shifted at random (seeded, ±0.02 Å, below the coordinate
error); a Kekulé form that changes while the InChI does not is G6. Not part of the benchmark.

Each case found becomes a small quick-lane fixture in `src/_test/` (CCD ideal coordinates, H stripped), written as
`it.fails` where the current code gets it wrong, so phase 2 turns it into `it`. One white-box test: a failed
`perceiveAromaticSet` leaves `State` unchanged (export it from `huckel.ts`; the module is not in the package exports).

## 3. Phase 2: the change

1. `feasibleRingKekule` returns the matching it found (pairs of set atoms), or `null`.
2. `resolveAromaticDOF` keeps the unweighted feasibility check per mask; the matching of the accepted mask is kept.
3. `perceiveAromaticSet` writes, in one place after every check: the `Aromatic` flags, the 0e/2e enforcements, the
   matched doubles as assigned bonds, Trigonal geometry for every set atom (ambiguity cleared), and the step 8d charge.
   No write before the set is accepted.
4. `perceiveAromaticRings` calls `computeOpenValence` after each accepted set.
5. Above 60 demand atoms, `feasibleRingKekule` returns `null` ("not aromatic"), so the rings go through the per-ring
   fallback, then step 9a. No set is accepted unchecked.
6. `assignBondOrders` drops step 8 (`matchDoubleBondsBy(..., isAromaticBond)`) and `isAromaticBond`. `mustCoverDouble`'s
   TODO goes: aromatic atoms are covered by construction.
7. Docs: `huckel.ts` and `kekule.ts` headers, the `_spec` file table (§2: kekule no longer kekulizes aromatic rings),
   `_spec §4` (the new tests), this item out of `TODO.md`, the baseline updated.

Kekulé choice: the first perfect matching found, a function of the topology and atom order only (G6). No length
weighting for aromatic sets.

## 4. Acceptance

- The benchmark's compare report against the baseline: no regression in any tier, or each explained. A moved Kekulé
  double inside an aromatic set (L0 only) is expected; a change at L2 is not, unless explained.
- The 4n+2 / 1e diagnostic reports no aromatic set.
- The `it.fails` cases pass as `it`.

## 5. Open questions

- A ring pair for G3 and a case for G5: none known yet.
- The diagnostic's π reading of unusual motifs (`*N(=O)*` N-oxide, `*[N+](*)=*`): follow `AromaticAtomConfigs` exactly,
  or share the table with the benchmark.
