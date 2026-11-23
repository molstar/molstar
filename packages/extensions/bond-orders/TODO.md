# Bond orders: to do

What is missing, and what was considered and deferred; an entry leaves when it is done (its `_spec` section says so
then). The spec: [\_spec](_spec/README.md).

## Next

- **Benchmark: bonding by distance** (`--connectivity distance`, plan §3.3): `chem_comp_bond` removed, the bond graph
  compared with the CCD's before the orders, a `connectivity` bucket; in the viewer too (plan §7). The rest of the
  benchmark is done (`_spec §1.10`).
- **Kekulé assignment of aromatic sets inside Hückel**: the Hückel typing (which atoms take a ring double) does not
  reach step 8 of `kekule.ts`; the benchmark's baseline first, then the change. [Plan](plans/aromatic-kekulization.md).
- **Flag `OrderUnknown` in the bond-table readers** where the file gives no order: mol2 `un`, LAMMPS, topologies (PSF,
  PRMTOP, TOP), CIF core without `_geom_bond_valence` (core `formats/structure/`). Then `auto` can perceive those
  tables; the table fast path in `perceive.ts` (`_spec §1.9`) stays, for the tables none of whose bonds is
  `OrderUnknown` (checked once per model).

## Later

- **Doubles on borderline CH2 in saturated rings**: with H removed, a ring CH2 (degree 2, ~110°) stays ambiguous when
  one of its bonds passes the dihedral test (5-rings, or just under 15° in a half-chair), and the general Kekulé step
  can place a double on it (59H, 8DA, 1QB: 40 tier A and 78 tier B regressions once the stranded-atom demotion checks
  for unambiguous singles, `hybridization.ts` `markAmbiguous`). The false aromatic set they used to fall into hid it.
  Idea: count the ambiguity status in the Kekulé matching's scoring (`edgeWeight`, `kekule.ts`), so a double between
  trigonal atoms is preferred over one on an ambiguous atom.
- **Benchmark report: no highlight for Kekulé swaps**: the L0 diff highlights every bond whose order differs, the
  swapped doubles of an aromatic ring and of a carboxylate included (197 in tier C: three phenyls lit, the failure is
  the thioacid's C71-O72). Highlight only what the canonical SMILES does not absorb.
- **Perception over the whole PDB**: every entry, not one instance per component; on the benchmark's pipeline.
- **Protonation and charge** (deferred: decide on the benchmark, which shows how often they cause errors; no
  case-by-case rule before then):
  - Functional groups check one multiple bond per neighbour, not its open valence, so a stated H does not keep a
    pattern's double off its atom. Right for a protonated guanidinium or amidinium (the H make the N+). A protonated
    carboxylic acid is left to the geometry (C=O shorter than C–OH): the test with the H on the nearer O was removed, a
    geometry that does not occur.
  - The Kekulé steps and Hückel's demand use the octet valence, so an sp2 N with three σ bonds and no stated charge
    (iminium, N-alkylpyridinium as NAD's N1, imidazolium with stated H) takes no double bond. The evidence for a cation
    would be a stated charge, a pattern, or the ring's electron count.
- **Perceive bonds across residues and units** (`_spec §1.11`): today a bond leaving a residue keeps its order (counted
  in the valence, never written), so a multiple bond between residues (a Schiff base, a crosslink) is not found, and
  `forceCompute` drops a stated one between residues of a unit. Writing them needs `assignedBonds` per unit (a decision
  must survive the next residue), the links back in the pattern cache's key with those bonds in the pattern
  (`_spec §1.5`: left out while they are never written), and a check of the partner's state before a multiple bond is
  placed (a result independent of residue order).
- **Bonds between symmetry copies** (`_spec §1.3`): core-CIF files whose `geom_bond` uses `site_symmetry` codes
  (molecules on special positions, polymeric networks) bond across copies; the model structure does not have those
  bonds. Detect and warn.
