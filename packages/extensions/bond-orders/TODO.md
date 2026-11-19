# Bond orders: to do

What is missing, and what was considered and deferred; an entry leaves when it is done (its `_spec` section says so
then). The spec: [\_spec](_spec/README.md).

## Next

- Check aromaticity propagation to kekulization step
- **Flag `OrderUnknown` in the bond-table readers** where the file gives no order: mol2 `un`, LAMMPS, topologies (PSF,
  PRMTOP, TOP), CIF core without `_geom_bond_valence` (core `formats/structure/`). Then `auto` can perceive those
  tables; the table fast path in `perceive.ts` (`_spec §1.9`) stays, for the tables none of whose bonds is
  `OrderUnknown` (checked once per model).

## Later

- **Benchmark**: perception over a large set (the Labute 2005 set, then the whole PDB), the data fetched on demand,
  never committed; connectivity and order comparison possibly with RDKit (RDKit-esm). A project of its own
  (`_spec §1.10`). Perception accross the PDB
- **Perceive bonds across residues and units** (`_spec §1.11`): today a bond leaving a residue keeps its order (counted
  in the valence, never written), so a multiple bond between residues (a Schiff base, a crosslink) is not found, and
  `forceCompute` drops a stated one between residues of a unit. Writing them needs `assignedBonds` per unit (a decision
  must survive the next residue), the links back in the pattern cache's key with those bonds in the pattern
  (`_spec §1.5`: left out while they are never written), and a check of the partner's state before a multiple bond is
  placed (a result independent of residue order).
- **Bonds between symmetry copies** (`_spec §1.3`): core-CIF files whose `geom_bond` uses `site_symmetry` codes
  (molecules on special positions, polymeric networks) bond across copies; the model structure does not have those
  bonds. Detect and warn.
