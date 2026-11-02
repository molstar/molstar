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
- **Bonds between symmetry copies** (`_spec §1.3`): core-CIF files whose `geom_bond` uses `site_symmetry` codes
  (molecules on special positions, polymeric networks) bond across copies; the model structure does not have those
  bonds. Detect and warn.
