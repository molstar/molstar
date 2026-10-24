## Fixes from the perf audit

- Fix the symmetry mates of a structure whose units carry an operator (a moved structure, an assembly): the centre was
  tested in the wrong frame
- Fix `Loci.toExpression` merging the elements of different models of a trajectory structure (shared operator names)
- Fix a hang of the loci stats on units with the `Partitioned` trait
- Fix `StructureSubsetBuilder.getStructureDeduplicate` giving every unit the last unit's elements
- Fix MolScript dropping the first object-form argument of `ring-fingerprint`
- Fix the DSSP backbone H-bond graph storing its energies in input order instead of edge-slot order
- Performance: volume readers (CCP4, DX, Cube, DSN6) compute min, max, mean and rms in one pass. 256³: 47 → 20 ms
- Performance: segmentation volumes bound their sets in memory order (160³: 61 → 30 ms); the CCP4 writer reads rows
  by strides (256³: 41 → 27 ms, x slowest 381 → 174 ms)
- Performance: the PDB unobserved-residue scan keys residues, not atoms (100k atoms: 8.1 → 4.1 ms); kinemage points
  share their master tag arrays (232 → 124 B per point)
- Performance: inter-unit distance bonds are searched at the largest pairing threshold (2.95 Å), not 4 Å. 1.1M atoms:
  136 → 107 ms
- Performance: the per-model cache of intra-unit bonds by element set keeps the 64 sets used last (2000 sets: 98 →
  3 MB retained)
- Performance: H-bond geometry checks allocate nothing per candidate (186k atoms: 354 → 286 MB collected during the
  computation)
