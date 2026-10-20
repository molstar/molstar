## Model: formats and alignment (perf audit 04)

- One Needleman–Wunsch for Mol\* and Seq\*: `@molstar/core/data/util/needleman-wunsch` (`globalAlign`, affine gaps
  over integer codes, two rows of scores and a byte of traceback a cell). `@molstar/model`'s `align` (PDB SEQRES
  mapping, sequence superposition) and Seq\*'s `needleman_wunsch` aligner run on it with the same alignments, ties and
  scores as before; the two differed in leading end gaps and gap ties, now options of the shared DP. 3000 residues:
  `align` 210 ms / 1 GB peak → 48 ms / 32 MB, `needleman_wunsch` 82 → 45 ms; 10000 residues: `align` 3.1 s / 3 GB →
  0.5 s / 97 MB. PDB SEQRES: chains with the same SEQRES and observed residues align once (#2017)
- Streamlines compute volume gradients on demand: nothing is cached on the `Grid` (12 bytes a voxel before;
  `Grid.getGradients` no longer caches either). 256³: 201 MB → 0 retained, 532 → 139 ms at the default seed density
  (+15 % at the maximum density) (#2017)
- `Grid.getIsosurfaceBoundingSphere` (CPU isosurface, wireframe and dots) scans in memory order and caches the last iso
  value: 256³ 33–122 → 1–9 ms (#2017)
- Multi-model mmCIF/BCIF: an in-order model reads its typed `atom_site` columns in place instead of copying them.
  multi50.bcif (7qpd × 50): 131 → 120 ms, 35 MB less allocated, 66 → 29 MB retained (#2017)
- `guessElementSymbolString` (GRO, PSF, TOP, MOL2, CIF-core) without a regex: 2M atoms 83 → 31 ms (#2017)
- `IndexPairBonds` edge props in typed arrays: −28 MB per 1M bonds (#2017)
- struct_conn partners are looked up in their asym's entity only: 10k rows over 100 entities 38 → 5 ms (#2017)
- IHM: model group names from one pass (5000 models 66 → 0.7 ms); the residue source index skips its sort when residues
  are in order (#2017)
