## Fixes from the perf audit

- Fix: volume and particle representations keep their state (visibility, alpha, overpaint, transparency, clipping,
  transform) on a new render object: each streamed volume box, new data that is not equivalent, a recreated visual
- Fix: segmentation surfaces touching a face of the grid close at the face; reads past it wrapped into the next row and
  added geometry outside the grid
- Fix: ellipsoids (anisotropic ADPs) for every atom: the last atom and atom ids that do not start near 0 had none, and
  model 2 onwards showed model 1's (where atom ids are unique over the models)
- Fix: a palette set with offsets (volume-value) has a color everywhere; integer offsets and the part after the last
  stop were black
- Fix: the external-structure color theme colors from the nearest atom, and blob surfaces assign atoms to the nearest
  centroid (`approxNearest` gave the nearest last while callers read the first)
- Fix: ionic interactions no longer take a ligand's nitrogens (an amide N) for negative charges
- Fix: PDB `CONECT` bonds listed on the higher-serial atom only are read; a triple bond is written as `trip` and reads
  back as a triple bond
- Fix: resolving an asset that fails to read (an HTML error page as JSON, a corrupt file) no longer keeps it held
- Fix: `applyOrUpdateTagged` with new params and extra tags sets the tags too
- Fix: the picking granularity select lists residue, chain, entity, model and operator as plain options, not under
  'With Symmetry'
- Fix: mirrored aromatic bonds as dashed mesh cylinders keep their top caps with an odd dash count and no stub
- Performance: impostor cylinders (ball-and-stick, backbone, interactions) are written into the buffers directly.
  Ball-and-stick of 1.12M atoms: 2.32–2.45 → 2.14–2.25 s, 2.0 → 1.2 GB allocated
- Performance: the molecular surface's torus pass runs on scalars (the same field). 186k atoms: 6.42 → 6.10 s
