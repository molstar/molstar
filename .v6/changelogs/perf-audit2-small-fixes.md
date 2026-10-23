## Fixes from the perf audit

- Fix: gzip files with an extra field (`FLG.FEXTRA`: single BGZF blocks, some archivers) open; the field's 2-byte
  length was not skipped
- Fix: the interactions extension no longer throws with weak hydrogen bonds on (the weak-H-bond refiner read the
  intra-unit contacts it skips)
- Fix: accessible surface area, its color theme and `is-buried`/`is-accessible` work on structures with coarse units
  (IHM, 8zzh): coarse elements have no value; ASA is offered only for structures with atomic units
- Fix: the label representation works on coarse units: coarse chains are labelled by their asym id, residue and atom
  labels skip coarse units
- Fix: the static components of a moved or superposed structure stay where the structure is during trajectory playback
  (a whole unit took the untransformed root's unit on every frame)
- Fix: visuals made after a representation was moved (a new group, or a GPU/CPU or backend switch) get its transform;
  they stayed at the identity
- Fix: volume-value and other trilinear reads of a periodic map (MTZ, structure factors, a whole-cell CCP4) read the
  periodic image for positions outside the cell, not the wrong voxel or NaN
- Fix: volume-value with `isRelative` colors by the value relative to mean and sigma over the whole palette; a custom
  domain is in sigma units too; direct volume rendering and Prism's transfer axis place the domain right
- Fix: frames of a coordinate trajectory (XTC, DCD, ... over an SDF, MOL, MOL2, CIF core or LAMMPS model) keep the
  model's explicit bonds and bond orders
- Fix: the entity-source color theme works when a non-polymer entity comes before a polymer one, or a source row names
  a ligand
- Performance: the interactions size slider updates a uniform instead of rebuilding the impostor geometry, and the
  first update after building interactions no longer rebuilds them. 1tup ×32: size step 44–46 → 2–3.5 ms, first update
  51–52 → 5 ms
- Performance: structure boundaries no longer copy the extrema of every moved unit. asm-192 5.5–5.7 → 4.2 ms and
  33 → 1.5 MB allocated; asm-864 29 → 19 ms and 144 → 0.3 MB
- Performance: `Column.view` (sorted `atom_site`) gathers each column kind in its own loop
  (`arrayPickIndices`). Model stage 1.1M atoms: 411–417 → 382–387 ms, 437 → 373 MB allocated
- Prism: showing a frame (`trajectory.frame`, `particles.frame`) logs no transform timing lines; the lines of each seek
  pushed every other line out of the Log's Transforms
