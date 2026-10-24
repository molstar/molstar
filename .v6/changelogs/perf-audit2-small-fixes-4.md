## Fixes from the perf audit

- Fix: `Structure.entityIndices` of coarse units reads their own chains: a coarse-only structure (IHM) had no entities,
  and its mmCIF export wrote empty `_entity` / `_entity_poly`
- Fix: coarse sequences made from sphere ranges (no `_entity_poly_seq`) start at the first seq id and keep the chain's
  last element; they were off by one
- Fix: MOL2 bonds name their atoms by id, so files whose atom ids are not 1..N in file order bond the right atoms; a
  bond to an id no atom has is dropped
- Fix: a segmentation volume with a voxel label in no set loads; it threw
- Fix: accessible surface area reaches the neighbours of atoms with radii above the standard ones (Se, metal ions);
  areas next to them were overestimated
- Fix: custom shapes with group ids at or above their item count map their groups (focus and loci bounds)
- Fix: isosurface wireframes of more than 16.7M vertices join the right vertices
- Fix: views, windows and array copies of a masked BinaryCIF number column read missing rows as 0 (an encoder that
  does not zero masked slots gave garbage `label_seq_id` or formal charges in sorted models)
- Fix: segment highlight and selection over a whole segmentation (slice, direct volume) mark the instance the loci
  names and each element's own segments
- Fix: the plugin log keeps its last 1000 entries; it grew by every transform's timing line for the whole session
- Fix (Prism): a trajectory or model-index step with the sequence linked keeps the 3D hover and records nothing; each
  step cleared the hover and recorded a `select.set`
- Performance: interaction features are reused when only contact params change (a distance slider). 186k atoms: a
  recompute 284–299 → 236–244 ms, and no ~4 MB kept per change
- Performance: gzip files decode with the platform's gzip decoder, without a CRC-32 pass in JS. 100 MB of mmCIF:
  334–346 → 129–135 ms
- Performance: the volume-segment color theme looks segments up in a map. 50k segments, a recolor: 470 → 5 ms
