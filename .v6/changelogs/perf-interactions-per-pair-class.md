## Model: interactions per symmetry-equivalent unit pair

- Interactions param `symmetry` (named as ASA's): 'all' (default, as before) or 'per-pair-class', which searches
  inter-unit contacts once per class of symmetry-equivalent unit pairs (same symmetry groups and relative operator) and
  checks only line of sight per pair. Units bonded to other units, with operator-scoped bonds, non-rigid operators or
  atom geometries that differ from their group's are searched in full. Equal to 'all' except a contact on a threshold
  within the operators' precision (~1e-4 Å): no difference on 1tup assemblies of 32–864 copies or a closed C6. 5M-atom
  assembly (864 copies): interactions 10.7 → 3.3 s; 1.1M atoms 2.3 → 0.76 s
- `computeInteractions(…, { stats })` counts the close unit pairs, classes and replayed pairs
