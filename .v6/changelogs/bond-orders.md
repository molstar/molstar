## Bond providers

- `@molstar/bond-orders-extension`: bond-order perception from 3D coordinates (Sayle 2001, Labute 2005, Zhang 2012) as
  a bond provider, `'bond-order-perception'` with modes `auto` (orders the file does not state, outside canonical
  residues) and `forceCompute` (every bond); `BondOrders` is its registry entry. `examples/bond-orders`: a Prism viewer
  with it
- `@molstar/model`: `computeResidueRings(unit, start, end, bonds?)`, the rings of one residue as `computeRings` gives
  them, for a bond provider that needs a few residues' rings, not the whole unit's
