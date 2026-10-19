## Model: faster structures and components on large assemblies

- Child units (components, selections, focus) restrict their root unit's intra-unit bonds when those exist instead of
  computing them; the edges equal a recompute's, in the order of the root's, and bonds the root carried over from
  another frame (`remapModel`) carry over to the child. 1tup 2×2×2 (186k atoms): every second
  atom of the 416 units 78 → 1.8 ms; of the 288 units in 9 symmetry groups 2.7 → 0.1 ms; as one unit, 2000
  one-residue children 62 → 5.5 ms, every second atom 86 → 1.7 ms, the first half 111 → 2.6 ms
- [Breaking] `unit.props` is read-only; `unit.parent` is the root unit a child was made from. `remapModel` drops a
  unit's rings and resonance together with its bonds, so they are recomputed from the new bonds
