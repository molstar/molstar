## Model: faster structures and components on large assemblies

Measured on 1tup expanded to 2×2×2 cell blocks (186k atoms): as 288 units in 9 symmetry groups, or as one model of
416 units.

- Inter-unit bonds take the atoms of a unit near its partner from the unit's grid instead of transforming every atom
  for every unit pair; edges and their order are unchanged. 288 units: 84 → 44 ms, 91 → 12 MB allocated; 416 units:
  60 → 47 ms, 51 → 8 MB
- `Bundle.toStructure` and `Loci.toStructure` keep symmetry sharing: copies of a child unit share one set of props
  (lookups, bonds, boundary), as the subset builder's children do. Every second atom of the 288 units: 195 → 9 distinct
  props, child boundary 44 → 4 ms, intra-unit bonds of all child units 115 → 4 ms
- Child structures (components, selections, focus) restrict the parent's inter-unit bonds when those exist instead of
  computing them: polymer component 79 → 0.7 ms (288 units), 65 → 0.1 ms (416 units). The edges are a recompute's in
  the order of the parent's; a child shows what its parent shows, so bonds a trajectory frame carried over from an
  earlier one (`remapModel`) carry over to components created on that frame too
- Structure and loci boundary spheres above 10k points (the unit spheres' extrema, or the selected elements) use 14
  directions instead of 98; unit spheres and smaller structures and selections keep 98: structure sphere 14 → 2.6 ms,
  loci sphere of every second atom 32 → 11 ms (416 units). The spheres are ~2 % tighter and enclose the atoms as
  before; a sphere merged from several of them (the camera focusing several large loci) can miss atoms by up to ~2 Å
  (0.6 Å before)
- `Mat4.isRigid`: whether a matrix keeps distances (orthonormal linear part, affine)
- `Box3D.fromPositionArray` (interleaved xyz), `Box3D.intersect`, `Box3D.isEmpty`
