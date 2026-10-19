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
  computing them: polymer component 79 → 0.7 ms (288 units), 65 → 0.1 ms (416 units). The edges are a recompute's on
  every pinned fixture, in the order of the parent's. A child shows what its parent shows: where a recompute could
  differ (a struct_conn partner outside the child, water and ion pairing decided by a unit's first element) and for
  bonds a trajectory frame carried over from an earlier one (`remapModel`)
- Structure and loci boundary spheres pick their EPOS directions by size: 98 up to 10k points (the unit spheres'
  extrema, or the selected elements), 14 above; before, structures used 14 above 500 units and loci always 98. Unit
  spheres keep 98. 416 units: structure sphere 14 → 2.6 ms, loci sphere of every second atom 32 → 11 ms. Radii move
  by a few % either way (−4 % to +1 % on the fixtures) and enclose the atoms as before; a sphere merged from several
  14-direction ones (the camera focusing several large loci) can miss atoms by up to ~2 Å (0.6 Å before)
- `Mat4.isRigid`: whether a matrix keeps distances (orthonormal linear part, affine)
- `Box3D.fromPositionArray` (interleaved xyz), `Box3D.intersect`, `Box3D.isEmpty`
- [Breaking] `Cage.vertices` and the fan and strip vertices of `MeshBuilder` are `Float32Array`, see the
  [API change ledger](../plans/breaking-v6-changes.md)
- `Vec3.fromFloat32Array` and `Vec3.toFloat32Array`: `fromArray` / `toArray` restricted to Float32Array, for hot loops
  (the generic ones are megamorphic app-wide: 2x slower per element); mesh fans, strips and cages use them
