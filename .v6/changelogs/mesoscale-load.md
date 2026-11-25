## Graphics and model: faster mesoscale loads

- Loading a mesoscale particle file in the Viewer (WebGPU; parse, build, first full frame): PetWorld
  (`postsynstubby_2`, 517k particles over 3.7 M target atoms) 3.0–3.7 s → about 2.0 s, MG_1184 cellpack (2.1 M
  particles, 2.9 M target atoms) 5.1–5.7 s → about 3.1 s.
- BinaryCIF models of millions of atoms 1.6–1.75× faster to build: the element symbols kept a dictionary (not a string
  per atom), the atomic numbers once per element, a residue's trace and direction atoms found by their name's
  dictionary index.
- mmCIF particle lists 1.4–2.1× faster: the operators read as whole columns, operator expressions parsed to operator
  indices (a range is not a string per id), a chain's lookups once per chain, the atom rows read once.
- The bounding sphere of 1000 or more instances is the sphere around their box (center and largest distance), for
  every renderable: as tight or tighter than before (a crystal lattice of 3PQR 8% smaller, the 1M4X capsid 1%), many
  times cheaper (MG_1184: 616 → about 60 ms).
- Instance transforms of millions of instances 2.3× faster to set up (the reflection test and the instance grid
  without calls per instance, no copies before the grid; `createTransform` takes `adopt` to keep a new array).
- GPU culling's sync at the commit about 2× faster (934 → 455 ms on MG_1184): the instance copy reads the instances
  from the identity region it writes, not from a packed copy of them (80 bytes an instance less to upload and hold).
- Structure targets' spheres: the serial index of an element without a search, the sizes written unit by unit, the
  level of detail order in one pass.
