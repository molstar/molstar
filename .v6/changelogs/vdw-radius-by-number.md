## Model: vdW radii

- `VdwRadiusOfNumber`: the vdW radius by atomic number (`AtomicHierarchy.derived.atom.atomicNumber`), a table lookup;
  the physical size theme, `Model.getAtomicRadii` and units' radii use it instead of a lookup of the element symbol per
  atom (PetWorld's targets, 3.7 M atoms: the size lookup 157 → about 65 ms).
- `VdwRadius` gives the default radius for the elements without one in its table (110 to 118, Ds to Og), not
  `undefined`.
