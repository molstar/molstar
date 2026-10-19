## Graphics: faster marching cubes

- Marching cubes (molecular, Gaussian and blob surfaces, isosurfaces) reads the field at its strides
  (`Tensor.strides`, new) and each cell reads 4 corners instead of 8, with the same mesh: 220³ field, z fastest
  151 → 73 ms, x fastest 130 → 61 ms; plain fields of different layouts meshed in one session no longer slow each
  other down (519 → 61 ms); floodfilled fields (still read through `get`) 1.5–2× faster (#1987)
