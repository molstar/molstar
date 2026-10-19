## Volumes: field reads at the strides

- Marching cubes (molecular, Gaussian and blob surfaces, isosurfaces) reads the field at its strides
  (`Tensor.strides`, new) and each cell reads 4 corners instead of 8, with the same mesh: 220³ field, z fastest
  151 → 73 ms, x fastest 130 → 61 ms; plain fields of different layouts meshed in one session no longer slow each
  other down (519 → 61 ms); floodfilled fields (still read through `get`) 1.5–2× faster (#1987)
- Volume gradients (`Grid.getGradients`, streamlines) read at the strides too: 160³ grid 39 → 27 ms, and 143 → 18 ms when
  a second layout was seen in the session; same gradients and streamlines (#1987)
- Fix the vertex groups of volume isosurfaces with floodfill `inside`: cell ids below the iso value inside the
  surface became twice the iso value (#1987)
