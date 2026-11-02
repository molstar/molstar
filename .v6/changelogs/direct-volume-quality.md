## Graphics: direct volume quality

- Direct volumes store their values only (were values and gradients): a quarter of the GPU memory (512³ map: 537 →
  134 MB); normals by central differences of the filtered values, with the grid's transform: lit alike from every side,
  also mirrored grids
- Direct volumes follow a backend switch with the kind of grid the new context has (2D on WebGL1, 3D elsewhere)
- `preIntegration` param of direct volumes (on by default): opacity integrated between consecutive samples, so narrow
  transfer peaks and steep ramps no longer show slabs, rings or noise (also at one step per cell); no measurable cost
