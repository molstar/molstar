## Graphics: direct volume quality

- Direct volumes store their values only (were values and gradients): a quarter of the GPU memory (512³ map: 537 →
  134 MB), float data 15-20 % faster, byte data up to 10 % slower; normals computed from the filtered values
- `preIntegration` param of direct volumes (on by default): opacity integrated between consecutive samples, so narrow
  transfer peaks and steep ramps no longer show slabs, rings or noise (also at one step per cell); no measurable cost
