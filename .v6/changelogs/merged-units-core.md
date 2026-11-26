## Graphics: merged units and GPU-culled views

- [Breaking] The geometry param `renderMerged` is removed: units representations always merge their groups (one group
  too), at small counts measured no worse than unmerged (2 to 30 units, WebGL2 and WebGPU) and cheaper to update at
  30; the particle targets merge two or more objects of a type, or one culled on the GPU.
- Overpaint, transparency, emissive and substance layers on merged surfaces are color-smoothed where they reach (they
  were left unsmoothed on merged ones): the reached units drawn as a merged object of their own, smoothed per vertex on
  the CPU. 60 copies of 7qpd (420 units), a layer set on one unit 18.8 ms against 12.9 ms for per-unit grids, on all
  2.75 s against 2.04 s (WebGPU), the same memory (35 MB). Smoothing colors per vertex is faster: 3.67 to 2.75 s of
  that.
- On WebGPU the smoothed layers of merged surfaces are sampled per vertex by a compute kernel into the storage they are
  drawn from, as the per-unit grids sample them, with no CPU copy: 60 copies of 7qpd, a layer set on all 420 units
  2.75 to 2.13 s (per-unit grids 2.04), on one unit 18.8 to 11.3 ms. The geometry export reads them back as drawn.
- Smoothed colors and layers (overpaint, transparency, emissive, substance) of surfaces are sampled at their grid's
  cell corners, where the smoothing puts the values: they were drawn half a cell off in x and y (WebGL and WebGPU).
- The kernel sampling the smoothed layers of merged surfaces is compiled in the background once a surface smooths, and
  by `Canvas3D.precompile` for lists with merged meshes, not at the first layer.
- `Structure.instances` (the structure-instances transform) numbers every copy: with the identity among its transforms
  it gave units of one id (that copy kept its unit's, the other copies' numbering repeating it).
- On WebGPU the smoothed colors and layers of CPU meshes (Gaussian and molecular surfaces built on the CPU) are
  splatted into their grid by compute kernels, over the cells of the instances they reach: 60 copies of 7qpd, a layer
  set on all 420 units 1.9 s to 20 ms a frame per unit and to 33 ms merged, on one unit 13 ms to 2–3 ms.
- On WebGPU merged surfaces smooth their colors on the GPU too (splatted and sampled per vertex, as their layers): 60
  copies of 7qpd colored by element-index, the surfaces built in 2.7 s to 0.24 s (as without smoothing); WebGL as
  before.
