## Graphics

- WebGPU: small opaque spheres can be drawn by compute kernels instead of impostor quads (renderer props `sphereRaster`,
  `auto` by default, and `sphereRasterRadius`, the largest radius on screen in CSS pixels, 24 by default): a visibility buffer of
  64-bit atomics (`atomic-vec2u-min-max`) resolved in the opaque pass, for perspective views; other spheres, passes and
  devices draw impostors as before. `auto` draws the spheres of objects culled on the GPU (merged, many instances),
  where the raster was measured faster (up to 7× zoomed into large particle scenes); `on` also those of other objects.
  Objects clipped by instance or by pixel are drawn too (the caps of solid interiors as the impostors draw them, a
  clipping layer's flags per sphere), not those clipped by instance with a clipping layer
- Sphere impostors reaching behind the camera no longer draw a quad over the whole view: the quad of a sphere through
  the near plane is the bounds of its part beyond it (zoomed into 300k atoms as spacefill at 2560×1440: 68.6 → 45.6 ms
  on WebGPU; PetWorld zoomed 3.23 → 2.12 ms)
- A geometry of a single sphere is bounded by its size: its bounding radius could be exactly 0, by where it lies, culling
  it as a point
