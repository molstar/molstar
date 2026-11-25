## Graphics: levels of detail by size on screen

- Levels of detail of spheres and points by their size on screen (`lodPixelSize`, the smallest sphere drawn in CSS
  pixels, 0 for none): the same look at any resolution, zoom and field of view, from a chain of strides 1, 10, 100,
  1000. The graphics modes (`Spheres.LodPixelPresets`: `performance` 6, `balanced` 4, `quality` 3, `ultra` 1.5) set it
  in the mesoscale preset and the Mesoscale Explorer; particles default to 3. PetWorld zoomed at 720p (WebGPU, opaque pass): 20.5 → 7.2 ms.
- [Breaking] The levels by distance are removed: the `lodLevels` param of spheres, points, spacefill, volume dots and
  particles, `Spheres.LodLevelsPresets` and `getLodLevelsParam`; `Spheres.LodPixelPresets` are pixel sizes.
