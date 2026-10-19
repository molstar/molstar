## Graphics: faster surface recoloring

- Faster mesh color smoothing, which reruns on every color, overpaint, transparency, emissive and substance change of
  a surface, with identical output: 7qpd gaussian surface 277 → 85 ms, molecular surface 211 → 69 ms, 1tup-2x2x2
  (186k atoms, 416 surface meshes) 807 → 308 ms
- The smoothing grid spans the vertices instead of the bounding sphere: 16% fewer texels on these surfaces (−19 MB of
  grid textures on 1tup-2x2x2) with the colors at the vertices unchanged, except where the grid is coarsened to fit
  the maximum texture size, where its resolution can now stay finer
