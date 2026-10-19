## Graphics: faster surface recoloring

- Overpaint, transparency, emissive and substance on an instanced mesh surface (assemblies, symmetry mates) smooth
  only the instances the layer reaches and their neighbours, with identical output. One chain on 1tup ×864 (5M atoms,
  `low` quality): gaussian surface overpaint 0.88 s → 38 ms, transparency 0.66 s → 50 ms, molecular surface 0.81 s →
  63 ms; the smoothing's accumulators 2.2 GB → 1.4 MB (peak RSS 1.8 → 0.3 GB). ×192: 0.23 s → 23 ms; ×32 (186k atoms):
  114 → 23 ms. The smoothing grid still spans every instance
