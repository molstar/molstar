## Graphics: faster surface recoloring

- Overpaint, transparency, emissive and substance on an instanced mesh surface (assemblies, symmetry mates) smooth
  only the instances the layer reaches and their neighbours, with identical output. One chain on 1tup ×864 (5M atoms,
  `low` quality): gaussian surface overpaint 794 → 42 ms, transparency 604 → 54 ms, molecular surface 926 → 50 ms; the
  smoothing's accumulators 2.2 GB → 1.4 MB (peak RSS 1.8 → 0.34 GB). ×192: 168 → 16 ms; ×32 (186k atoms): 45 →
  12 ms. The smoothing grid still spans every instance.
