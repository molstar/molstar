## Graphics: interactive marking and layers on large structures

- Hover, select and the overpaint, transparency, emissive and substance layers of structure representations cost what
  the loci reaches, not the number of units. 1tup 2×2×2 as one model (186k atoms, 416 units): one-chain overpaint of a
  gaussian surface 227 → 0.07 ms (merged), 212 → 1.4 ms (unmerged); marking all carbons in ball-and-stick 22 → 6.8 ms,
  every other residue in cartoon 8.5 → 0.7 ms; a residue hover 0.08 → 0.01 ms, no longer growing with the unit count.
  Marker and layer arrays are unchanged (#1995)
- Layers are color-smoothed only on the surfaces they reach, and not on members of a merged surface (`renderMerged`, by
  default above 100 units): smoothing gave members a layer type the merged object cannot show per member (#1995)
