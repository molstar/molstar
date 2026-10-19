## Accessible surface area

- Shrake-Rupley ASA takes a `symmetry` parameter: `'all'` (the default, unchanged) or `'per-group'`, which computes one
  unit of each symmetry group (the one nearest the structure's centre) and copies its residue areas to the other units
  of the group. Approximate where copies have other surroundings (crystal-block edges, open assemblies); in a closed
  point-group assembly it differs from `'all'` only by the sphere points' orientation noise. 5M-atom assembly of 864
  copies (asm-864): 27 s → 0.4 s
