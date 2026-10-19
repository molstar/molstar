## Accessible surface area

- Shrake-Rupley ASA takes a `symmetry` parameter: `'all'` (the default, unchanged) or `'per-group'`, which computes one
  unit of each symmetry group (the one nearest the structure's centre) and copies its residue areas to the other units
  of the group. Approximate where copies have other surroundings (crystal-block edges, helical segments, chains under
  different operators); when every chain is under one closed point group it differs from `'all'` only by the sphere
  points' orientation noise. Icosahedral 60-mer (350k atoms): 4.3 s → 0.1 s; 1tup crystal block of 864 copies (5M
  atoms, approximate there): 27 s → 0.4 s (#2003)
