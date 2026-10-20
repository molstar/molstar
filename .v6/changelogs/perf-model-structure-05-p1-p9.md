## Model: query unions and carbohydrate links (perf audit 05)

- Structure unions (a Sequence selection to one structure: components from expressions, loci of selections, the
  targets of `within`, `intersectBy`, `exceptBy` and `isConnectedTo`; all-singleton results of `pick`, `within` and the
  PyMOL/VMD `around`) merge each unit once instead of once per input: linear in a unit's size, not quadratic. Same
  structures (units, element sets, symmetry groups). big-bcif (1tup 2×2×2 as one model: 186k atoms, 416 units):
  residues 37 → 4–6 ms, every atom as a singleton ~300 → 15–19 ms, `within(all, all, 1 Å)` 700–835 → 240–280 ms;
  singletons in one 186k-atom unit (MD and PDB files without chain ids): 10k 103 → 0.8 ms, 100k 8.5 s → 9–14 ms;
  100 overlapping halves of that unit 16 → 9 ms (#2018)
- Carbohydrate links within a chain come from one bond search per sugar ring instead of one per pair of rings: same
  links in the same order. Hexoses in one chain: 1330 sugars 312 → 6.5 ms, 5320 sugars 5.0 s → 41 ms (#2018)
