## Performance: small allocations and repeated work (from #1932, Taylor Hoffmann)

- CIF: a database with aliases (CIF core, RELION STAR) flattens its frame once instead of once per table (12 tables
  for CIF core)
- Chain-id, entity-id and operator-hkl color themes cache their serial maps on the structure: the themes of every
  representation of a structure share them, and operator-hkl's params getter no longer recomputes the map
- `StructureElement.Stats` and `StructureElement.Loci.extendToRadius` reuse their maps between calls
