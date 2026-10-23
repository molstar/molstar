# Model: to do

The decisions still open, and defects found in the code. One entry each, with a pointer; an entry leaves when it is
decided or fixed. Spec: [\_spec](_spec/README.md). Performance items and bugs of `model/structure` are tracked in
`.v6/perf-audit/05-model-structure.md` (P1, P7–P12; B1–B8).

## Decisions open

- **Rings of a child from its parent.** The parent's rings whose atoms are all in the child are likely exact for
  children of whole residues; not shown for partial ones ([unit.md](_spec/unit.md) §4).
- **The parent of a remapped unit.** `remapModel` drops it; once `Structure.remapModel` receives the remapped root (B1),
  pass the root's unit with the same id ([unit.md](_spec/unit.md) §3).
- **Per-operator index-pair bonds and shared props.** Symmetry copies share `props`, so index-pair bonds listed per
  operator (`molstar_bond_site` with operator ids) are computed by whichever copy asks first and then serve all of them,
  though [unit.md](_spec/unit.md) §1 and §4 treat them as depending on the operator. No reader sets operator ids today
  (05-model-structure P8); key the bonds by operator, or keep such copies from sharing props, once one does.
- **`model/` imports `props/`** (`model.ts`, `structure.ts`, `properties.ts`, mmCIF export, SIFTS superposition) while
  `props/` imports `model/`: one layer or two?
- **Format properties set by readers are not registered.** `FormatPropertyProvider.set` writes the bag but not
  `customProperties`, so `cifExport` skips values a reader set (CCD `ComponentBond`); the volume provider registers.
- **`focus-first-residue.ts` reads `unit.props.polymerElements`** (plugin core), so it sees only units whose polymer
  elements happen to be computed; `unit.polymerElements` was likely meant.

## Defects

- `Structure.instances` gives the copies new ids from 0 and keeps the original unit, with its id, for an identity
  transform: ids collide (1cbs with `[I, T]`: `0 0 1 1 2 2`), and without an identity the copies' ids differ from their
  parent structure's, which children and loci assume (README invariant 3). Reached from the plugin (`hierarchy.ts`).
- `Model.getInvertedAtomSourceIndex` (`invertIndex`) sizes its array by the model's atom count but indexes it by
  `atomSourceIndex`, which holds absolute rows of the file's `atom_site` (`sortAtomSite`): for the second and later
  models of a file whose hierarchies differ, the writes fall past the end and the inverted index is zeros. Consumers:
  index-pair inter-unit bonds (`inter-compute.ts`). Which index space `IndexPairBonds` uses (source rows of the file, of
  the model, or element indices for `molstar_bond_site`) has to be settled first.
