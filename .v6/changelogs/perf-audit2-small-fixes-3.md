## Fixes from the perf audit

- Fix: state selections `subtree()` and `children()` skip nodes without a cell (below a Null or failed node); chained
  filters threw
- Fix: `StateSelection.findTagInSubtree` returns the first match in pre-order and stops there; it returned the last
  match after walking the whole tree
- Fix: `Bundle.toStructure` maps each copy of a unit through its own elements: a bundle on a child structure whose
  copies hold different atoms took one copy whole
- Fix: MolScript `include-surroundings` with `atom-radius` finds small units (ions, waters, ligands) within the radii's
  reach
- Fix: the out-of-plane H-bond check also runs when the only heavy neighbour is the first atom of its unit
- Fix: water bridges stay within a model; an ensemble shown as all models bridged one model's atoms through another's
  water
- Fix: impostor double, triple and aromatic bonds follow a size theme change; their offsets kept the old size
- Fix: CPU color smoothing (no graphics context, headless capture) colors vertices on grid planes; they were black
- Fix: molecular surfaces and CPU gaussian densities yield and can be aborted for small radii (a size factor below
  ~0.55 ran the whole structure in one chunk)
- Fix: a nested MolScript `pick` keeps the outer query's current structure
- Fix (Prism): a whole-chain selection or highlight (past the size gates) takes a ligand's hydrogens
- Performance: particle-hierarchy and particle-entity colors index the particles once per list. 10M particles, a theme
  update: 0.97–1.12 s → 11 ms and 0.23–0.35 s → 0.1 ms
- Performance: the focus representation updates its selections without a mutative walk of the new params and its
  console warning per update
- Performance: the unit molecular surface uses arrays of the unit's size; they were sized to the model per unit.
  1.12M atoms in 2,496 units: position prep 122 → 38 ms, 8.1 GB → 31 MB allocated per build
