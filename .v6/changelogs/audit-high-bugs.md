## Fixes from the perf audit

- Fix: a visibility (or other state) flag set while a state builder is outstanding (made but not yet committed, or a
  commit waiting in the queue) is no longer lost from the state tree: saved sessions, undo and `setSnapshot` keep it.
  A builder can also be edited again after `getTree()` (it threw `Could not find node`)
- Fix: models 2..N of a multi-model mmCIF or PDB file have their assemblies, spacegroup and symmetry mates again, and
  explicit bonds from `_molstar_bond_site` (the dynamic format properties of every model after the first one asked were
  `undefined`); each model reads only its own rows of `_molstar_bond_site`
- Fix: attaching a dynamic custom model property (e.g. the RCSB validation report) to one model of a multi-model file
  no longer disposes the assets of the same property on the other models
