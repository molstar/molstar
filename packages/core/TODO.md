# Core: to do

The decisions still open. One entry each, with a pointer; an entry leaves when it is decided (the spec says so then).

## State (`src/state`, [\_spec](src/state/_spec/README.md))

- **Cells for every node of the tree.** The diff does not descend below Null objects and errors, so nodes added there
  have no cell; the dependents' cell pass and `DependencyWaiters` exist for that
  ([dependencies.md](src/state/_spec/dependencies.md) §3). Cells for all nodes (pending, no object) would make invariant
  1 unconditional and remove both.
- **A node evaluated below a Null parent** reports `created` with a Null object, but its cell keeps no object
  (`updateNode`, the parent-Null branch); one evaluated before gets the Null object. Pick one.
- **Snapshots share state flags.** `TransientTree.assignState` merges flags in place into a transform it already
  replaced (`_stateUpdates`), and undo, transaction and builder snapshots taken since share that object: they see later
  `updateCellState` changes.
- **Dependents one level deep.** A dependent's own dependents are not scheduled in the same pass (`initCells`, "check if
  dependent cells are all proper roots").
