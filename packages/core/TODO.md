# Core: to do

The decisions still open. One entry each, with a pointer; an entry leaves when it is decided (the spec says so then).

## State (`src/state`, [\_spec](src/state/_spec/README.md))

- **Cells for every node of the tree.** The diff does not descend below Null objects and errors, so nodes added there
  have no cell; the dependents' cell pass and `DependencyWaiters` exist for that
  ([dependencies.md](src/state/_spec/dependencies.md) §3). Cells for all nodes (pending, no object) would make invariant
  1 unconditional and remove both.
- **A node evaluated below a Null parent** reports `created` with a Null object, but its cell keeps no object
  (`updateNode`, the parent-Null branch); one evaluated before gets the Null object. Pick one.
- **A builder's removal can leak into the state's tree.** The first random generator of `_test/diff.test.ts` with
  generated refs instead of explicit ones (which changes the nodes it picks) hit it at seed 5 of 1000 on the code before
  the diff from edits: a child added by a builder appears in a child set of `state.tree` before that builder commits,
  and the next whole-tree diff fails (`findDeletes`: a child without a transform). Only with removals; likely a mutable
  child set shared between a `TransientTree` and a frozen tree. A bug to fix, with that generator as the test.
- **Dependents one level deep.** A dependent's own dependents are not scheduled in the same pass (`initCells`, "check if
  dependent cells are all proper roots").
