# Core: to do

The decisions still open. One entry each, with a pointer; an entry leaves when it is decided (the spec says so then).

## State (`src/state`, [\_spec](src/state/_spec/README.md))

- **Drop the Immutable workaround once Immutable is fixed.** `TransientTree` gives a map a new owner on the first write
  after deletes (`renewOwner`, [tree.md](src/state/_spec/tree.md) §2) because Immutable's `HashCollisionNode.update`
  shares the survivor's `[key, value]` entry with the source tree, and a later `set` writes into it (seen in 4.3.9,
  5.1.9, 6.0.0-beta.1; known upstream, fixed by the open
  [immutable-js#2279](https://github.com/immutable-js/immutable-js/pull/2279)). When an upgrade with that fix makes the
  canary in `src/state/_test/transient.test.ts` fail, remove `renewOwner`, the `deleted*` flags and the canary, and keep
  the other tests there.

- **Cells for every node of the tree.** The diff does not descend below Null objects and errors, so nodes added there
  have no cell; the dependents' cell pass and `DependencyWaiters` exist for that
  ([dependencies.md](src/state/_spec/dependencies.md) §3). Cells for all nodes (pending, no object) would make invariant
  1 unconditional and remove both.
- **A node evaluated below a Null parent** reports `created` with a Null object, but its cell keeps no object
  (`updateNode`, the parent-Null branch); one evaluated before gets the Null object. Pick one.
- **Dependents one level deep.** A dependent's own dependents are not scheduled in the same pass (`initCells`, "check if
  dependent cells are all proper roots").

## Data utilities (`src/data/util`)

- **A common home for sequence algorithms.** `data/util/needleman-wunsch.ts` (`globalAlign`) is here only because both
  `@molstar/model` (`align`: SEQRES mapping, superposition) and `@molstar/seq-core` (the `needleman_wunsch` aligner) use
  it, and core is the one package both depend on. The other aligners (biomsa, myers) and their limits live in
  `packages/sequences/core/src/algorithms`. Decide where sequence algorithms shared by Mol\* and Seq\* belong (a package
  or entry both may import) and move it there with its test (perf-audit 04/P1).
