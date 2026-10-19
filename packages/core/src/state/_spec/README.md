# Mol\* state: the specification

`@molstar/core/state` (`packages/core/src/state/`): the state tree and its reconciler. A `State` holds a persistent tree
of transforms (what to compute) and a map of cells (what was computed from them); `updateTree` brings the cells in line
with a new tree. A plugin has two: `plugin.state.data` and `plugin.state.behaviors`. This spec is for working on the
code without reading all of it: the model, the invariants the code relies on, and the decisions (each file's
"Decisions": the rule, then why). Code comments cite sections (`reconcile.md §2`). The code is the truth for how. Open
items: [../../../TODO.md](../../../TODO.md) (State).

| File                               | What it covers                                                                      |
| ---------------------------------- | ----------------------------------------------------------------------------------- |
| [tree.md](tree.md)                 | `StateTree`, `TransientTree`, versions, what is mutated in place, the builder       |
| [reconcile.md](reconcile.md)       | `updateTree`: the queue, the three diffs, cells, evaluation, events, current object |
| [dependencies.md](dependencies.md) | dependencies between cells: where they come from, links, waiters, deferral, cycles  |
| [testing.md](testing.md)           | the tests, the differential test, how to test a change, measuring                   |

## The model

| Piece       | Where                                    | What it is                                                                                                                                        |
| ----------- | ---------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| transform   | `transform.ts`                           | a node: `ref`, `parent`, `transformer`, `params`, `state` flags, `tags`, `dependsOn`, `version`; a value (tree.md §3)                             |
| transformer | `transformer.ts`                         | how a node computes: `from`/`to` types, `params(a)`, `apply`, `update`, `dispose`, `getDependencies`, `isDecorator`                               |
| tree        | `tree/immutable.ts`, `tree/transient.ts` | transforms by ref, ordered children, explicit dependencies; persistent, edited through a `TransientTree` (tree.md)                                |
| cell        | `object.ts` (`StateObjectCell`)          | per ref: the transform it was last evaluated with, `status`, `obj`, `params` (definition and values), `cache`, `state` (a copy), dependency links |
| object      | `object.ts` (`StateObject`)              | a typed result; `StateObject.Null`: no value, nothing below it is evaluated                                                                       |
| builder     | `state/builder.ts`                       | edits on a transient copy of a tree, committed through `updateTree` (tree.md §4)                                                                  |
| selection   | `state/selection.ts`                     | queries over cells (`state.select`, `StateSelection.Generators`)                                                                                  |
| spine       | `tree/spine.ts`                          | a transformer's view of the ancestors of the node being evaluated                                                                                 |
| actions     | `action.ts`, `action/manager.ts`         | tasks applied to a cell (`state.applyAction`), registered per state                                                                               |

Statuses: `pending` (scheduled, not evaluated yet), `processing`, `ok`, `error` (`errorText`).

## Invariants

The reconciler's shortcuts rely on these; a change that breaks one must change the shortcut too.

1. A cell exists for every ref of the tree, except below an error cell or a Null object, where the diff does not descend
   (reconcile.md §2; nodes added there get no cell until their ancestor changes). After an update completes, every
   cell's ref is in the tree.
2. A cell's `transform` has the version of the tree's transform with its ref, except below an error cell or a Null
   object, and after an update that threw (the cells it did not reach; the next whole-tree diff finds them, reconcile.md
   §2).
3. A child set is a new object exactly when the children changed (structural sharing): Prism's change signatures add
   child sets (`prism/viewer/src/sources/scene.ts`).
4. A transform is shared by every tree and snapshot that has it: never mutate one, except what tree.md §3 lists.
5. Dependency links (`dependencies.dependsOn`, `dependentBy`) hold live cells only; a dependency with no cell is a wait
   (dependencies.md §3).

## Decisions

- **A persistent tree beside mutable cells.** Trees are values: builders, undo history, transactions and snapshots keep
  old ones for free and compare by identity; cells hold the expensive results, updated in place. The diff between two
  trees decides what to evaluate.
- **Versions decide, not params.** A node is evaluated again when its version differs from its cell's; the reconciler
  never compares params. Cheap and exact, provided a change of params comes with a new version (tree.md §3).
- **Immutable.js is the persistent map.** It keeps a commit O(edits · log n) and a snapshot free; its lookups (string
  hashing in JS) are most of the cost of walking a whole tree.
- **Two states in a plugin, one reconciler.** The behaviors' state is a `State` like the data's: the same rules apply.
