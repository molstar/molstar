# The tree

`StateTree` (`tree/immutable.ts`) is persistent; `TransientTree` (`tree/transient.ts`) edits one; `StateBuilder`
(`state/builder.ts`) is the API over a transient tree.

## 1. Structure

- `transforms`: ref → transform.
- `children`: ref → child set, in insertion order (a `ChildSet`: `size`, `has`, `first`, `values`, `forEach`, `map`,
  `toArray`). Every node has one, possibly empty.
- `dependencies`: ref → the refs whose explicit `dependsOn` names it. Explicit only: dependencies derived from params
  are not here (dependencies.md §1), so removing a node cascades along `dependsOn` alone.
- The root (`StateTransform.RootRef`) is its own parent.

Traversal: `StateTree.doPreOrder` (a visitor returning `false` skips the subtree) and `doPostOrder` (leaves first),
children in order. Serialization: `toJSON` in pre-order; `fromJSON` needs every transformer id registered
(`StateTransformer.get`).

## 2. Editing: `TransientTree`

`tree.asTransient()` makes one; `asImmutable()` returns a new `StateTree`, or the source tree when nothing changed. Each
map turns mutable once, on its first write; a child set too.

The `State` keeps its current tree as a `TransientTree` too: `updateCellState` writes state flags into it
(`assignState`), so `state.tree` changes between updates without a new version.

## 3. Versions, and what is mutated in place

A transform is a value. `withParams`, `withTags`, `withDependsOn` and `withParent` make a new object with a new version;
`withState` keeps the version (state flags never trigger an evaluation). Versions are UUIDs (`UUID.create22`), or a hash
of the params with `getTree({ useHashVersion: true })` and in MVS (`setParamsHashVersion`), so equal params give equal
versions across snapshots.

Mutated in place, on objects trees share:

- `resolveParams` (reconcile.md §5) replaces `transform.params` with the normalized values once per version
  (`_normalized_param_version`), later fills missing defaults into them, and sets `getValue` on `ValueRef`/`DataRef`
  values.
- `setParamHashVersion` writes `version` on the transforms a builder touched.
- `TransientTree.assignState` merges state flags into a transform it already replaced in this transient
  (`_stateUpdates`); snapshots taken before share that object (an open item: [../../../TODO.md](../../../TODO.md)).

Anything else that mutates a transform breaks invariant 4 (README) and the version rule.

## 4. The builder

`state.build()` makes a `StateBuilder.Root` on `state.tree`: a transient copy, and `editInfo`:

| Field        | Meaning                                                                              |
| ------------ | ------------------------------------------------------------------------------------ |
| `sourceTree` | the tree it was made on                                                              |
| `count`      | the edits so far; with `lastUpdate`, one edit takes the fast track (reconcile.md §2) |
| `lastUpdate` | the ref of the last edit                                                             |
| `applied`    | set on commit: a builder commits once                                                |

`apply` adds a node (a decorator transformer is inserted instead: between the node and its children, which move to it);
`apply` on a decorated node applies to the end of its decorator chain. `update`, `updateState`, `tag`, `dependsOn`,
`delete`, `insert`.

`getTree()`: when the state's tree is no longer `sourceTree` (another update ran since the builder was made), the
builder's actions are applied again to the current tree, and `sourceTree` becomes the builder's own tree (so the commit
takes the whole-tree diff).

## Decisions

- **State flags without versions.** Hiding or collapsing a node must not evaluate it; `withState` keeps the version, and
  the flags reach cells by `updateCellState` or, without a builder, `syncNewStates`.
- **A builder rebases by replaying its actions.** Edits made on an older tree are replayed on the current one rather
  than merged; an action that no longer applies throws (`Parent ... must be present in the tree`).
