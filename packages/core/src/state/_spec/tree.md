# The tree

`StateTree` (`tree/immutable.ts`) is persistent; `TransientTree` (`tree/transient.ts`) edits one; `StateBuilder`
(`state/builder.ts`) is the API over a transient tree.

## 1. Structure

- `transforms`: ref → transform.
- `children`: ref → child set, in insertion order (a `ChildSet`: `size`, `has`, `first`, `values`, `forEach`, `map`,
  `toArray`). Every node has one, possibly empty.
- `dependencies`: ref → the refs whose explicit `dependsOn` names it, kept by `add`, `fromJSON` and `setDependsOn` (a
  replaced `dependsOn` drops its old edges). Explicit only: dependencies derived from params are not here
  (dependencies.md §1), so removing a node cascades along `dependsOn` alone.
- The root (`StateTransform.RootRef`) is its own parent.

Traversal: `StateTree.doPreOrder` (a visitor returning `false` skips the subtree) and `doPostOrder` (leaves first),
children in order. Serialization: `toJSON` in pre-order; `fromJSON` needs every transformer id registered
(`StateTransformer.get`).

## 2. Editing: `TransientTree`

`tree.asTransient()` makes one; `asImmutable()` returns a new `StateTree`, or the source tree when nothing changed. Each
map turns mutable once, on its first write; a child set too. `asImmutable()` freezes those maps in place (Immutable's
`asImmutable` on a mutable map), so the result also becomes the transient's new base: later writes start from new
mutable copies and never reach the snapshot, nor are they lost to a frozen map.

A map takes a new owner on the first write after deletes (`renewOwner`): Immutable keeps the survivor of two keys with
equal hashes in a node of the map's owner that shares its `[key, value]` entry with the source tree, and a later `set`
of that key writes into the shared entry, so a builder's child, transform or dependent would appear in `state.tree`
before it commits. Generated refs collide this way often. Every write to a map goes through its gate (`changeNodes`,
`changeChildren`, `changeDependencies`) and every delete through its `delete*`, which marks the map; deletes in a row
stay in place. A child set stays mutable after the root is cleared, so later children of the root are kept.

The `State` keeps its current tree as a `TransientTree` too: `updateCellState` writes state flags into it
(`assignState`), so `state.tree` changes between updates without a new version or a new identity.

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
- `TransientTree.assignState` merges state flags into a transform it already replaced since the last `asImmutable()`
  (`_stateUpdates`); no snapshot shares that object, since `asImmutable()` starts a new set (§2).

Anything else that mutates a transform breaks invariant 4 (README) and the version rule.

## 4. The builder

`state.build()` makes a `StateBuilder.Root` on `state.tree`: a transient copy, and `editInfo`:

| Field        | Meaning                                                                              |
| ------------ | ------------------------------------------------------------------------------------ |
| `sourceTree` | the tree it was made on                                                              |
| `count`      | the edits so far; with `lastUpdate`, one edit takes the fast track (reconcile.md §2) |
| `lastUpdate` | the ref of the last edit                                                             |
| `edits`      | what the builder's tree set (refs), moved and removed: the diff from edits           |
| `applied`    | set on commit: a builder commits once                                                |

`edits` is recorded by the builder's `TransientTree` (`TransientTree.Edits`), not by the builder's methods, so a write
to `currentTree` directly counts too; it is not reset by `asImmutable()`. The reconciler reads it (reconcile.md §2).

`apply` adds a node (a decorator transformer is inserted instead: between the node and its children, which move to it);
`apply` on a decorated node applies to the end of its decorator chain. `update`, `updateState`, `tag`, `dependsOn`,
`delete`, `insert`.

`getTree()`: when the state's tree is no longer `sourceTree` (another update ran since the builder was made), the
builder's actions (adds, inserts, deletes, and the params, tags, `dependsOn` and flags each edit set) are applied again
to the current tree, and `sourceTree` becomes the builder's own tree (so the commit takes the whole-tree diff).

State flags: a builder's tree is a snapshot, and flag writes (`updateCellState`) keep `state.tree`'s identity, so they
do not rebase it. The `State` logs its flag writes (`flagLog`: ref → flag → write number and value); a builder commit
writes into the new tree the flags written after `build()` (`carryFlagWrites`, which finds them in a log in write order,
not by visiting every ref ever flagged), and only those, so flags the builder set itself (`updateState`) stay. They keep
the version, so the fast track still holds. Without it, hiding a node while a builder is outstanding (a queued commit,
an abandoned builder) loses the flag from the tree, the session and undo. A plain `StateTree` passed to `updateTree`
(e.g. `getTree()` committed later) is a snapshot whose flags win (reconcile.md §6).

## Decisions

- **State flags without versions.** Hiding or collapsing a node must not evaluate it; `withState` keeps the version, and
  the flags reach cells by `updateCellState` or, without a builder, `syncNewStates`.
- **Immutable's shared entries are worked around, not patched.** The bug is in immutable 4.3.9, 5.1.9 and 6.0.0-beta.1
  (the versions checked), and a patch of the dependency would not reach the packages' consumers. The cost: the first
  write after deletes copies the paths it takes. Revisit once Immutable is fixed (immutable-js#2279, open; TODO; the
  canary in `_test/transient.test.ts` fails then).
- **A builder rebases by replaying its actions.** Edits made on an older tree are replayed on the current one rather
  than merged; an action that no longer applies throws (`Parent ... must be present in the tree`).
