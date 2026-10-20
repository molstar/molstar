# Reconciliation: `State.updateTree`

`updateTree(tree | builder, options)` (`state.ts`) makes the state's tree the given one and brings the cells in line.
The steps below are `update(ctx)`; an `UpdateContext` holds one pass.

## 1. The queue, transactions, undo

- Updates run one at a time (`updateQueue`, an `AsyncQueue`); a task aborted before its turn is dropped.
- `transaction(edits)`: updates inside it share one `isUpdating` span; on an error (or an update with an error inside
  it) the tree from before is restored with another update. Nested transactions rethrow to the outer one.
- `canUndo` keeps the tree from before in `history` (`historyCapacity`, 5); `undo` updates to it. An update without
  `canUndo` clears the history.
- `revertOnError` / `revertIfAborted`: a second update back to the old tree.
- `inUpdate` is true inside one update; `behaviors.isUpdating` spans a transaction.

## 2. The diff: deletes and update roots

**Deletes**: refs in the old tree, not in the new one, with a cell; in post-order (leaves first). **Update roots**: the
topmost nodes of the new tree whose cell is missing or whose version differs from their cell's
(`findUpdateRootsVisitor`): pre-order, not descending below a root, an error cell or a Null object. Two ways to find
them:

1. **Fast track**: a builder made on the current tree with one edit (`count === 1`, `lastUpdate`): no deletes, roots =
   `[lastUpdate]`, evaluated even below a Null object or an error.
2. **From the edits** (`state/diff.ts`: `findEditedDeletes`, `findEditedRoots`): a builder made on the current tree that
   moved no node and removed at most one (`editInfo.edits`, tree.md §4), while the cells are in line with the tree. The
   same deletes and roots as the whole trees, in the same order, without walking them: the deletes are the removed
   node's subtree in the old tree (post-order of a subtree is a run of the whole tree's); a root is an edited ref whose
   cell is missing or has another version, with no ancestor the whole-tree diff stops at, sorted in pre-order.
3. **Whole trees** (`findDeletes` over the old tree, `findUpdateRoots` over the new one): everything else (a builder
   that moved nodes, removed more than one, or was rebased; a `StateTree` passed in: snapshots, undo, transactions;
   after cells fell behind).

An update that throws leaves cells behind the tree (missing, or with an older transform): the next whole-tree diff finds
them; the fast track does not. So does a fast track whose builder's tree changed outside the edited node (a write to
`currentTree` directly, not counted). After either, the state takes the whole trees (`cellsInSync`) until a pass of them
completes.

## 3. Deletes

The current object moves off a deleted ref (`findNewCurrent`: a sibling, else up). Transformers' `dispose` runs in
reverse order (the top of a deleted subtree first); then each cell is unlinked (dependencies.md §3) and removed.

## 4. Cells

`initCells`:

1. Cells for every node of the roots' subtrees that has none (pre-order, `addCellsVisitor`), linked to their
   dependencies, and the cells that waited for them linked; cycles checked (dependencies.md).
2. Dependents: cells outside those subtrees that depend on a cell inside, and the cells whose dependency was deleted in
   this pass (dependencies.md §3); each becomes a root too. Their subtrees get cells as well (a dependent below a Null
   object may have children with none), and a cell that waited for one of those is a dependent too.

Then `cell.created` for the new cells, `object.removed` and `cell.removed` for the deletes. `initCellStatus`: every cell
of the roots' and dependents' subtrees takes the tree's transform, is relinked, and is set `pending`; these refs are
`ctx.scheduled`.

## 5. Evaluation

Roots, then dependents, in order, one at a time (`updateSubtree`, awaited):

- deferred if a dependency is still pending (dependencies.md §4);
- `processing`, `updateNode`, then `ok`; then the children, in order; below a Null result, the children are cleared
  instead (`doError`, silent);
- an exception: the cell is `error` (`errorText`), its object and its subtree's objects are removed (`doError`), and it
  and its subtree are no longer scheduled. The first error of the pass is the candidate current object (§6).

`updateNode`:

1. The root: nothing.
2. A parent with a Null object: a node evaluated before is disposed and gets a Null object (`updated`); a new one
   reports `created` with a Null object but its cell keeps none (open: [../../../TODO.md](../../../TODO.md)).
3. The source object: the closest ancestor whose object is of a `from` type (the parent first; any parent when `from` is
   empty).
4. Params (`resolveParams`): the definition is `params(source)`; values normalized against it once per version, then
   missing defaults filled; refs resolved (tree.md §3).
5. No cell params or a new ref: `apply`. Else `update` (no `update`: Recreate): `Unchanged` → nothing, `Updated` → in
   place, `Recreate` → `dispose` and `apply`, `Null` → `dispose`, a Null object.
6. The object gets the transform's tags.

## 6. After the pass

- Without a builder (a tree passed in), cells take the tree's state flags (`syncNewStates`, `cell.stateUpdated`).
- With a builder, cells keep theirs; the tree took the flags written since `build()` (tree.md §4). A tree passed in
  resets flags written after it was taken.
- `object.created` / `object.updated` (`in-place`, `recreate`) in evaluation order.
- The current object: the pass's candidate is the first error, else the last object created (not a ghost, not Null); it
  becomes current only with `doNotUpdateCurrent: false` (the default is `true`). Only when there is no candidate, a
  current that became Null or failed for a Null parent moves to a sibling or up (`findNewCurrent`).
- `changed`, when anything was deleted, scheduled or changed.

## 7. Events, in order

| When           | Events                                                                                     |
| -------------- | ------------------------------------------------------------------------------------------ |
| deletes found  | `behaviors.currentObject` (if the current was deleted)                                     |
| cells made     | `cell.created`, then per delete `object.removed`, `cell.removed`                           |
| scheduled      | `cell.stateUpdated` (to `pending`; a new cell is already pending)                          |
| per node       | `cell.stateUpdated` (`processing`, `ok` / `error`), `log`, `object.removed` (cleared)      |
| after the pass | `cell.stateUpdated` (state flags), `object.created` / `object.updated`, current, `changed` |

Listeners run synchronously in the middle of a pass: they may read the state and call `updateCellState`; a new
`updateTree` queues behind this one.

## 8. Cost

A builder commit with one edit costs the edited subtree; a builder commit diffed from its edits costs the edits, their
ancestors and the removed subtree (and, to sort several roots, one pass over the children of each parent where their
paths part); any other update walks both trees (O(n)), plus the evaluated subtrees, and a tree passed in syncs every
cell's state flags. With many objects in the state (hundreds of structures), the whole-tree diff on every commit made
loading O(n²) (perf-audit 14). Per evaluated node: the transformer, three status events, `params(a)`, and Immutable.js
lookups (README, Decisions).

## Decisions

- **One edit takes a fast track.** It is the common case (a UI control): the edited node is the root, with no walk of
  either tree. It differs from the whole-tree diff: it evaluates the node even below a Null object or an error, and it
  never finds cells an earlier update left behind.
- **More edits are diffed from the edits, exactly.** Presets commit a few edits at a time (a component and its
  representation, an empty component deleted), one commit per object; walking the whole state for each made loading N
  structures O(N²). The result must equal the whole-tree diff (a differential test, testing.md), so it applies only
  while the cells are in line with the tree, and the edits are recorded by the builder's tree itself, writes to
  `currentTree` included.
- **Evaluate sequentially.** Transformers are tasks that may take seconds; one at a time keeps a pass deterministic
  (events, the current object, dispose order) and lets a dependency defer its dependent (dependencies.md §4).
- **Nothing below a Null object is evaluated.** Null is "no value here" (a representation of an empty selection): its
  subtree is cleared, not errored, and the diff does not descend there.
