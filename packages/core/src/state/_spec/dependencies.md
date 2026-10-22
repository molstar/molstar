# Dependencies between cells

A transform may depend on nodes other than its ancestors (sibling-like: never on an ancestor or a descendant). Its
transformer gets their objects as `dependencies` (by ref) in `apply` and `update`, and it is evaluated again when one of
them is.

## 1. Where they come from

`StateTransform.getEffectiveDependsOn` (`transform.ts`), de-duplicated, without the node itself and the root:

1. explicit `dependsOn` (`apply(..., { dependsOn })`, `StateBuilder.To.dependsOn`);
2. `transformer.definition.getDependencies(params)`;
3. the refs of `PD.ValueRef` / `PD.DataRef` values in the params (`PD.collectRefs` over the transformer's default
   schema; without one, a structural scan for `{ ref, getValue }` values). The default schema is
   `params(undefined, undefined)`, built once per transformer and kept (`tryGetDefaultParamDefinition`): without a
   source or a context it must not read anything that changes (a registry), and a transformer's must not depend on when
   it is first asked. The transforms of a snapshot's frames are new objects on every frame, so a cache on the transform
   did not hold.

Only (1) is in the tree (`tree.dependencies`): removing a node removes the nodes that name it in `dependsOn`, not those
that reach it through params.

## 2. Links

`cell.dependencies.dependsOn` and `dependentBy`: cells, both ways. A new cell is linked when cells are made
(reconcile.md §4); a scheduled cell is relinked to its transform's current set (`relinkCells`). A new dependency on a
ref that is not in the tree throws (`Cannot depend on a non-existent transform.`), and the update with it; a deleted
dependency the cell already waits for does not (§3). Cycles are checked whenever links are added (`StateTreeCycleError`,
the update throws).

## 3. Waiters: a dependency with no cell

A ref in the tree can have no cell (below a Null object or an error: reconcile.md §2), and a dependency derived from
params outlives the cell it named. Then the dependent **waits** for the ref (`DependencyWaiters`, per `State`):

- it is not linked, and `resolveDependencies` reports it unresolved: the cell is in error (`Unresolved dependency.`);
- when a cell of that ref is made, the waiters are linked to it, and are dependents in that pass (reconcile.md §4);
- a deleted cell's dependents wait for its ref (`unlinkCell`) and are evaluated again in the same update (their objects
  were made from it), so they are in error until the ref comes back, when they are linked again;
- a waiter that stops depending on the ref (`relinkCells`) or is deleted stops waiting.

## 4. Deferral

A cell is evaluated after its dependencies when they are evaluated in the same pass, whatever the order of the roots:
`updateSubtree` defers a cell (and its subtree) while a dependency is scheduled and still `pending` (`pendingBlockers`).
Deferred cells are retried after the roots, until a round makes no progress; then the update throws
`Unresolved dependency: <ref> (waiting on: <refs>)`. A dependency evaluated without an object (an error, below a Null
object) does not block: the dependent evaluates, and fails alone (`Unresolved dependency.`).

## Decisions

- **Dependencies are links between cells, not tree edges.** The tree stays a tree (parents, order, decorators); the
  links live with the cells and are rebuilt from the transforms, so params can name dependencies.
- **A missing cell is a wait, not a throw.** A dependency below a Null object, or a deleted dependency kept by params,
  is a normal state of an interactive session: the dependent alone is in error, and recovers when the cell is back. A
  ref that is not in the tree at all is a mistake in the edit: the update throws.
- **Only pending dependencies defer.** A dependency that was evaluated and has no object will not get one in this pass;
  deferring for it would only fail the whole update.
