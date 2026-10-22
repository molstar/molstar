/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The diff of a state update (reconcile.md §2): deletes and update roots over the whole trees, and the same found from a
// builder's edits. Internal to the reconciler (state.ts); not exported by the package's state index.

import { StateObject, type StateObjectCell } from '../object.js';
import { StateTransform } from '../transform.js';
import { StateTree } from '../tree/immutable.js';

type Ref = StateTransform.Ref;
type Cells = ReadonlyMap<Ref, StateObjectCell>;

/** refs of the old tree that are not in the new one and have a cell, in post-order (leaves first) */
export function findDeletes(oldTree: StateTree, newTree: StateTree, cells: Cells): Ref[] {
  const deleteCtx: FindDeletesCtx = { newTree, cells, deletes: [] };
  StateTree.doPostOrder(oldTree, oldTree.root, deleteCtx, checkDeleteVisitor);
  return deleteCtx.deletes;
}

type FindDeletesCtx = { newTree: StateTree; cells: Cells; deletes: Ref[] };
function checkDeleteVisitor(n: StateTransform, _: any, ctx: FindDeletesCtx) {
  if (!ctx.newTree.transforms.has(n.ref) && ctx.cells.has(n.ref)) ctx.deletes.push(n.ref);
}

/**
 * `findDeletes` when the builder only removed the nodes `removed` (`TransientTree.Edits.removed`, dependents removed
 * with one included; reconcile.md §2): only their subtrees in the old tree can be gone. Each topmost one's subtree, in
 * the old tree's pre-order of those roots: post-order of disjoint subtrees is a run of the whole tree's each, in that
 * order.
 */
export function findEditedDeletes(oldTree: StateTree, newTree: StateTree, cells: Cells, removed: readonly Ref[]) {
  const deleteCtx: FindDeletesCtx = { newTree, cells, deletes: [] };
  if (removed.length === 0) return deleteCtx.deletes;
  const removedSet = new Set<Ref>();
  for (const r of removed) if (oldTree.transforms.has(r)) removedSet.add(r);
  const paths = new Map<Ref, Ref[]>();
  for (const r of removedSet) {
    const path = pathOf(oldTree, r);
    // below another removed node: in that one's subtree
    let nested = false;
    for (let i = 0; i < path.length - 1; i++) {
      if (removedSet.has(path[i])) {
        nested = true;
        break;
      }
    }
    if (!nested) paths.set(r, path);
  }
  const roots = Array.from(paths.keys());
  if (roots.length > 1) {
    const order = childOrder(oldTree);
    roots.sort((a, b) => comparePreOrder(order, paths.get(a)!, paths.get(b)!));
  }
  for (const r of roots) StateTree.doPostOrder(oldTree, oldTree.transforms.get(r), deleteCtx, checkDeleteVisitor);
  return deleteCtx.deletes;
}

/**
 * The roots `findUpdateRoots` finds when the tree differs from the cells' only at `edited` (reconcile.md §2): an edited
 * node whose cell is missing or has another version, below no node the diff would stop at, in pre-order.
 */
export function findEditedRoots(cells: Cells, tree: StateTree, edited: readonly Ref[]) {
  const roots: Ref[] = [];
  const seen = new Set<Ref>();
  for (const ref of edited) {
    if (seen.has(ref)) continue;
    seen.add(ref);
    const n = tree.transforms.get(ref);
    if (!n) continue;
    const cell = cells.get(ref);
    if (cell && cell.transform.version === n.version) continue;
    if (reachedByDiff(cells, tree, n)) roots.push(ref);
  }
  if (roots.length > 1) {
    const paths = new Map(roots.map((r) => [r, pathOf(tree, r)]));
    const order = childOrder(tree);
    roots.sort((a, b) => comparePreOrder(order, paths.get(a)!, paths.get(b)!));
  }
  return roots;
}

/** whether every ancestor of `n` passes `findUpdateRootsVisitor`: the whole-tree diff gets to `n` */
function reachedByDiff(cells: Cells, tree: StateTree, n: StateTransform) {
  for (let t = n; t.ref !== StateTransform.RootRef; ) {
    t = tree.transforms.get(t.parent);
    const cell = cells.get(t.ref);
    if (!cell || cell.transform.version !== t.version) return false;
    if (cell.status === 'error' || cell.obj === StateObject.Null) return false;
  }
  return true;
}

/** refs from the root to `ref` */
function pathOf(tree: StateTree, ref: Ref) {
  const path: Ref[] = [];
  for (let t = tree.transforms.get(ref); ; t = tree.transforms.get(t.parent)) {
    path.push(t.ref);
    if (t.ref === StateTransform.RootRef) break;
  }
  return path.reverse();
}

/** a child's position among its parent's children, each parent's indexed once (one pass over its children) */
function childOrder(tree: StateTree) {
  const byParent = new Map<Ref, Map<Ref, number>>();
  return (parent: Ref, child: Ref) => {
    let index = byParent.get(parent);
    if (!index) {
      index = new Map();
      let i = 0;
      for (const c of tree.children.get(parent).values()) index.set(c, i++);
      byParent.set(parent, index);
    }
    return index.get(child)!;
  };
}

/** pre-order of two nodes: their paths from the root first differ at siblings, ordered as the children are */
function comparePreOrder(order: (parent: Ref, child: Ref) => number, a: Ref[], b: Ref[]) {
  let i = 0;
  while (i < a.length && i < b.length && a[i] === b[i]) i++;
  if (i === a.length || i === b.length) return a.length - b.length;
  return order(a[i - 1], a[i]) - order(a[i - 1], b[i]);
}

export function findUpdateRoots(cells: Cells, tree: StateTree) {
  const findState = { roots: [] as Ref[], cells };
  StateTree.doPreOrder(tree, tree.root, findState, findUpdateRootsVisitor);
  return findState.roots;
}

function findUpdateRootsVisitor(n: StateTransform, _: any, s: { roots: Ref[]; cells: Cells }) {
  const cell = s.cells.get(n.ref);
  if (!cell || cell.transform.version !== n.version) {
    s.roots.push(n.ref);
    return false;
  }
  if (cell.status === 'error') return false;

  // nothing below a Null object can be an update root
  if (cell && cell.obj === StateObject.Null) return false;
  return true;
}
