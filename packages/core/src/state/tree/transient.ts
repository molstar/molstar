/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { OrderedSet } from 'immutable';
import { StateTransform } from '../transform.js';
import { StateTree } from './immutable.js';
import { shallowEqual } from '@molstar/core/util/object';
import { arrayEqual } from '@molstar/core/util/array';

export { TransientTree };

class TransientTree implements StateTree {
  transforms = this.tree.transforms as StateTree.MutableTransforms;
  children = this.tree.children as StateTree.MutableChildren;
  dependencies = this.tree.dependencies as StateTree.MutableDependencies;

  private changedNodes = false;
  private changedChildren = false;
  private changedDependencies = false;
  /** a map that had deletes takes a new owner before its next write (`renewOwner`, tree.md §2) */
  private deletedNodes = false;
  private deletedChildren = false;
  private deletedDependencies = false;

  private _childMutations: Map<StateTransform.Ref, OrderedSet<StateTransform.Ref>> | undefined = void 0;
  private _dependencyMutations: Map<StateTransform.Ref, StateTree.MutableChildSet> | undefined = void 0;
  private _stateUpdates: Set<StateTransform.Ref> | undefined = void 0;

  /**
   * Set by a builder (`StateBuilder.EditInfo.edits`): every ref whose transform this tree sets, and whether a node was
   * moved or removed. The reconciler finds the update roots of such edits without walking the trees (reconcile.md §2).
   */
  edits: TransientTree.Edits | undefined = void 0;

  private edited(ref: StateTransform.Ref) {
    if (this.edits) this.edits.refs.push(ref);
  }

  private moved(ref: StateTransform.Ref) {
    if (this.edits) this.edits.moved = true;
    this.edited(ref);
  }

  private removed(ref: StateTransform.Ref) {
    if (this.edits) this.edits.removed.push(ref);
  }

  private get childMutations() {
    if (this._childMutations) return this._childMutations;
    this._childMutations = new Map();
    return this._childMutations;
  }

  private get dependencyMutations() {
    if (this._dependencyMutations) return this._dependencyMutations;
    this._dependencyMutations = new Map();
    return this._dependencyMutations;
  }

  // Every write to a map goes through its `change*` gate and every delete through its `delete*`: the map turns
  // mutable on the first write or delete, and a write after deletes renews its owner (tree.md §2); deletes in a row
  // stay in place.

  private changeNodes() {
    if (this.deletedNodes) {
      this.deletedNodes = false;
      this.transforms = renewOwner(this.transforms);
    } else if (!this.changedNodes) {
      this.changedNodes = true;
      this.transforms = this.transforms.asMutable();
    }
  }

  private changeChildren() {
    if (this.deletedChildren) {
      this.deletedChildren = false;
      this.children = renewOwner(this.children);
    } else if (!this.changedChildren) {
      this.changedChildren = true;
      this.children = this.children.asMutable();
    }
  }

  private changeDependencies() {
    if (this.deletedDependencies) {
      this.deletedDependencies = false;
      this.dependencies = renewOwner(this.dependencies);
    } else if (!this.changedDependencies) {
      this.changedDependencies = true;
      this.dependencies = this.dependencies.asMutable();
    }
  }

  private deleteNode(ref: StateTransform.Ref) {
    if (!this.deletedNodes) this.changeNodes();
    this.transforms.delete(ref);
    this.deletedNodes = true;
  }

  private deleteChildren(ref: StateTransform.Ref) {
    if (!this.deletedChildren) this.changeChildren();
    this.children.delete(ref);
    this.deletedChildren = true;
  }

  private deleteDependencies(ref: StateTransform.Ref) {
    if (!this.deletedDependencies) this.changeDependencies();
    this.dependencies.delete(ref);
    this.deletedDependencies = true;
  }

  get root() {
    return this.transforms.get(StateTransform.RootRef)!;
  }

  asTransient() {
    return this.asImmutable().asTransient();
  }

  private addChild(parent: StateTransform.Ref, child: StateTransform.Ref) {
    if (this.childMutations.has(parent)) {
      this.childMutations.get(parent)!.add(child);
    } else {
      const set = (this.children.get(parent) as OrderedSet<StateTransform.Ref>).asMutable();
      set.add(child);
      this.changeChildren();
      this.children.set(parent, set);
      this.childMutations.set(parent, set);
    }
  }

  private removeChild(parent: StateTransform.Ref, child: StateTransform.Ref) {
    if (this.childMutations.has(parent)) {
      this.childMutations.get(parent)!.remove(child);
    } else {
      const set = (this.children.get(parent) as OrderedSet<StateTransform.Ref>).asMutable();
      set.remove(child);
      this.changeChildren();
      this.children.set(parent, set);
      this.childMutations.set(parent, set);
    }
  }

  private clearRoot() {
    const parent = StateTransform.RootRef;
    if (this.children.get(parent).size === 0) return;

    this.changeChildren();

    // mutable: later `addChild`s to the root write into it
    const set = OrderedSet<StateTransform.Ref>().asMutable();
    this.children.set(parent, set);
    this.childMutations.set(parent, set);
  }

  private mutateDependency(parent: StateTransform.Ref, child: StateTransform.Ref, action: 'add' | 'remove') {
    let set: StateTree.MutableChildSet | undefined = this.dependencyMutations.get(parent);

    if (!set) {
      const src = this.dependencies.get(parent);
      if (!src && action === 'remove') return;
      this.changeDependencies();
      set = src ? src.asMutable() : OrderedSet<string>().asMutable();
      this.dependencyMutations.set(parent, set);
      this.dependencies.set(parent, set);
    }

    if (action === 'add') {
      set.add(child);
    } else {
      set.remove(child);
    }
  }

  changeParent(ref: StateTransform.Ref, newParent: StateTransform.Ref) {
    ensurePresent(this.transforms, ref);

    const old = this.transforms.get(ref);
    this.removeChild(old.parent, ref);
    this.addChild(newParent, ref);
    this.changeNodes();
    this.transforms.set(ref, StateTransform.withParent(old, newParent));
    this.moved(ref);
  }

  /** Keeps `dependencies` in step with a node's `dependsOn` going from `from` to `to` */
  private relinkDependencies(
    ref: StateTransform.Ref,
    from: readonly string[] | undefined,
    to: readonly string[] | undefined,
  ) {
    if (from) {
      for (const d of from) if (!to?.includes(d)) this.mutateDependency(d, ref, 'remove');
    }
    if (to) {
      for (const d of to) if (!from?.includes(d)) this.mutateDependency(d, ref, 'add');
    }
  }

  add(transform: StateTransform) {
    const ref = transform.ref;

    if (this.transforms.has(transform.ref)) {
      const node = this.transforms.get(transform.ref);
      if (node.parent !== transform.parent) alreadyPresent(transform.ref);
      // a replaced node keeps no dependency it no longer has (perf-audit-2 03/B4)
      this.relinkDependencies(ref, node.dependsOn, transform.dependsOn);
    }

    const children = this.children.get(transform.parent);
    if (!children) parentNotPresent(transform.parent);

    if (!children.has(transform.ref)) {
      this.addChild(transform.parent, transform.ref);
    }

    if (!this.children.has(transform.ref)) {
      this.changeChildren();
      this.children.set(transform.ref, OrderedSet());
    }

    this.changeNodes();
    this.transforms.set(ref, transform);
    this.edited(ref);

    if (transform.dependsOn) {
      for (const d of transform.dependsOn) {
        this.mutateDependency(d, ref, 'add');
      }
    }

    return this;
  }

  /** Calls Transform.definition.params.areEqual if available, otherwise uses shallowEqual to check if the params changed */
  setParams(ref: StateTransform.Ref, params: any) {
    ensurePresent(this.transforms, ref);

    const transform = this.transforms.get(ref)!;
    // TODO: should this be here?
    if (shallowEqual(transform.params, params)) {
      return false;
    }

    this.changeNodes();

    this.transforms.set(transform.ref, StateTransform.withParams(transform, params));
    this.edited(ref);
    return true;
  }

  /** Calls Transform.definition.params.areEqual if available, otherwise uses shallowEqual to check if the params changed */
  setTags(ref: StateTransform.Ref, tags: string | string[] | undefined) {
    ensurePresent(this.transforms, ref);

    const transform = this.transforms.get(ref)!;

    const withTags = StateTransform.withTags(transform, tags);
    // TODO: should this be here?
    if (arrayEqual(transform.tags, withTags.tags)) {
      return false;
    }

    this.changeNodes();

    this.transforms.set(transform.ref, withTags);
    this.edited(ref);
    return true;
  }

  setDependsOn(ref: StateTransform.Ref, dependsOn: string | string[] | undefined) {
    ensurePresent(this.transforms, ref);

    const transform = this.transforms.get(ref)!;

    const withDependsOn = StateTransform.withDependsOn(transform, dependsOn);
    if (arrayEqual(transform.dependsOn, withDependsOn.dependsOn)) {
      return false;
    }

    this.changeNodes();

    this.transforms.set(transform.ref, withDependsOn);
    // deleting a new dependency cascades to this node (perf-audit-2 03/B4)
    this.relinkDependencies(ref, transform.dependsOn, withDependsOn.dependsOn);
    this.edited(ref);
    return true;
  }

  assignState(ref: StateTransform.Ref, state?: Partial<StateTransform.State>) {
    ensurePresent(this.transforms, ref);

    const old = this.transforms.get(ref);
    if (this._stateUpdates && this._stateUpdates.has(ref)) {
      StateTransform.assignState(old.state, state);
      return old;
    } else {
      if (!this._stateUpdates) this._stateUpdates = new Set();
      this._stateUpdates.add(old.ref);
      this.changeNodes();
      const updated = StateTransform.withState(old, state);
      this.transforms.set(ref, updated);
      this.edited(ref);
      return updated;
    }
  }

  remove(ref: StateTransform.Ref): StateTransform[] {
    const node = this.transforms.get(ref);
    if (!node) return [];

    const st = StateTree.subtreePostOrder(this, node);
    if (ref === StateTransform.RootRef) {
      st.pop();
      if (st.length === 0) return st;
      this.clearRoot();
    } else {
      if (st.length === 0) return st;
      this.removeChild(node.parent, node.ref);
    }

    this.removed(ref);

    for (const n of st) {
      this.deleteNode(n.ref);
      this.deleteChildren(n.ref);
      if (this._childMutations) this._childMutations.delete(n.ref);
    }

    const depRemoves: StateTransform[] = [];
    for (const n of st) {
      if (n.dependsOn) {
        for (const d of n.dependsOn) {
          if (!this.transforms.has(d)) continue;
          this.mutateDependency(d, n.ref, 'remove');
        }
      }

      if (this.dependencies.has(n.ref)) {
        const deps = this.dependencies.get(n.ref).toArray();
        this.deleteDependencies(n.ref);
        if (this._dependencyMutations) this._dependencyMutations.delete(n.ref);

        for (const dep of deps) {
          if (!this.transforms.has(dep)) continue;
          for (const del of this.remove(dep)) depRemoves[depRemoves.length] = del;
        }
      }
    }

    for (const dep of depRemoves) st[st.length] = dep;

    return st;
  }

  /**
   * The result becomes this tree's new base (tree.md §2): `asImmutable()` freezes the mutable maps in place, so later
   * writes start from fresh mutable copies, and flags merged in place (`_stateUpdates`) never reach the snapshot.
   */
  asImmutable() {
    if (!this.changedNodes && !this.changedChildren && !this.changedDependencies && !this._childMutations) {
      return this.tree;
    }
    if (this._childMutations) {
      this.changeChildren();
      fixChildMutations(this._childMutations, this.children);
    }
    if (this._dependencyMutations) {
      this.changeDependencies();
      fixDependencyMutations(this._dependencyMutations, this.dependencies);
    }
    const tree = StateTree.create(
      this.changedNodes ? this.transforms.asImmutable() : this.transforms,
      this.changedChildren ? this.children.asImmutable() : this.children,
      this.changedDependencies ? this.dependencies.asImmutable() : this.dependencies,
    );

    this.tree = tree;
    this.transforms = tree.transforms as StateTree.MutableTransforms;
    this.children = tree.children as StateTree.MutableChildren;
    this.dependencies = tree.dependencies as StateTree.MutableDependencies;
    this.changedNodes = false;
    this.changedChildren = false;
    this.changedDependencies = false;
    this.deletedNodes = false;
    this.deletedChildren = false;
    this.deletedDependencies = false;
    this._childMutations = void 0;
    this._dependencyMutations = void 0;
    this._stateUpdates = void 0;

    return tree;
  }

  constructor(private tree: StateTree) {}
}

namespace TransientTree {
  export interface Edits {
    /** refs whose transform was set (added, params, tags, dependencies, state flags, parent), in order, with repeats */
    refs: StateTransform.Ref[];
    /** a node was moved to another parent */
    moved: boolean;
    /** the refs `remove` was called with that removed nodes, in order (dependents removed with one come after it) */
    removed: StateTransform.Ref[];
  }
}

function fixChildMutations(
  mutations: Map<StateTransform.Ref, OrderedSet<StateTransform.Ref>>,
  children: StateTree.MutableMap<StateTree.ChildSet>,
) {
  mutations.forEach((m, k) => children.set(k, m.asImmutable()));
}

/** the sets, then the deletes: the map is frozen right after them, so no write follows a delete (`renewOwner`) */
function fixDependencyMutations(
  mutations: Map<StateTransform.Ref, StateTree.MutableChildSet>,
  dependencies: StateTree.MutableMap<StateTree.ChildSet>,
) {
  const empty: StateTransform.Ref[] = [];
  mutations.forEach((m, k) => {
    if (m.size === 0) empty.push(k);
    else dependencies.set(k, m.asImmutable());
  });
  for (const k of empty) dependencies.delete(k);
}

/**
 * A new owner for a mutable map, before a write that follows deletes (tree.md §2). Deleting one of two keys with equal
 * hashes leaves the other in a node of the map's owner that shares its `[key, value]` entry with the source tree's
 * node, and Immutable writes a later `set` of that key into the entry in place: into the source tree
 * (`HashCollisionNode.update`; seen in immutable 4.3.9, 5.1.9 and 6.0.0-beta.1; fixed by immutable-js#2279, see the
 * core TODO). Under a new owner, that `set` copies the node instead; the writes after it copy their paths once more.
 * Generated refs collide this way often.
 */
function renewOwner<M extends StateTree.MutableMap<any>>(map: M): M {
  return map.asImmutable().asMutable() as M;
}

function alreadyPresent(ref: StateTransform.Ref) {
  throw new Error(`Transform '${ref}' is already present in the tree.`);
}

function parentNotPresent(ref: StateTransform.Ref) {
  throw new Error(`Parent '${ref}' must be present in the tree.`);
}

function ensurePresent(nodes: StateTree.Transforms, ref: StateTransform.Ref) {
  if (!nodes.has(ref)) {
    throw new Error(`Node '${ref}' is not present in the tree.`);
  }
}
