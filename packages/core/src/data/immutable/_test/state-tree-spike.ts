/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `StateTree` spike (README §9): type-check only, never run, no code of the state tree changes. Part 1 assigns
// these collections to `StateTree`'s interfaces as they are: each `@ts-expect-error` is a misfit. Part 2 restates the
// interfaces on this API and writes `StateTree`'s and `TransientTree`'s collection code against it. Checked with a
// scratch tsconfig (`tsc --noEmit -p`), since `tsc -b` leaves `_test/` out.

import type { StateTransform } from '../../../state/transform.js';
import type { StateTree } from '../../../state/tree/immutable.js';
import {
  HashMap,
  type MutableHashMap,
  type MutableOrderedHashSet,
  OrderedHashSet,
  type ReadonlyHashMap,
  type ReadonlyOrderedHashSet,
} from '../../immutable.js';

type Ref = string;
interface Transform {
  ref: Ref;
  parent: Ref;
  dependsOn?: Ref[];
}

// --- Part 1: as the interfaces are ---

// `ChildSet.forEach` is immutable.js's (returns the count, its callback gets the collection), and `map` is missing;
// `asMutable()` gives a set with `delete`, not `remove`.
// @ts-expect-error forEach's return type, map, remove
export const childSet: StateTree.ChildSet = OrderedHashSet.empty<Ref>();

// `_Map.get` returns `T` (not `T | undefined`), and a persistent `_Map` has `asImmutable()`.
// @ts-expect-error get's undefined, asImmutable on a persistent map
export const transforms: StateTree.Transforms = HashMap.empty<Ref, any>();

// member by member
// @ts-expect-error `map` is not in the agreed surface
export const mapFits: Pick<StateTree.ChildSet, 'map'> = OrderedHashSet.empty<Ref>();
// @ts-expect-error immutable.js's `forEach` returns the number of calls; ours returns nothing
export const forEachFits: Pick<StateTree.ChildSet, 'forEach'> = OrderedHashSet.empty<Ref>();
export const readsFit: Pick<StateTree.ChildSet, 'size' | 'has' | 'values' | 'toArray' | 'first'> =
  OrderedHashSet.empty<Ref>();
// @ts-expect-error `remove` is `delete` here
export const removeFits: Pick<StateTree.MutableChildSet, 'remove'> = OrderedHashSet.empty<Ref>().asMutable();
// @ts-expect-error `add` returns this set, which lacks `remove` and `map`
export const addFits: Pick<StateTree.MutableChildSet, 'add'> = OrderedHashSet.empty<Ref>().asMutable();
export const mutableSetReads: Pick<StateTree.MutableChildSet, 'has' | 'size' | 'toArray'> = OrderedHashSet.empty<Ref>()
  .asMutable()
  .add('x');
// @ts-expect-error `get` returns `T | undefined` here
export const getFits: Pick<StateTree.Transforms, 'get'> = HashMap.empty<Ref, StateTransform>();
export const hasFits: Pick<StateTree.Transforms, 'has' | 'size'> = HashMap.empty<Ref, StateTransform>();
// @ts-expect-error a persistent map has no `asImmutable()` here
export const asImmutableFits: Pick<StateTree.Transforms, 'asImmutable'> = HashMap.empty<Ref, StateTransform>();
// @ts-expect-error `set(...)` returns this map, whose `asImmutable()` returns a map without `asImmutable()`
export const mutableTransforms: Pick<StateTree.MutableTransforms, 'set'> = HashMap.empty<
  Ref,
  StateTransform
>().asMutable();
export const mutableReads: Pick<StateTree.MutableTransforms, 'has' | 'size'> = HashMap.empty<
  Ref,
  StateTransform
>().asMutable();

// behaviour and code outside the interfaces
export function behaviourMisfits(set: ReturnType<typeof OrderedHashSet.empty<Ref>>) {
  // immutable.js passes the collection as a third argument; ours passes two
  // @ts-expect-error the callback's third parameter
  set.forEach((_v: Ref, _k: Ref, _iter: unknown) => {});
  // immutable.js stops a `forEach` when its callback returns `false`; ours visits every entry (no type error: a
  // `false`-returning visitor in `StateTree` would walk on silently; none does today, `doPreOrder` stops by itself)
  set.forEach(() => false);
  // `dump` (state/tree/immutable.ts) reads `keySeq()` and `valueSeq()`: `keys()` and `values()` here
  // @ts-expect-error no lazy sequences
  HashMap.empty<Ref, Transform>().keySeq();
  // `transient.ts` casts `children.get(parent)` to immutable.js's `OrderedSet` (lines 128, 140): the casts move to
  // `OrderedHashSet`, or go once `ChildSet` is this API's type
  const t = OrderedHashSet.empty<Ref>().asMutable();
  // immutable.js's `withMutations` on a mutable collection edits it and keeps it mutable; a transient here has none,
  // and the implementation's (reached through a persistent type) would freeze the handle (README §3)
  // @ts-expect-error no withMutations on a transient
  t.withMutations(() => {});
}

// --- Part 2: the interfaces restated on this API ---

namespace Spike {
  export type ChildSet = ReadonlyOrderedHashSet<Ref> & { asMutable(): MutableChildSet };
  export type MutableChildSet = MutableOrderedHashSet<Ref>;
  export type Transforms = ReadonlyHashMap<Ref, Transform> & { asMutable(): MutableHashMap<Ref, Transform> };
  export type Children = ReadonlyHashMap<Ref, ChildSet> & { asMutable(): MutableHashMap<Ref, ChildSet> };

  export interface Tree {
    transforms: Transforms;
    children: Children;
    dependencies: Children;
  }

  // StateTree.createEmpty
  export function createEmpty(root: Transform): Tree {
    return {
      transforms: HashMap.of<Ref, Transform>([[root.ref, root]]),
      children: HashMap.of<Ref, ChildSet>([[root.ref, OrderedHashSet.empty<Ref>()]]),
      dependencies: HashMap.empty<Ref, ChildSet>(),
    };
  }

  // the visitors: `forEach` with a context
  type Ctx = { tree: Tree; f: (t: Transform) => void };
  function _postOrderFunc(this: Ctx, c: Ref) {
    _doPostOrder(this, this.tree.transforms.get(c)!);
  }
  function _doPostOrder(ctx: Ctx, root: Transform) {
    const children = ctx.tree.children.get(root.ref);
    if (children && children.size) children.forEach(_postOrderFunc, ctx);
    ctx.f(root);
  }

  // subtreeHasRef: `values()` (allocates; a cursor would not)
  export function subtreeHasRef(tree: Tree, root: Ref, ref: Ref): boolean {
    if (root === ref) return true;
    const it = tree.children.get(root)!.values();
    for (let n = it.next(); !n.done; n = it.next()) if (subtreeHasRef(tree, n.value, ref)) return true;
    return false;
  }

  // getDecoratorRoot: `first()`
  export function decoratorChild(tree: Tree, ref: Ref) {
    const children = tree.children.get(ref)!;
    return children.size === 1 ? tree.transforms.get(children.first()!) : undefined;
  }

  // the `ChildSet.map` callers (`children.map(...)`) become a cursor or `toArray().map`
  export function mapChildren<T>(set: ChildSet, f: (r: Ref) => T): T[] {
    return set.toArray().map(f);
  }

  // fromJSON: transients of maps and sets, frozen after
  export function fromJSON(data: Transform[]): Tree {
    const nodes = HashMap.empty<Ref, Transform>().asMutable();
    const children = HashMap.empty<Ref, MutableChildSet | ChildSet>().asMutable();
    const dependencies = HashMap.empty<Ref, MutableChildSet | ChildSet>().asMutable();
    for (const t of data) {
      nodes.set(t.ref, t);
      if (!children.has(t.ref)) children.set(t.ref, OrderedHashSet.empty<Ref>().asMutable());
      if (t.ref !== t.parent) (children.get(t.parent) as MutableChildSet).add(t.ref);
    }
    for (const t of data) {
      children.update(t.ref, (s) => (s as MutableChildSet).asImmutable());
      for (const d of t.dependsOn ?? []) {
        const s = dependencies.get(d) as MutableChildSet | undefined;
        if (s) s.add(t.ref);
        else dependencies.set(d, OrderedHashSet.of([t.ref]).asMutable());
      }
    }
    // `ReadonlyHashMap` is not a variance-safe view of a map of either set: the values are frozen above
    return {
      transforms: nodes.asImmutable(),
      children: children.asImmutable() as unknown as Children,
      dependencies: dependencies.asImmutable() as unknown as Children,
    };
  }

  // TransientTree: the maps turn mutable once; `asImmutable()` returns a new map and the transient handle is dead, so
  // the tree reassigns from the result and calls `asMutable()` before its next write. `renewOwner` goes: deletes write
  // no shared node (README §3).
  export class Transient {
    transforms: Transforms | MutableHashMap<Ref, Transform>;
    children: Children | MutableHashMap<Ref, ChildSet>;
    private changedChildren = false;
    private childMutations = new Map<Ref, MutableChildSet>();

    constructor(private tree: Tree) {
      this.transforms = tree.transforms;
      this.children = tree.children;
    }

    private changeChildren(): MutableHashMap<Ref, ChildSet> {
      if (!this.changedChildren) {
        this.changedChildren = true;
        this.children = this.children.asMutable();
      }
      return this.children as MutableHashMap<Ref, ChildSet>;
    }

    addChild(parent: Ref, child: Ref) {
      const m = this.childMutations.get(parent);
      if (m) {
        m.add(child);
        return;
      }
      const set = this.children.get(parent)!.asMutable();
      set.add(child);
      // a transient set stored in the map until `asImmutable()` fixes it: the map's value type has to admit it
      this.changeChildren().set(parent, set as unknown as ChildSet);
      this.childMutations.set(parent, set);
    }

    removeChild(parent: Ref, child: Ref) {
      const set = this.childMutations.get(parent) ?? this.children.get(parent)!.asMutable();
      set.delete(child);
      this.changeChildren().set(parent, set as unknown as ChildSet);
      this.childMutations.set(parent, set);
    }

    asImmutable(): Tree {
      if (!this.changedChildren) return this.tree;
      const children = this.changeChildren();
      this.childMutations.forEach((m, k) => children.set(k, m.asImmutable()));
      const tree: Tree = {
        transforms: this.tree.transforms,
        children: children.asImmutable(),
        dependencies: this.tree.dependencies,
      };
      this.tree = tree;
      this.children = tree.children;
      this.changedChildren = false;
      this.childMutations = new Map();
      return tree;
    }
  }
}

export { Spike };
