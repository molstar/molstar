/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { deepEqual } from './index.js';

/**
 * Repeated immutable updates of the same paths of a plain JSON-like object, e.g. the frames of an animation
 * (MVS: `packages/mvs/runtime/_spec/animation.md §2`).
 *
 * The paths are declared once, in an `ObjectMutator.Plan` (a trie of keys, its leaves are "slots"); then each update
 * (a frame) only sets values of slots and commits. A commit copies the objects on the paths whose values changed since
 * the previous commit and keeps the identity of everything else, the previous commit's included: a commit with no
 * change returns the previous root itself, and unchanged subtrees are the same objects frame after frame.
 *
 * ```
 * const plan = ObjectMutator.plan(source);
 * plan.pushPrefix('children', 0, 'params');
 * const opacity = plan.slot('opacity');
 * plan.popPrefix();
 * const mutator = ObjectMutator.create(plan);
 * mutator.set(opacity, 0.5);
 * const next = mutator.commit(); // a copy of source with children[0].params.opacity = 0.5
 * ```
 *
 * A slot not set before a commit has the source's value at that path (the key is removed if the source has none). A
 * slot is "dead" when an object on its path is missing in the source: setting it does nothing. When a slot lies inside
 * another slot's value (the plan is not disjoint), the sets of a frame are applied in order, each copying from the
 * source, and the frame is compared with the previous commit as a whole (`deepEqual`).
 */
export interface ObjectMutator<T> {
  readonly plan: ObjectMutatorPlan<T>;
  /** The value of `slot` in the frame being built: the last value set since the previous commit, else the source's. */
  get(slot: number): unknown;
  set(slot: number, value: unknown): void;
  /** The object with this frame's values; the previous commit's root (or the source) when nothing differs from it. */
  commit(): T;
}

type Key = string | number;

const Absent = Symbol('absent');

export class ObjectMutatorPlan<T> {
  /** Trie nodes: 0 is the root; `keys[i]` is the key of node `i` in the object of its parent. */
  readonly keys: Key[] = [''];
  readonly parents: number[] = [-1];
  readonly children: number[][] = [[]];
  readonly childByKey: Map<string, number>[] = [new Map()];
  /** Slot of a trie node, -1 for inner nodes. */
  readonly slotOfNode: number[] = [-1];
  /** Trie node of a slot. */
  readonly slotNodes: number[] = [];
  /** Keys from the root to the slot. */
  readonly slotPaths: Key[][] = [];
  private prefix: Key[] = [];
  private prefixLengths: number[] = [];

  constructor(readonly source: T) {}

  pushPrefix(...keys: Key[]) {
    this.prefixLengths.push(keys.length);
    this.prefix.push(...keys);
  }

  popPrefix() {
    const n = this.prefixLengths.pop();
    if (n === undefined) throw new Error('ObjectMutator.Plan: popPrefix without pushPrefix');
    this.prefix.length -= n;
  }

  /** The slot at the current prefix followed by `path`; the same path gives the same slot. */
  slot(path: Key | Key[]): number {
    const keys = [...this.prefix, ...(Array.isArray(path) ? path : [path])];
    let node = 0;
    for (const key of keys) {
      let child = this.childByKey[node].get(String(key));
      if (child === undefined) {
        child = this.keys.length;
        this.keys.push(key);
        this.parents.push(node);
        this.children.push([]);
        this.childByKey.push(new Map());
        this.slotOfNode.push(-1);
        this.children[node].push(child);
        this.childByKey[node].set(String(key), child);
      }
      node = child;
    }
    if (this.slotOfNode[node] >= 0) return this.slotOfNode[node];
    const slot = this.slotNodes.length;
    this.slotOfNode[node] = slot;
    this.slotNodes.push(node);
    this.slotPaths.push(keys);
    return slot;
  }

  /** No slot lies inside another slot's value. */
  get isDisjoint() {
    for (const node of this.slotNodes) {
      if (this.children[node].length > 0) return false;
    }
    return true;
  }
}

function shallowCopy(obj: any): any {
  return Array.isArray(obj) ? obj.slice() : { ...obj };
}

function isObject(x: unknown): x is Record<Key, unknown> {
  return typeof x === 'object' && x !== null;
}

/** The value at `path` in `root`, `Absent` when the last key is not there, `undefined` when the slot is dead. */
function resolve(root: unknown, path: Key[]): { value: unknown; dead: boolean } {
  let f: unknown = root;
  for (let i = 0; i < path.length - 1; i++) {
    if (!isObject(f)) return { value: undefined, dead: true };
    f = f[path[i]];
  }
  if (!isObject(f)) return { value: undefined, dead: true };
  const key = path[path.length - 1];
  return { value: key in f ? f[key] : Absent, dead: false };
}

/** Slots disjoint: values are kept per slot, a commit rebuilds the paths of the slots that changed. */
class DisjointMutator<T> implements ObjectMutator<T> {
  private readonly sourceValues: unknown[];
  private readonly dead: boolean[];
  private readonly current: unknown[];
  private readonly pending: unknown[];
  private readonly isSet: Uint8Array;
  private readonly dirty: Uint8Array;
  private root: T;

  constructor(
    readonly plan: ObjectMutatorPlan<T>,
    private readonly equals: (a: unknown, b: unknown) => boolean,
  ) {
    const n = plan.slotNodes.length;
    this.sourceValues = new Array(n);
    this.dead = new Array(n);
    for (let s = 0; s < n; s++) {
      const r = resolve(plan.source, plan.slotPaths[s]);
      this.sourceValues[s] = r.value;
      this.dead[s] = r.dead;
    }
    this.current = this.sourceValues.slice();
    this.pending = new Array(n);
    this.isSet = new Uint8Array(n);
    this.dirty = new Uint8Array(plan.keys.length);
    this.root = plan.source;
  }

  get(slot: number) {
    const v = this.isSet[slot] ? this.pending[slot] : this.sourceValues[slot];
    return v === Absent ? undefined : v;
  }

  set(slot: number, value: unknown) {
    this.pending[slot] = value;
    this.isSet[slot] = 1;
  }

  commit(): T {
    const { plan, dirty } = this;
    let changed = false;
    for (let s = 0, n = plan.slotNodes.length; s < n; s++) {
      const value = this.isSet[s] ? this.pending[s] : this.sourceValues[s];
      this.isSet[s] = 0;
      this.pending[s] = undefined;
      if (this.dead[s] || value === this.current[s] || this.equals(value, this.current[s])) continue;
      this.current[s] = value;
      changed = true;
      for (let node = plan.slotNodes[s]; node >= 0 && !dirty[node]; node = plan.parents[node]) dirty[node] = 1;
    }
    if (!changed) return this.root;
    this.root = this.rebuild(0, this.root);
    return this.root;
  }

  private rebuild(node: number, obj: any): any {
    const { plan, dirty } = this;
    dirty[node] = 0;
    const copy = shallowCopy(obj);
    for (const child of plan.children[node]) {
      if (!dirty[child]) continue;
      const key = plan.keys[child];
      const slot = plan.slotOfNode[child];
      if (slot >= 0) {
        dirty[child] = 0;
        const value = this.current[slot];
        if (value === Absent) delete copy[key];
        else copy[key] = value;
      } else {
        copy[key] = this.rebuild(child, obj[key]);
      }
    }
    return copy;
  }
}

/** Slots inside other slots' values: each frame applies its sets in order to a copy-on-write of the source. */
class SequentialMutator<T> implements ObjectMutator<T> {
  private working: any;
  private copied = new Set<object>();
  private root: T;

  constructor(readonly plan: ObjectMutatorPlan<T>) {
    this.working = plan.source;
    this.root = plan.source;
  }

  get(slot: number) {
    let f: any = this.working;
    for (const key of this.plan.slotPaths[slot]) {
      if (!f) break;
      f = f[key];
    }
    return f;
  }

  set(slot: number, value: unknown) {
    const path = this.plan.slotPaths[slot];
    if (!this.copied.has(this.working)) {
      this.working = shallowCopy(this.working);
      this.copied.add(this.working);
    }
    let f = this.working;
    for (let i = 0; i < path.length - 1; i++) {
      const next = f[path[i]];
      if (!isObject(next)) return;
      if (!this.copied.has(next)) {
        const copy = shallowCopy(next);
        this.copied.add(copy);
        f[path[i]] = copy;
        f = copy;
      } else {
        f = next;
      }
    }
    f[path[path.length - 1]] = value;
  }

  commit(): T {
    const next = this.copied.size > 0 ? this.working : this.plan.source;
    this.working = this.plan.source;
    this.copied.clear();
    if (next !== this.root && !deepEqual(next, this.root)) this.root = next;
    return this.root;
  }
}

export const ObjectMutator = {
  plan<T>(source: T): ObjectMutator.Plan<T> {
    return new ObjectMutatorPlan(source);
  },
  /** `equals` compares a slot's value with the previous one (`deepEqual` by default); the plan not disjoint, whole trees */
  create<T>(plan: ObjectMutator.Plan<T>, options?: { equals?: (a: unknown, b: unknown) => boolean }): ObjectMutator<T> {
    return plan.isDisjoint ? new DisjointMutator(plan, options?.equals ?? deepEqual) : new SequentialMutator(plan);
  },
};

export namespace ObjectMutator {
  export type Plan<T> = ObjectMutatorPlan<T>;
}
