/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A persistent vector (README §2): Clojure's 32-way trie with a tail of up to 32. `get`, `set`, `push` and `pop` are
// effectively O(1); `shift` and `remove(index)` move every later element (O(n)) until the RRB vector.

import {
  type CollectionConfig,
  type Cursor,
  deadHandle,
  hashOrdered,
  newOwner,
  type Owner,
  type Resolved,
  resolveConfig,
  sameConfig,
} from './config.js';
import { asCursor } from './cursor.js';

export interface ReadonlyVector<T> {
  readonly size: number;
  get(index: number): T | undefined;
  first(): T | undefined;
  last(): T | undefined;
  /** The first index of `value` by the vector's `valueEquals` (default `Object.is`), or −1. */
  indexOf(value: T): number;
  find(pred: (value: T, index: number) => boolean): T | undefined;
  forEach<C>(fn: (this: C, value: T, index: number) => void, ctx?: C): void;
  /** Resets `c` (or a new cursor) to the first element (README §5). */
  cursor(c?: Cursor<T>): Cursor<T>;
  /** Allocates. */
  toArray(): T[];
  equals(other: ReadonlyVector<T>): boolean;
  hashCode(): number;
}

export interface Vector<T> extends ReadonlyVector<T> {
  /** Replaces the element at `index` in `[0, size)`; throws a `RangeError` outside. */
  set(index: number, value: T): Vector<T>;
  push(value: T): Vector<T>;
  pop(): Vector<T>;
  /** Without the first element: O(n). */
  shift(): Vector<T>;
  /** Without the element at `index`: O(n − index). */
  remove(index: number): Vector<T>;
  asMutable(): MutableVector<T>;
  withMutations(fn: (v: MutableVector<T>) => void): Vector<T>;
}

export interface MutableVector<T> extends ReadonlyVector<T> {
  set(index: number, value: T): this;
  push(value: T): this;
  pop(): this;
  shift(): this;
  remove(index: number): this;
  asMutable(): this;
  /** A new persistent vector; this handle is dead (README §3). */
  asImmutable(): Vector<T>;
}

export type VectorConfig<T> = Pick<CollectionConfig<unknown, T>, 'valueEquals' | 'valueHash'>;

class VNode {
  constructor(
    public owner: Owner | null,
    public array: any[],
  ) {}
}

const EmptyNode = new VNode(null, []);

function editable(node: VNode, owner: Owner | null) {
  return owner !== null && node.owner === owner;
}

function tailOffset(size: number) {
  return size < 32 ? 0 : ((size - 1) >>> 5) << 5;
}

function newPath(owner: Owner | null, level: number, node: VNode): VNode {
  if (level === 0) return node;
  return new VNode(owner, [newPath(owner, level - 5, node)]);
}

export class VectorImpl<T> {
  /** The cached `hashCode` of a persistent vector. */
  private _hash: number | undefined = undefined;
  /** The tail is this transient's own array. */
  private tailOwned = false;
  private dead = false;

  constructor(
    public size: number,
    public shiftBits: number,
    public root: VNode,
    public tail: any[],
    public cfg: Resolved,
    public owner: Owner | null,
  ) {}

  /** The leaf array that holds `index` (in range). */
  arrayFor(index: number): any[] {
    if (index >= tailOffset(this.size)) return this.tail;
    let node = this.root;
    for (let level = this.shiftBits; level > 0; level -= 5) node = node.array[(index >>> level) & 31];
    return node.array;
  }

  get(index: number): T | undefined {
    if (index < 0 || index >= this.size || index !== (index | 0)) return undefined;
    return this.arrayFor(index)[index & 31];
  }

  first(): T | undefined {
    return this.size > 0 ? this.get(0) : undefined;
  }

  last(): T | undefined {
    return this.size > 0 ? this.tail[this.tail.length - 1] : undefined;
  }

  indexOf(value: T): number {
    const eq = this.cfg.valueEquals;
    for (let i = 0; i < this.size; i += 32) {
      const a = this.arrayFor(i);
      for (let j = 0; j < a.length; j++) if (eq(a[j], value)) return i + j;
    }
    return -1;
  }

  find(pred: (value: T, index: number) => boolean): T | undefined {
    for (let i = 0; i < this.size; i += 32) {
      const a = this.arrayFor(i);
      for (let j = 0; j < a.length; j++) if (pred(a[j], i + j)) return a[j];
    }
    return undefined;
  }

  forEach<C>(fn: (this: C, value: T, index: number) => void, ctx?: C): void {
    for (let i = 0; i < this.size; i += 32) {
      const a = this.arrayFor(i);
      for (let j = 0; j < a.length; j++) fn.call(ctx as C, a[j], i + j);
    }
  }

  cursor(c?: Cursor<T>): Cursor<T> {
    const cur = asCursor(c);
    cur.resetVector(this);
    return cur;
  }

  toArray(): T[] {
    const out = new Array<T>(this.size);
    for (let i = 0; i < this.size; i += 32) {
      const a = this.arrayFor(i);
      for (let j = 0; j < a.length; j++) out[i + j] = a[j];
    }
    return out;
  }

  equals(other: ReadonlyVector<T>): boolean {
    if ((other as unknown) === this) return true;
    if (!(other instanceof VectorImpl) || other.size !== this.size || !sameConfig(this.cfg, other.cfg)) return false;
    const eq = this.cfg.valueEquals;
    for (let i = 0; i < this.size; i += 32) {
      const a = this.arrayFor(i);
      const b = other.arrayFor(i);
      if (a === b) continue;
      for (let j = 0; j < a.length; j++) if (!eq(a[j], b[j])) return false;
    }
    return true;
  }

  hashCode(): number {
    if (this._hash !== undefined) return this._hash;
    const vh = this.cfg.valueHash;
    let acc = this.size;
    if (vh) {
      for (let i = 0; i < this.size; i += 32) {
        const a = this.arrayFor(i);
        for (let j = 0; j < a.length; j++) acc = hashOrdered(acc, vh(a[j]));
      }
    }
    if (this.owner === null) this._hash = acc;
    return acc;
  }

  // --- edits ---

  private result(size: number, shiftBits: number, root: VNode, tail: any[]): any {
    if (this.owner === null) return new VectorImpl<T>(size, shiftBits, root, tail, this.cfg, null);
    if (this.dead) deadHandle();
    this.size = size;
    this.shiftBits = shiftBits;
    this.root = root;
    this.tail = tail;
    return this;
  }

  private check() {
    if (this.dead) deadHandle();
  }

  set(index: number, value: T): any {
    this.check();
    if (index < 0 || index >= this.size || index !== (index | 0)) {
      throw new RangeError(`Index ${index} out of range [0, ${this.size}).`);
    }
    const owner = this.owner;
    if (index >= tailOffset(this.size)) {
      const t = this.tail;
      const j = index & 31;
      if (this.cfg.valueEquals(t[j], value)) return this;
      if (owner !== null && this.tailOwned) {
        t[j] = value;
        return this;
      }
      const nt = t.slice();
      nt[j] = value;
      if (owner !== null) this.tailOwned = true;
      return this.result(this.size, this.shiftBits, this.root, nt);
    }
    if (this.cfg.valueEquals(this.arrayFor(index)[index & 31], value)) return this;
    const root = assoc(owner, this.shiftBits, this.root, index, value);
    return this.result(this.size, this.shiftBits, root, this.tail);
  }

  push(value: T): any {
    this.check();
    const owner = this.owner;
    const size = this.size;
    if (size - tailOffset(size) < 32) {
      if (owner !== null && this.tailOwned) {
        this.tail.push(value);
        this.size++;
        return this;
      }
      const nt = this.tail.slice();
      nt.push(value);
      if (owner !== null) this.tailOwned = true;
      return this.result(size + 1, this.shiftBits, this.root, nt);
    }
    // the full tail goes into the tree
    // a transient's tail it has not copied yet is shared: copied before the tree may write it in place
    const tailNode = new VNode(owner, owner !== null && !this.tailOwned ? this.tail.slice() : this.tail);
    let shiftBits = this.shiftBits;
    let root: VNode;
    if (size >>> 5 > 1 << shiftBits) {
      root = new VNode(owner, [this.root, newPath(owner, shiftBits, tailNode)]);
      shiftBits += 5;
    } else {
      root = pushTail(owner, size, shiftBits, this.root, tailNode);
    }
    if (owner !== null) this.tailOwned = true;
    return this.result(size + 1, shiftBits, root, [value]);
  }

  pop(): any {
    this.check();
    const owner = this.owner;
    const size = this.size;
    if (size === 0) return this;
    if (size === 1) {
      if (owner !== null) this.tailOwned = true;
      return this.result(0, 5, EmptyNode, []);
    }
    if (size - tailOffset(size) > 1) {
      if (owner !== null && this.tailOwned) {
        this.tail.pop();
        this.size--;
        return this;
      }
      const nt = this.tail.slice(0, this.tail.length - 1);
      if (owner !== null) this.tailOwned = true;
      return this.result(size - 1, this.shiftBits, this.root, nt);
    }
    // the tail empties: the last leaf becomes the tail (copied, it may be shared)
    const newTail = this.arrayFor(size - 2).slice();
    let root = popTail(owner, size, this.shiftBits, this.root) ?? EmptyNode;
    let shiftBits = this.shiftBits;
    if (shiftBits > 5 && root.array.length === 1) {
      root = root.array[0];
      shiftBits -= 5;
    }
    if (owner !== null) this.tailOwned = true;
    return this.result(size - 1, shiftBits, root, newTail);
  }

  remove(index: number): any {
    this.check();
    if (index < 0 || index >= this.size || index !== (index | 0)) {
      throw new RangeError(`Index ${index} out of range [0, ${this.size}).`);
    }
    if (index === this.size - 1) return this.pop();
    const m: VectorImpl<T> = this.owner === null ? (this.asMutable() as any) : this;
    for (let i = index; i < m.size - 1; i++) m.setUnchecked(i, m.get(i + 1)!);
    m.pop();
    return this.owner === null ? m.asImmutable() : m;
  }

  /** Sets `index` and `index + 1` (an even `index`: one leaf) in one path copy, an ordered map's slot (README §2). */
  setPair(index: number, a: any, b: any): any {
    this.check();
    const owner = this.owner;
    if (index >= tailOffset(this.size)) {
      let t = this.tail;
      if (owner === null || !this.tailOwned) {
        t = t.slice();
        if (owner !== null) this.tailOwned = true;
      }
      t[index & 31] = a;
      t[(index & 31) + 1] = b;
      return this.result(this.size, this.shiftBits, this.root, t);
    }
    const root = assocPair(owner, this.shiftBits, this.root, index, a, b);
    return this.result(this.size, this.shiftBits, root, this.tail);
  }

  /** `set` without the no-op check: `remove` moves values that may be equal. */
  private setUnchecked(index: number, value: T) {
    if (index >= tailOffset(this.size)) {
      if (!this.tailOwned) {
        this.tail = this.tail.slice();
        this.tailOwned = true;
      }
      this.tail[index & 31] = value;
      return;
    }
    this.root = assoc(this.owner, this.shiftBits, this.root, index, value);
  }

  shift(): any {
    return this.remove(0);
  }

  asMutable(): any {
    if (this.owner !== null) {
      this.check();
      return this;
    }
    return new VectorImpl<T>(this.size, this.shiftBits, this.root, this.tail, this.cfg, newOwner());
  }

  asImmutable(): Vector<T> {
    if (this.owner === null) return this as any;
    this.check();
    this.dead = true;
    return new VectorImpl<T>(this.size, this.shiftBits, this.root, this.tail, this.cfg, null) as any;
  }

  withMutations(fn: (v: MutableVector<T>) => void): Vector<T> {
    const m = this.asMutable();
    fn(m);
    return m.asImmutable();
  }
}

function assocPair(owner: Owner | null, level: number, node: VNode, index: number, a: any, b: any): VNode {
  const n = editable(node, owner) ? node : new VNode(owner, node.array.slice());
  if (level === 0) {
    n.array[index & 31] = a;
    n.array[(index & 31) + 1] = b;
  } else {
    const i = (index >>> level) & 31;
    n.array[i] = assocPair(owner, level - 5, n.array[i], index, a, b);
  }
  return n;
}

function assoc(owner: Owner | null, level: number, node: VNode, index: number, value: any): VNode {
  const n = editable(node, owner) ? node : new VNode(owner, node.array.slice());
  if (level === 0) {
    n.array[index & 31] = value;
  } else {
    const i = (index >>> level) & 31;
    n.array[i] = assoc(owner, level - 5, n.array[i], index, value);
  }
  return n;
}

function pushTail(owner: Owner | null, size: number, level: number, parent: VNode, tailNode: VNode): VNode {
  const i = ((size - 1) >>> level) & 31;
  const n = editable(parent, owner) ? parent : new VNode(owner, parent.array.slice());
  let child: VNode;
  if (level === 5) {
    child = tailNode;
  } else {
    const sub: VNode | undefined = parent.array[i];
    child = sub ? pushTail(owner, size, level - 5, sub, tailNode) : newPath(owner, level - 5, tailNode);
  }
  n.array[i] = child;
  return n;
}

function popTail(owner: Owner | null, size: number, level: number, node: VNode): VNode | undefined {
  const i = ((size - 2) >>> level) & 31;
  if (level > 5) {
    const child = popTail(owner, size, level - 5, node.array[i]);
    if (child === undefined && i === 0) return undefined;
    const n = editable(node, owner) ? node : new VNode(owner, node.array.slice());
    if (child === undefined) n.array.length = i;
    else n.array[i] = child;
    return n;
  }
  if (i === 0) return undefined;
  const n = editable(node, owner) ? node : new VNode(owner, node.array.slice());
  n.array.length = i;
  return n;
}

const EmptyVector = new VectorImpl<any>(0, 5, EmptyNode, [], resolveConfig(undefined), null);

export const Vector = {
  empty<T>(config?: VectorConfig<T>): Vector<T> {
    const cfg = resolveConfig(config);
    if (cfg === EmptyVector.cfg) return EmptyVector as any;
    return new VectorImpl<T>(0, 5, EmptyNode, [], cfg, null) as any;
  },
  of<T>(values: Iterable<T>, config?: VectorConfig<T>): Vector<T> {
    const m = Vector.empty<T>(config).asMutable();
    for (const v of values) m.push(v);
    return m.asImmutable();
  },
};
