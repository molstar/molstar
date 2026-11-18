/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A persistent set in insertion order (README §2): a `HashMap` from key to slot and a `Vector` of keys, holes
// compacted as in the ordered map.

import { EmptyFlat, NotFound, rootGet } from './champ.js';
import {
  type Cursor,
  deadHandle,
  hashOrdered,
  newOwner,
  type Owner,
  type Resolved,
  resolveConfig,
  sameConfig,
} from './config.js';
import { asCursor, Hole } from './cursor.js';
import { HashMapImpl } from './hash-map.js';
import type { ReadonlyHashSet, SetConfig } from './hash-set.js';
import { CompactMin, indexConfig } from './ordered-hash-map.js';
import { Vector, type VectorImpl } from './vector.js';

export interface ReadonlyOrderedHashSet<K> extends Omit<ReadonlyHashSet<K>, 'equals'> {
  /** The first key in insertion order. */
  first(): K | undefined;
  /** The last key in insertion order. */
  last(): K | undefined;
  /** Order is part of equality (README §4). */
  equals(other: ReadonlyOrderedHashSet<K>): boolean;
}

export interface OrderedHashSet<K> extends ReadonlyOrderedHashSet<K> {
  add(key: K): OrderedHashSet<K>;
  delete(key: K): OrderedHashSet<K>;
  asMutable(): MutableOrderedHashSet<K>;
  withMutations(fn: (s: MutableOrderedHashSet<K>) => void): OrderedHashSet<K>;
}

export interface MutableOrderedHashSet<K> extends ReadonlyOrderedHashSet<K> {
  add(key: K): this;
  delete(key: K): this;
  asMutable(): this;
  /** A new persistent set; this handle is dead (README §3). */
  asImmutable(): OrderedHashSet<K>;
}

const EmptyKeys = Vector.empty<any>() as unknown as VectorImpl<any>;

export class OrderedHashSetImpl<K> {
  private _hash: number | undefined = undefined;
  private dead = false;

  constructor(
    public index: HashMapImpl<K, number>,
    public keys: VectorImpl<any>,
    public holes: number,
    public cfg: Resolved,
    public owner: Owner | null,
  ) {}

  get size() {
    return this.index.size;
  }

  has(key: K): boolean {
    return rootGet(this.index.root, key, this.index.cfg) !== NotFound;
  }

  first(): K | undefined {
    const v = this.keys;
    for (let i = 0; i < v.size; i++) {
      const k = v.get(i);
      if (k !== Hole) return k;
    }
    return undefined;
  }

  last(): K | undefined {
    // trailing holes are trimmed on delete
    return this.keys.last();
  }

  forEach<C>(fn: (this: C, value: K, key: K) => void, ctx?: C): void {
    const v = this.keys;
    for (let i = 0; i < v.size; i += 32) {
      const a = v.arrayFor(i);
      for (let j = 0; j < a.length; j++) if (a[j] !== Hole) fn.call(ctx as C, a[j], a[j]);
    }
  }

  cursor(c?: Cursor<K>): Cursor<K> {
    const cur = asCursor(c);
    cur.resetKeys(this.keys);
    return cur;
  }

  *values(): IterableIterator<K> {
    const c = this.cursor();
    while (c.next()) yield c.value;
  }

  [Symbol.iterator](): IterableIterator<K> {
    return this.values();
  }

  toArray(): K[] {
    const out: K[] = [];
    this.forEach((k) => out.push(k));
    return out;
  }

  equals(other: ReadonlyOrderedHashSet<K>): boolean {
    if ((other as unknown) === this) return true;
    if (!(other instanceof OrderedHashSetImpl) || other.size !== this.size || !sameConfig(this.cfg, other.cfg)) {
      return false;
    }
    if (this._hash !== undefined && other._hash !== undefined && this._hash !== other._hash) return false;
    const a = this.cursor();
    const b = other.cursor();
    while (a.next()) {
      b.next();
      if (!this.cfg.equals(a.value, b.value)) return false;
    }
    return true;
  }

  hashCode(): number {
    if (this._hash !== undefined) return this._hash;
    const hash = this.cfg.hash;
    let acc = this.size;
    this.forEach((k) => {
      acc = hashOrdered(acc, hash(k));
    });
    if (this.owner === null) this._hash = acc;
    return acc;
  }

  private edited(index: HashMapImpl<K, number>, keys: VectorImpl<any>, holes: number): any {
    if (this.owner === null) {
      if (index === this.index && keys === this.keys) return this;
      if (index.size === 0) return emptyOrderedSet<K>(this.cfg);
      return new OrderedHashSetImpl<K>(index, keys, holes, this.cfg, null);
    }
    this.index = index;
    this.keys = keys;
    this.holes = holes;
    return this;
  }

  add(key: K): any {
    if (this.dead) deadHandle();
    if (rootGet(this.index.root, key, this.index.cfg) !== NotFound) return this;
    return this.edited(this.index.set(key, this.keys.size), this.keys.push(key), this.holes);
  }

  delete(key: K): any {
    if (this.dead) deadHandle();
    const s = rootGet(this.index.root, key, this.index.cfg);
    if (s === NotFound) return this;
    const index = this.index.delete(key);
    let holes = this.holes;
    let keys: VectorImpl<any>;
    if (s + 1 === this.keys.size) {
      keys = this.keys.pop();
      while (keys.size > 0 && keys.last() === Hole) {
        keys = keys.pop();
        holes--;
      }
    } else {
      keys = this.keys.set(s, Hole);
      holes++;
    }
    if (index.size === 0) {
      if (this.owner === null) return emptyOrderedSet<K>(this.cfg);
      return this.edited(index, EmptyKeys.asMutable() as any, 0);
    }
    if (keys.size >= CompactMin && 2 * holes > keys.size) return this.compact(index, keys);
    return this.edited(index, keys, holes);
  }

  private compact(index: HashMapImpl<K, number>, keys: VectorImpl<any>): any {
    const ni = new HashMapImpl<K, number>(EmptyFlat, 0, index.cfg, newOwner());
    const nv = EmptyKeys.asMutable() as unknown as VectorImpl<any>;
    for (let i = 0; i < keys.size; i++) {
      const k = keys.get(i);
      if (k === Hole) continue;
      ni.set(k, nv.size);
      nv.push(k);
    }
    if (this.owner === null) return this.edited(ni.asImmutable() as any, nv.asImmutable() as any, 0);
    return this.edited(ni, nv, 0);
  }

  asMutable(): any {
    if (this.owner !== null) {
      if (this.dead) deadHandle();
      return this;
    }
    return new OrderedHashSetImpl<K>(this.index.asMutable(), this.keys.asMutable(), this.holes, this.cfg, newOwner());
  }

  asImmutable(): OrderedHashSet<K> {
    if (this.owner === null) return this as any;
    if (this.dead) deadHandle();
    this.dead = true;
    if (this.index.size === 0) return emptyOrderedSet<K>(this.cfg);
    return new OrderedHashSetImpl<K>(
      this.index.asImmutable() as any,
      this.keys.asImmutable() as any,
      this.holes,
      this.cfg,
      null,
    ) as any;
  }

  withMutations(fn: (s: MutableOrderedHashSet<K>) => void): OrderedHashSet<K> {
    const m = this.asMutable();
    fn(m);
    return m.asImmutable();
  }
}

const DefaultResolved = resolveConfig(undefined);
const EmptyOrderedSet = new OrderedHashSetImpl<any>(
  new HashMapImpl<any, number>(EmptyFlat, 0, indexConfig(DefaultResolved), null),
  EmptyKeys,
  0,
  DefaultResolved,
  null,
);

function emptyOrderedSet<K>(cfg: Resolved): OrderedHashSet<K> {
  if (cfg === DefaultResolved) return EmptyOrderedSet as any;
  return new OrderedHashSetImpl<K>(
    new HashMapImpl<K, number>(EmptyFlat, 0, indexConfig(cfg), null),
    EmptyKeys,
    0,
    cfg,
    null,
  ) as any;
}

export const OrderedHashSet = {
  empty<K>(config?: SetConfig<K>): OrderedHashSet<K> {
    return emptyOrderedSet(resolveConfig(config));
  },
  of<K>(keys: Iterable<K>, config?: SetConfig<K>): OrderedHashSet<K> {
    const m = OrderedHashSet.empty<K>(config).asMutable();
    for (const k of keys) m.add(k);
    return m.asImmutable();
  },
};
