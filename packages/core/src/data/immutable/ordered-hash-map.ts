/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A persistent map in insertion order (README §2): a `HashMap` from key to slot, and a `Vector` of keys and values
// interleaved (slot `s` at `2s` and `2s + 1`). A removal leaves a hole; the vector is compacted when holes pass half
// its slots and it holds at least `CompactMin` slots. Setting a present key keeps its place.

import { EmptyFlat, NotFound, rootGet } from './champ.js';
import {
  type CollectionConfig,
  deadHandle,
  hashEntry,
  hashOrdered,
  type MapCursor,
  newOwner,
  type Owner,
  type Resolved,
  resolveConfig,
  sameConfig,
} from './config.js';
import { asCursor, Hole } from './cursor.js';
import { HashMapImpl, type ReadonlyHashMap } from './hash-map.js';
import { Vector, VectorImpl } from './vector.js';

/** Compaction waits for this many slots (README §2, as immutable.js's `OrderedMap`). */
export const CompactMin = 32;

export interface ReadonlyOrderedHashMap<K, V> extends Omit<ReadonlyHashMap<K, V>, 'equals'> {
  /** The first entry's value in insertion order. */
  first(): V | undefined;
  /** The last entry's value in insertion order. */
  last(): V | undefined;
  /** Order is part of equality (README §4). */
  equals(other: ReadonlyOrderedHashMap<K, V>): boolean;
}

export interface OrderedHashMap<K, V> extends ReadonlyOrderedHashMap<K, V> {
  set(key: K, value: V): OrderedHashMap<K, V>;
  delete(key: K): OrderedHashMap<K, V>;
  update(key: K, fn: (value: V | undefined) => V): OrderedHashMap<K, V>;
  asMutable(): MutableOrderedHashMap<K, V>;
  withMutations(fn: (m: MutableOrderedHashMap<K, V>) => void): OrderedHashMap<K, V>;
}

export interface MutableOrderedHashMap<K, V> extends ReadonlyOrderedHashMap<K, V> {
  set(key: K, value: V): this;
  delete(key: K): this;
  update(key: K, fn: (value: V | undefined) => V): this;
  asMutable(): this;
  /** A new persistent map; this handle is dead (README §3). */
  asImmutable(): OrderedHashMap<K, V>;
}

const indexConfigs = new WeakMap<Resolved, Resolved>();

/** The slot index's config: the keys' hash and equality, slots compared by identity. */
export function indexConfig(cfg: Resolved): Resolved {
  let r = indexConfigs.get(cfg);
  if (!r) {
    r = { ...cfg, valueEquals: Object.is, valueHash: undefined };
    indexConfigs.set(cfg, r);
  }
  return r;
}

const EmptySlots = Vector.empty<any>() as unknown as VectorImpl<any>;

export class OrderedHashMapImpl<K, V> {
  private _hash: number | undefined = undefined;
  private dead = false;

  constructor(
    public index: HashMapImpl<K, number>,
    public slots: VectorImpl<any>,
    public holes: number,
    public cfg: Resolved,
    public owner: Owner | null,
  ) {}

  get size() {
    return this.index.size;
  }

  get(key: K): V | undefined {
    const s = rootGet(this.index.root, key, this.index.cfg);
    return s === NotFound ? undefined : this.slots.get(2 * s + 1);
  }

  has(key: K): boolean {
    return rootGet(this.index.root, key, this.index.cfg) !== NotFound;
  }

  first(): V | undefined {
    const v = this.slots;
    for (let i = 0; i < v.size; i += 2) if (v.get(i) !== Hole) return v.get(i + 1);
    return undefined;
  }

  last(): V | undefined {
    // trailing holes are trimmed on delete
    return this.slots.size > 0 ? this.slots.get(this.slots.size - 1) : undefined;
  }

  forEach<C>(fn: (this: C, value: V, key: K) => void, ctx?: C): void {
    const v = this.slots;
    for (let i = 0; i < v.size; i += 32) {
      const a = v.arrayFor(i);
      for (let j = 0; j < a.length; j += 2) if (a[j] !== Hole) fn.call(ctx as C, a[j + 1], a[j]);
    }
  }

  cursor(c?: MapCursor<K, V>): MapCursor<K, V> {
    const cur = asCursor(c);
    cur.resetPairs(this.slots);
    return cur;
  }

  *keys(): IterableIterator<K> {
    const c = this.cursor();
    while (c.next()) yield c.key;
  }

  *values(): IterableIterator<V> {
    const c = this.cursor();
    while (c.next()) yield c.value;
  }

  *[Symbol.iterator](): IterableIterator<[K, V]> {
    const c = this.cursor();
    while (c.next()) yield [c.key, c.value];
  }

  toArray(): [K, V][] {
    const out: [K, V][] = [];
    this.forEach((v, k) => out.push([k, v]));
    return out;
  }

  equals(other: ReadonlyOrderedHashMap<K, V>): boolean {
    if ((other as unknown) === this) return true;
    if (!(other instanceof OrderedHashMapImpl) || other.size !== this.size || !sameConfig(this.cfg, other.cfg)) {
      return false;
    }
    if (this._hash !== undefined && other._hash !== undefined && this._hash !== other._hash) return false;
    const a = this.cursor();
    const b = other.cursor();
    const cfg = this.cfg;
    while (a.next()) {
      b.next();
      if (!cfg.equals(a.key, b.key) || !cfg.valueEquals(a.value, b.value)) return false;
    }
    return true;
  }

  hashCode(): number {
    if (this._hash !== undefined) return this._hash;
    const cfg = this.cfg;
    const vh = cfg.valueHash;
    let acc = this.size;
    this.forEach((v, k) => {
      acc = hashOrdered(acc, vh ? hashEntry(cfg.hash(k), vh(v)) : cfg.hash(k));
    });
    if (this.owner === null) this._hash = acc;
    return acc;
  }

  private edited(index: HashMapImpl<K, number>, slots: VectorImpl<any>, holes: number): any {
    if (this.owner === null) {
      if (index === this.index && slots === this.slots) return this;
      if (index.size === 0) return emptyOrdered<K, V>(this.cfg);
      return new OrderedHashMapImpl<K, V>(index, slots, holes, this.cfg, null);
    }
    this.index = index;
    this.slots = slots;
    this.holes = holes;
    return this;
  }

  set(key: K, value: V): any {
    if (this.dead) deadHandle();
    const s = rootGet(this.index.root, key, this.index.cfg);
    if (s !== NotFound) {
      if (this.cfg.valueEquals(this.slots.get(2 * s + 1), value)) return this;
      return this.edited(this.index, this.slots.set(2 * s + 1, value), this.holes);
    }
    const index = this.index.set(key, this.slots.size >>> 1);
    return this.edited(index, this.slots.push(key).push(value), this.holes);
  }

  delete(key: K): any {
    if (this.dead) deadHandle();
    const s = rootGet(this.index.root, key, this.index.cfg);
    if (s === NotFound) return this;
    const index = this.index.delete(key);
    let holes = this.holes;
    let slots: VectorImpl<any>;
    if (2 * s + 2 === this.slots.size) {
      slots = this.slots.pop().pop();
      // trailing holes go with it, so `last()` reads the end
      while (slots.size > 0 && slots.get(slots.size - 2) === Hole) {
        slots = slots.pop().pop();
        holes--;
      }
    } else {
      slots = this.slots.set(2 * s, Hole).set(2 * s + 1, Hole);
      holes++;
    }
    if (index.size === 0) {
      if (this.owner === null) return emptyOrdered<K, V>(this.cfg);
      return this.edited(index, EmptySlots.asMutable() as any, 0);
    }
    const n = slots.size >>> 1;
    if (n >= CompactMin && 2 * holes > n) return this.compact(index, slots);
    return this.edited(index, slots, holes);
  }

  private compact(index: HashMapImpl<K, number>, slots: VectorImpl<any>): any {
    const owner = this.owner;
    const ni = new HashMapImpl<K, number>(EmptyFlat, 0, index.cfg, newOwner());
    const nv = EmptySlots.asMutable() as unknown as VectorImpl<any>;
    for (let i = 0; i < slots.size; i += 2) {
      const k = slots.get(i);
      if (k === Hole) continue;
      ni.set(k, nv.size >>> 1);
      nv.push(k).push(slots.get(i + 1));
    }
    if (owner === null) return this.edited(ni.asImmutable() as any, nv.asImmutable() as any, 0);
    return this.edited(ni, nv, 0);
  }

  update(key: K, fn: (value: V | undefined) => V): any {
    return this.set(key, fn(this.get(key)));
  }

  asMutable(): any {
    if (this.owner !== null) {
      if (this.dead) deadHandle();
      return this;
    }
    return new OrderedHashMapImpl<K, V>(
      this.index.asMutable(),
      this.slots.asMutable(),
      this.holes,
      this.cfg,
      newOwner(),
    );
  }

  asImmutable(): OrderedHashMap<K, V> {
    if (this.owner === null) return this as any;
    if (this.dead) deadHandle();
    this.dead = true;
    if (this.index.size === 0) return emptyOrdered<K, V>(this.cfg);
    return new OrderedHashMapImpl<K, V>(
      this.index.asImmutable() as any,
      this.slots.asImmutable() as any,
      this.holes,
      this.cfg,
      null,
    ) as any;
  }

  withMutations(fn: (m: MutableOrderedHashMap<K, V>) => void): OrderedHashMap<K, V> {
    const m = this.asMutable();
    fn(m);
    return m.asImmutable();
  }
}

function emptyIndex<K>(cfg: Resolved) {
  return new HashMapImpl<K, number>(EmptyFlat, 0, indexConfig(cfg), null);
}

const EmptyOrdered = new OrderedHashMapImpl<any, any>(
  emptyIndex(resolveConfig(undefined)),
  EmptySlots,
  0,
  resolveConfig(undefined),
  null,
);

function emptyOrdered<K, V>(cfg: Resolved): OrderedHashMap<K, V> {
  if (cfg === EmptyOrdered.cfg) return EmptyOrdered as any;
  return new OrderedHashMapImpl<K, V>(emptyIndex(cfg), EmptySlots, 0, cfg, null) as any;
}

export const OrderedHashMap = {
  empty<K, V>(config?: CollectionConfig<K, V>): OrderedHashMap<K, V> {
    return emptyOrdered(resolveConfig(config));
  },
  of<K, V>(entries: Iterable<readonly [K, V]>, config?: CollectionConfig<K, V>): OrderedHashMap<K, V> {
    const m = OrderedHashMap.empty<K, V>(config).asMutable();
    for (const [k, v] of entries) m.set(k, v);
    return m.asImmutable();
  },
};
