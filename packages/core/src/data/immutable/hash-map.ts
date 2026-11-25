/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A persistent hash map (README §2): flat up to the threshold, a CHAMP trie above it (champ.ts).

import {
  EmptyFlat,
  Edit,
  everyEntry,
  forEachEntry,
  hashEntries,
  NotFound,
  type Root,
  rootDelete,
  rootEquals,
  rootGet,
  rootSet,
} from './champ.js';
import {
  type CollectionConfig,
  deadHandle,
  type MapCursor,
  newOwner,
  type Owner,
  type Resolved,
  resolveConfig,
  sameConfig,
} from './config.js';
import { asCursor } from './cursor.js';

export interface ReadonlyHashMap<K, V> {
  readonly size: number;
  get(key: K): V | undefined;
  has(key: K): boolean;
  forEach<C>(fn: (this: C, value: V, key: K) => void, ctx?: C): void;
  /** Resets `c` (or a new cursor) to the first entry; `next()` allocates nothing (README §5). */
  cursor(c?: MapCursor<K, V>): MapCursor<K, V>;
  /** Allocates. */
  keys(): IterableIterator<K>;
  /** Allocates. */
  values(): IterableIterator<V>;
  /** Allocates. */
  toArray(): [K, V][];
  /** Structural (README §4): same config, equal keys and values. */
  equals(other: ReadonlyHashMap<K, V>): boolean;
  /** Order-independent, cached on a persistent map. */
  hashCode(): number;
  /** Allocates. */
  [Symbol.iterator](): Iterator<[K, V]>;
}

export interface HashMap<K, V> extends ReadonlyHashMap<K, V> {
  set(key: K, value: V): HashMap<K, V>;
  delete(key: K): HashMap<K, V>;
  /** Stores what `fn` returns, `undefined` included (README §7). */
  update(key: K, fn: (value: V | undefined) => V): HashMap<K, V>;
  asMutable(): MutableHashMap<K, V>;
  withMutations(fn: (m: MutableHashMap<K, V>) => void): HashMap<K, V>;
}

export interface MutableHashMap<K, V> extends ReadonlyHashMap<K, V> {
  set(key: K, value: V): this;
  delete(key: K): this;
  update(key: K, fn: (value: V | undefined) => V): this;
  asMutable(): this;
  /** A new persistent map sharing this one's nodes; this handle is dead (README §3). */
  asImmutable(): HashMap<K, V>;
}

export class HashMapImpl<K, V> {
  private _hash: number | undefined = undefined;
  private dead = false;

  constructor(
    public root: Root,
    public size: number,
    public cfg: Resolved,
    public owner: Owner | null,
  ) {}

  get(key: K): V | undefined {
    const v = rootGet(this.root, key, this.cfg);
    return v === NotFound ? undefined : v;
  }

  has(key: K): boolean {
    return rootGet(this.root, key, this.cfg) !== NotFound;
  }

  forEach<C>(fn: (this: C, value: V, key: K) => void, ctx?: C): void {
    forEachEntry(this.root, fn, ctx);
  }

  cursor(c?: MapCursor<K, V>): MapCursor<K, V> {
    const cur = asCursor(c);
    cur.resetTrie(this.root);
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
    everyEntry(this.root, (k, v) => {
      out.push([k, v]);
      return true;
    });
    return out;
  }

  equals(other: ReadonlyHashMap<K, V>): boolean {
    if ((other as unknown) === this) return true;
    if (!(other instanceof HashMapImpl) || other.size !== this.size || !sameConfig(this.cfg, other.cfg)) return false;
    if (this._hash !== undefined && other._hash !== undefined && this._hash !== other._hash) return false;
    return rootEquals(this.root, other.root, this.cfg, true);
  }

  hashCode(): number {
    if (this._hash !== undefined) return this._hash;
    const h = (hashEntries(this.root, this.cfg, true) + Math.imul(this.size, 0x9e3779b1)) | 0;
    if (this.owner === null) this._hash = h;
    return h;
  }

  private edited(root: Root, size: number): any {
    if (this.owner === null) {
      return root === this.root && size === this.size ? this : new HashMapImpl<K, V>(root, size, this.cfg, null);
    }
    this.root = root;
    this.size = size;
    return this;
  }

  set(key: K, value: V): any {
    if (this.dead) deadHandle();
    const root = rootSet(this.root, this.size, key, value, this.cfg, this.owner);
    return this.edited(root, Edit.added ? this.size + 1 : this.size);
  }

  delete(key: K): any {
    if (this.dead) deadHandle();
    const root = rootDelete(this.root, this.size, key, this.cfg, this.owner);
    return this.edited(root, Edit.removed ? this.size - 1 : this.size);
  }

  update(key: K, fn: (value: V | undefined) => V): any {
    return this.set(key, fn(this.get(key)));
  }

  asMutable(): any {
    if (this.owner !== null) {
      if (this.dead) deadHandle();
      return this;
    }
    return new HashMapImpl<K, V>(this.root, this.size, this.cfg, newOwner());
  }

  asImmutable(): HashMap<K, V> {
    if (this.owner === null) return this as any;
    if (this.dead) deadHandle();
    this.dead = true;
    return new HashMapImpl<K, V>(this.root, this.size, this.cfg, null) as any;
  }

  withMutations(fn: (m: MutableHashMap<K, V>) => void): HashMap<K, V> {
    const m = this.asMutable();
    fn(m);
    return m.asImmutable();
  }
}

const EmptyMap = new HashMapImpl<any, any>(EmptyFlat, 0, resolveConfig(undefined), null);

export const HashMap = {
  empty<K, V>(config?: CollectionConfig<K, V>): HashMap<K, V> {
    const cfg = resolveConfig(config);
    if (cfg === EmptyMap.cfg) return EmptyMap as any;
    return new HashMapImpl<K, V>(EmptyFlat, 0, cfg, null) as any;
  },
  of<K, V>(entries: Iterable<readonly [K, V]>, config?: CollectionConfig<K, V>): HashMap<K, V> {
    const m = HashMap.empty<K, V>(config).asMutable();
    for (const [k, v] of entries) m.set(k, v);
    return m.asImmutable();
  },
};
