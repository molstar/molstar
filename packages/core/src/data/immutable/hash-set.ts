/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A persistent hash set (README §2): the map's trie with every key stored as its own value.

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
  type Cursor,
  deadHandle,
  newOwner,
  type Owner,
  type Resolved,
  resolveConfig,
  sameConfig,
} from './config.js';
import { asCursor } from './cursor.js';

export type SetConfig<K> = Pick<CollectionConfig<K>, 'hash' | 'equals'>;

export interface ReadonlyHashSet<K> {
  readonly size: number;
  has(key: K): boolean;
  forEach<C>(fn: (this: C, value: K, key: K) => void, ctx?: C): void;
  /** Resets `c` (or a new cursor) to the first key; `next()` allocates nothing (README §5). */
  cursor(c?: Cursor<K>): Cursor<K>;
  /** Allocates. */
  values(): IterableIterator<K>;
  /** Allocates. */
  toArray(): K[];
  equals(other: ReadonlyHashSet<K>): boolean;
  hashCode(): number;
  /** Allocates. */
  [Symbol.iterator](): Iterator<K>;
}

export interface HashSet<K> extends ReadonlyHashSet<K> {
  add(key: K): HashSet<K>;
  delete(key: K): HashSet<K>;
  asMutable(): MutableHashSet<K>;
  withMutations(fn: (s: MutableHashSet<K>) => void): HashSet<K>;
}

export interface MutableHashSet<K> extends ReadonlyHashSet<K> {
  add(key: K): this;
  delete(key: K): this;
  asMutable(): this;
  /** A new persistent set; this handle is dead (README §3). */
  asImmutable(): HashSet<K>;
}

const alwaysEqual = () => true;
const setConfigs = new WeakMap<Resolved, Resolved>();

/** The set's own config for the trie: a key's stored value is itself, so a present key is never rewritten. */
export function setConfig(cfg: Resolved): Resolved {
  let r = setConfigs.get(cfg);
  if (!r) {
    r = { ...cfg, valueEquals: alwaysEqual, valueHash: undefined };
    setConfigs.set(cfg, r);
  }
  return r;
}

export class HashSetImpl<K> {
  private _hash: number | undefined = undefined;
  private dead = false;

  constructor(
    public root: Root,
    public size: number,
    public cfg: Resolved,
    public owner: Owner | null,
    /** the trie's config (`setConfig`), passed on by edits: a `WeakMap` lookup an edit otherwise */
    private readonly trie: Resolved = setConfig(cfg),
  ) {}

  has(key: K): boolean {
    return rootGet(this.root, key, this.trie) !== NotFound;
  }

  forEach<C>(fn: (this: C, value: K, key: K) => void, ctx?: C): void {
    forEachEntry(this.root, fn, ctx);
  }

  cursor(c?: Cursor<K>): Cursor<K> {
    const cur = asCursor(c);
    cur.resetTrie(this.root);
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
    everyEntry(this.root, (k) => {
      out.push(k);
      return true;
    });
    return out;
  }

  equals(other: ReadonlyHashSet<K>): boolean {
    if ((other as unknown) === this) return true;
    if (!(other instanceof HashSetImpl) || other.size !== this.size || !sameConfig(this.cfg, other.cfg)) return false;
    if (this._hash !== undefined && other._hash !== undefined && this._hash !== other._hash) return false;
    return rootEquals(this.root, other.root, this.trie, false);
  }

  hashCode(): number {
    if (this._hash !== undefined) return this._hash;
    const h = (hashEntries(this.root, this.cfg, false) + Math.imul(this.size, 0x9e3779b1)) | 0;
    if (this.owner === null) this._hash = h;
    return h;
  }

  private edited(root: Root, size: number): any {
    if (this.owner === null) {
      return root === this.root && size === this.size
        ? this
        : new HashSetImpl<K>(root, size, this.cfg, null, this.trie);
    }
    this.root = root;
    this.size = size;
    return this;
  }

  add(key: K): any {
    if (this.dead) deadHandle();
    const root = rootSet(this.root, this.size, key, key, this.trie, this.owner);
    return this.edited(root, Edit.added ? this.size + 1 : this.size);
  }

  delete(key: K): any {
    if (this.dead) deadHandle();
    const root = rootDelete(this.root, this.size, key, this.trie, this.owner);
    return this.edited(root, Edit.removed ? this.size - 1 : this.size);
  }

  asMutable(): any {
    if (this.owner !== null) {
      if (this.dead) deadHandle();
      return this;
    }
    return new HashSetImpl<K>(this.root, this.size, this.cfg, newOwner(), this.trie);
  }

  asImmutable(): HashSet<K> {
    if (this.owner === null) return this as any;
    if (this.dead) deadHandle();
    this.dead = true;
    return new HashSetImpl<K>(this.root, this.size, this.cfg, null, this.trie) as any;
  }

  withMutations(fn: (s: MutableHashSet<K>) => void): HashSet<K> {
    const m = this.asMutable();
    fn(m);
    return m.asImmutable();
  }
}

const EmptySet = new HashSetImpl<any>(EmptyFlat, 0, resolveConfig(undefined), null);

export const HashSet = {
  empty<K>(config?: SetConfig<K>): HashSet<K> {
    const cfg = resolveConfig(config);
    if (cfg === EmptySet.cfg) return EmptySet as any;
    return new HashSetImpl<K>(EmptyFlat, 0, cfg, null) as any;
  },
  of<K>(keys: Iterable<K>, config?: SetConfig<K>): HashSet<K> {
    const m = HashSet.empty<K>(config).asMutable();
    for (const k of keys) m.add(k);
    return m.asImmutable();
  },
};
