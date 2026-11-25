/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Keys (_spec/README.md "Queries", Keys; the result type was *Keys*, renamed by the author, 2026-11-20): things
// data makes (an annotation's items, an alignment's rows, a structure's groups) are named by their owner's address and
// their key in it. A key is the data's own name, a string or a number, and the two never match (`1` is not `"1"`), as
// in an address. Keys are a value: sorted, each once, compared and hashed.

import { hash2, hashString } from '@molstar/core/data/util/hash-functions';
import { compareStrings } from './compare.js';

export interface Key {
  /** The owner's address, printed canonically (`/data/p53/domains`). */
  readonly owner: string;
  readonly key: string | number;
}

export interface Keys {
  /** Sorted by `compareKeys`, each once. */
  readonly keys: readonly Key[];
}

export const EmptyKeys: Keys = Object.freeze({ keys: Object.freeze([]) });

/** By owner, then numbers before strings, numbers by value, strings by code unit. */
export function compareKeys(a: Key, b: Key): number {
  const o = compareStrings(a.owner, b.owner);
  if (o) return o;
  const x = a.key;
  const y = b.key;
  if (typeof x === 'number')
    return typeof y === 'number' ? (x < y ? -1 : x > y ? 1 : Object.is(x, y) ? 0 : compareStrings(`${x}`, `${y}`)) : -1;
  return typeof y === 'number' ? 1 : compareStrings(x, y);
}

/** Normalized, frozen Keys. */
export function keySet(keys: Iterable<Key>): Keys {
  const sorted = [...keys].sort(compareKeys);
  const out: Key[] = [];
  for (const k of sorted)
    if (!out.length || compareKeys(out[out.length - 1], k) !== 0) out.push(Object.freeze({ ...k }));
  return out.length ? Object.freeze({ keys: Object.freeze(out) }) : EmptyKeys;
}

export function keysEqual(a: Keys, b: Keys): boolean {
  if (a === b) return true;
  if (a.keys.length !== b.keys.length) return false;
  for (let i = 0; i < a.keys.length; i++) if (compareKeys(a.keys[i], b.keys[i]) !== 0) return false;
  return true;
}

const hashes = new WeakMap<Keys, number>();

/** A 32-bit hash, equal for equal Keys; cached per value. */
export function hashKeys(e: Keys): number {
  let h = hashes.get(e);
  if (h !== undefined) return h;
  h = e.keys.length;
  for (const k of e.keys) {
    const key = typeof k.key === 'number' ? hashString(`#${k.key}`) : hashString(k.key);
    h = hash2(hash2(h, hashString(k.owner)), key);
  }
  hashes.set(e, h);
  return h;
}

/** Merges two sorted key lists by a predicate on (in a, in b). */
function combine(a: Keys, b: Keys, keep: (inA: boolean, inB: boolean) => boolean): Keys {
  const out: Key[] = [];
  let i = 0;
  let j = 0;
  while (i < a.keys.length || j < b.keys.length) {
    const c = i >= a.keys.length ? 1 : j >= b.keys.length ? -1 : compareKeys(a.keys[i], b.keys[j]);
    const k = c <= 0 ? a.keys[i] : b.keys[j];
    if (keep(c <= 0, c >= 0)) out.push(k);
    if (c <= 0) i++;
    if (c >= 0) j++;
  }
  return out.length ? Object.freeze({ keys: Object.freeze(out) }) : EmptyKeys;
}

export const unionKeys = (a: Keys, b: Keys) => combine(a, b, (x, y) => x || y);
export const intersectKeys = (a: Keys, b: Keys) => combine(a, b, (x, y) => x && y);
export const subtractKeys = (a: Keys, b: Keys) => combine(a, b, (x, y) => x && !y);
