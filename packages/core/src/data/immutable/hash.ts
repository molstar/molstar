/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hashString } from '../util/hash-functions.js';

/** MurmurHash3's 32-bit finalizer: every input bit reaches the low bits the trie reads first (README §4). */
export function fmix32(h: number): number {
  h ^= h >>> 16;
  h = Math.imul(h, 0x85ebca6b);
  h ^= h >>> 13;
  h = Math.imul(h, 0xc2b2ae35);
  h ^= h >>> 16;
  return h | 0;
}

const HashTrue = 0x2c1b3c6d;
const HashFalse = 0x297a2d39;
const HashNull = 0x5bd1e995 | 0;
const HashUndefined = 0x1b873593;
const HashNaN = 0x6b43a9b5;

const f64 = new Float64Array(1);
const i32 = new Int32Array(f64.buffer);

/**
 * An int32 by its value through `fmix32` (−0 as 0); any other number by its float64 bits, both words through `fmix32`;
 * every NaN alike (README §4).
 */
export function hashNumber(n: number): number {
  if ((n | 0) === n) return fmix32(n);
  if (Number.isNaN(n)) return HashNaN;
  f64[0] = n;
  return fmix32(i32[0] ^ fmix32(i32[1]));
}

let nextId = 1;
const weakIds = new WeakMap<WeakKey, number>();
/** Registered symbols cannot be `WeakMap` keys. */
const strongIds = new Map<symbol, number>();

function canBeWeak(key: object | symbol): boolean {
  return typeof key !== 'symbol' || Symbol.keyFor(key) === undefined;
}

function identityHash(key: object | symbol, assign: boolean): number {
  const weak = canBeWeak(key);
  const id = weak ? weakIds.get(key as WeakKey) : strongIds.get(key as symbol);
  if (id !== undefined) return id;
  if (!assign) return NaN;
  const h = fmix32(nextId++);
  if (weak) weakIds.set(key as WeakKey, h);
  else strongIds.set(key as symbol, h);
  return h;
}

/**
 * The default key hash (README §4): strings by core's `hashString` and `fmix32`, numbers by their bits, bigints by
 * their decimal string, fixed constants for booleans, `null` and `undefined`, and an identity id for symbols, objects
 * and functions, assigned on the first call.
 */
export function defaultHash(key: unknown): number {
  switch (typeof key) {
    case 'string':
      return fmix32(hashString(key));
    case 'number':
      return hashNumber(key);
    case 'boolean':
      return key ? HashTrue : HashFalse;
    case 'undefined':
      return HashUndefined;
    case 'bigint':
      return fmix32(hashString(key.toString()) ^ 0x3c6ef372);
    case 'object':
      if (key === null) return HashNull;
      return identityHash(key, true);
    default:
      return identityHash(key as object | symbol, true);
  }
}

/**
 * `defaultHash` for a lookup: an object or symbol without an id yet gives NaN (absent), and none is assigned, so a read
 * allocates nothing (README §5).
 */
export function lookupDefaultHash(key: unknown): number {
  const t = typeof key;
  if (t === 'string') return fmix32(hashString(key as string));
  if (t === 'number') return hashNumber(key as number);
  if (t === 'object' ? key !== null : t === 'function' || t === 'symbol') {
    return identityHash(key as object | symbol, false);
  }
  return defaultHash(key);
}

/** A hash's own lookup form (README §4): `resolveConfig` uses it as the config's `lookupHash`. */
export const LookupHash: unique symbol = Symbol('LookupHash') as any;

export interface MemoHash {
  (key: unknown): number;
  readonly [LookupHash]: (key: unknown) => number;
}

/**
 * An opt-in string-hash memo (README §4): `defaultHash`'s values, a string's hash kept in a bounded memo of its own.
 * Never the default: pass it as a config's `hash`, one memo shared by the collections that should compare equal (a
 * config is the same only with the same function). Holds at most `capacity` strings (default 4,096): two generations
 * of `capacity / 2`, the older dropped when the newer fills, a hit in the older moved up. A miss allocates (the memo's
 * entry); a hit allocates nothing. Pays off for strings hashed again and again (ids of a long-lived index); a string
 * hashed once only costs more.
 */
export function memoHash(capacity = 4096): MemoHash {
  const half = Math.max(1, capacity >> 1);
  let current = new Map<string, number>();
  let previous = new Map<string, number>();
  const stringHash = (s: string): number => {
    let h = current.get(s);
    if (h !== undefined) return h;
    h = previous.get(s);
    if (h === undefined) h = fmix32(hashString(s));
    if (current.size >= half) {
      previous = current;
      current = new Map();
    }
    current.set(s, h);
    return h;
  };
  const hash = ((key: unknown) => (typeof key === 'string' ? stringHash(key) : defaultHash(key))) as MemoHash;
  (hash as any)[LookupHash] = (key: unknown) => (typeof key === 'string' ? stringHash(key) : lookupDefaultHash(key));
  return hash;
}
