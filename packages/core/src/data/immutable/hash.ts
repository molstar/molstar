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

/** The float64's bits, both words through `fmix32`; −0 hashes as 0, every NaN alike (README §4). */
export function hashNumber(n: number): number {
  if (Number.isNaN(n)) return HashNaN;
  if (n === 0) return 0;
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
