/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { murmurHash3_128 } from '@molstar/core/data/util/hash-functions';

/**
 * A 128-bit hash of `JSON.stringify(value)`, built from the hashes of its parts (`load.md §3`): equal JSON, equal hash;
 * a part that is an object already in `cache` is not walked again, so values that share most of their parts with
 * earlier ones (animation frames) cost only what is new. Not the same value as `hashMurmur128o(value)`.
 *
 * The cache is keyed by object: valid while the objects in it are not changed (one load).
 */
export function jsonHash(value: unknown, cache: WeakMap<object, string>): string {
  // (small values: the native JSON text is faster than the parts; which way is decided by the JSON content alone, so
  // equal JSON still gives equal hashes)
  if (isSmallJson(value, SmallJsonEntries)) return murmurHash3_128(JSON.stringify(value) ?? 'undefined', 42);
  return murmurHash3_128(partOf(value, cache) ?? 'undefined', 42);
}

/** Values with at most this many JSON entries (keys and array items, at any depth) are hashed from their JSON text. */
const SmallJsonEntries = 256;

/** At most `budget` JSON-visible entries; stops counting as soon as it is exceeded. */
function isSmallJson(value: unknown, budget: number): boolean {
  let left = budget;
  const count = (v: unknown): boolean => {
    if (typeof v !== 'object' || v === null) return true;
    if (typeof (v as any).toJSON === 'function') return true;
    if (Array.isArray(v)) {
      left -= v.length;
      if (left < 0) return false;
      for (const x of v) if (!count(x)) return false;
      return true;
    }
    for (const k in v) {
      if (!Object.hasOwn(v, k)) continue;
      const x = (v as any)[k];
      if (x === undefined || typeof x === 'function' || typeof x === 'symbol') continue;
      if (--left < 0 || !count(x)) return false;
    }
    return true;
  };
  return count(value);
}

/**
 * Parts up to this length stay JSON text in their parent's text (hashing many small parts would cost more than the
 * text); longer ones become `#` + their hash, kept in the cache (no JSON text starts with `#`).
 */
const MaxLiteralLength = 1024;

/** The JSON text of `value`, or the hash of it; `undefined` when JSON leaves it out. */
function partOf(value: unknown, cache: WeakMap<object, string>): string | undefined {
  switch (typeof value) {
    case 'number':
      return Number.isFinite(value) ? String(value) : 'null';
    case 'string':
      return JSON.stringify(value);
    case 'boolean':
      return value ? 'true' : 'false';
    case 'undefined':
    case 'function':
    case 'symbol':
      return undefined;
    case 'bigint':
      throw new TypeError('Do not know how to serialize a BigInt');
  }
  if (value === null) return 'null';
  const obj = value as Record<string, unknown>;
  if (typeof obj.toJSON === 'function') return partOf(obj.toJSON(''), cache);
  const cached = cache.get(obj);
  if (cached !== undefined) return cached;

  let text: string;
  if (Array.isArray(obj)) {
    text = '[';
    for (let i = 0; i < obj.length; i++) {
      if (i > 0) text += ',';
      text += partOf(obj[i], cache) ?? 'null';
    }
    text += ']';
  } else {
    text = '{';
    let first = true;
    for (const key of Object.keys(obj)) {
      const part = partOf(obj[key], cache);
      if (part === undefined) continue;
      if (!first) text += ',';
      first = false;
      text += JSON.stringify(key) + ':' + part;
    }
    text += '}';
  }
  if (text.length <= MaxLiteralLength) return text;
  const result = '#' + murmurHash3_128(text, 42);
  cache.set(obj, result);
  return result;
}
