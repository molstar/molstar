/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { utf8ReadValue } from '@molstar/core/util/utf8';
import type { ReadonlyUint8Array } from '@molstar/core/util/type-helpers';

/** Slots of the cache, a power of 2 */
const Size = 256;
/** Lookups after which a cache that mostly missed turns off */
const Trial = 4096;

/**
 * The strings of short values of a column of a text read as bytes: the last string of up to 8 bytes of each of 256
 * hash slots (up to 4 bytes keyed by the bytes as an int, 5 to 8 by their start, compared by the bytes). The repeated
 * values of a column (atom, residue and chain names, element symbols, accessions) are made once, and a reader keying
 * maps by them gets the same string, whose hash V8 has computed. flat-192-cif model stage: allocation 385 MB uncached
 * past 4 bytes, 174 MB cached. A column makes one on its first string read (`ByteStringCache.read`).
 *
 * A column of unique values (ids, numbers read as strings) only pays the hash and the miss: 1M reads of 6 to 8 bytes,
 * 24 ms uncached, 41 ms cached; repeated names 25 → 17 ms. So the cache turns off after its first 4096 lookups when
 * more than 3 in 4 missed (`off`), and reads decode from then on: unique values 28 ms, repeated names 19 ms.
 */
export interface ByteStringCache {
  /** the bytes of a value of up to 4 bytes */
  keys: Int32Array;
  /** the start of a value of 5 to 8 bytes */
  starts: Uint32Array;
  lengths: Uint8Array;
  strings: string[];
  lookups: number;
  misses: number;
  /** decided after the first `Trial` lookups: the column's values repeat too rarely */
  off: boolean;
}

export function ByteStringCache(): ByteStringCache {
  return {
    keys: new Int32Array(Size),
    starts: new Uint32Array(Size),
    lengths: new Uint8Array(Size),
    strings: new Array(Size).fill(''),
    lookups: 0,
    misses: 0,
    off: false,
  };
}

export namespace ByteStringCache {
  /** Values of up to this many bytes go through the cache; longer ones are rarely repeated. */
  export const MaxLength = 8;
  /** Columns of fewer rows have no cache: its ~4 KB would outweigh their values (a dictionary's single items). */
  export const MinRows = 64;

  /** `utf8ReadValue(data, s, l)` for `l <= MaxLength`, through the cache. */
  export function read(cache: ByteStringCache, data: ReadonlyUint8Array, s: number, l: number): string {
    if (cache.off) return utf8ReadValue(data, s, l);
    if (++cache.lookups === Trial) cache.off = cache.misses > (3 * Trial) / 4;
    if (l <= 4) {
      const key = data[s] | (data[s + 1] << 8) | (data[s + 2] << 16) | (data[s + 3] << 24);
      // the bytes past the value are its neighbours
      const k = l === 4 ? key : key & ((1 << (8 * l)) - 1);
      const slot = (k ^ (k >>> 13) ^ (k >>> 22) ^ (l << 6)) & (Size - 1);
      if (cache.keys[slot] === k && cache.lengths[slot] === l) return cache.strings[slot];
      cache.misses++;
      const ret = utf8ReadValue(data, s, l);
      cache.keys[slot] = k;
      cache.lengths[slot] = l;
      cache.strings[slot] = ret;
      return ret;
    }
    let h = l;
    for (let j = 0; j < l; j++) h = (Math.imul(h, 31) + data[s + j]) | 0;
    const slot = (h ^ (h >>> 8)) & (Size - 1);
    if (cache.lengths[slot] === l) {
      const cs = cache.starts[slot];
      let j = 0;
      while (j < l && data[cs + j] === data[s + j]) j++;
      if (j === l) return cache.strings[slot];
    }
    cache.misses++;
    const ret = utf8ReadValue(data, s, l);
    cache.starts[slot] = s;
    cache.lengths[slot] = l;
    cache.strings[slot] = ret;
    return ret;
  }
}
