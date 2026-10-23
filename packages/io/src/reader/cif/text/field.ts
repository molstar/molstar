/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, ColumnHelpers } from '@molstar/core/data/db';
import { parseFloatBytes, parseIntBytes } from '@molstar/core/util/number-parser';
import { utf8ReadValue } from '@molstar/core/util/utf8';
import type { CifField } from '../data-model.js';

/**
 * The values of a text CIF loop (or of one value) as byte ranges of the file, field by field: a start per value (Uint32)
 * and its length (Uint8, 255 for a value of 255 bytes or more, whose end is in `long`). 5 bytes a value, where a start
 * and an end were 8; one pair of arrays for the loop, so the reading loop keeps them in locals, and a field's values
 * are contiguous for its readers. Internal to the CIF text reader; other formats keep `Tokens`.
 */
export interface ByteTokens {
  readonly data: Uint8Array;
  readonly fieldCount: number;
  /** rows of room: the value of field `f` in row `r` is at `f * capacity + r` */
  capacity: number;
  /** values read, row by row */
  count: number;
  starts: Uint32Array;
  lengths: Uint8Array;
  /** end of each value of 255 bytes or more, by index */
  long: Map<number, number> | undefined;
}

export namespace ByteTokens {
  export function create(data: Uint8Array, fieldCount: number, capacity: number): ByteTokens {
    capacity = Math.max(capacity, 1);
    return {
      data,
      fieldCount,
      capacity,
      count: 0,
      starts: new Uint32Array(capacity * fieldCount),
      lengths: new Uint8Array(capacity * fieldCount),
      long: void 0,
    };
  }

  /** Adds the next value, row by row. */
  export function add(tokens: ByteTokens, start: number, end: number) {
    const { count, fieldCount } = tokens;
    const row = Math.floor(count / fieldCount);
    if (row === tokens.capacity) grow(tokens, Math.ceil(1.61 * row));
    set(tokens, (count - row * fieldCount) * tokens.capacity + row, start, end);
    tokens.count = count + 1;
  }

  export function set(tokens: ByteTokens, index: number, start: number, end: number) {
    tokens.starts[index] = start;
    const length = end - start;
    if (length < 255) {
      tokens.lengths[index] = length;
    } else {
      tokens.lengths[index] = 255;
      (tokens.long ??= new Map()).set(index, end);
    }
  }

  /** Room for `capacity` rows (at least the rows read). */
  export function grow(tokens: ByteTokens, capacity: number) {
    const { fieldCount, capacity: old } = tokens;
    const rows = Math.ceil(tokens.count / fieldCount);
    capacity = Math.max(capacity, rows);
    if (capacity === old) return;
    const starts = new Uint32Array(capacity * fieldCount);
    const lengths = new Uint8Array(capacity * fieldCount);
    for (let f = 0; f < fieldCount; f++) {
      starts.set(tokens.starts.subarray(f * old, f * old + rows), f * capacity);
      lengths.set(tokens.lengths.subarray(f * old, f * old + rows), f * capacity);
    }
    if (tokens.long) {
      const long = new Map<number, number>();
      for (const [i, end] of tokens.long) long.set(Math.floor(i / old) * capacity + (i % old), end);
      tokens.long = long;
    }
    tokens.starts = starts;
    tokens.lengths = lengths;
    tokens.capacity = capacity;
  }

  /** Drops spare room of more than an eighth. */
  export function trim(tokens: ByteTokens) {
    const rows = Math.ceil(tokens.count / tokens.fieldCount);
    if (tokens.capacity - rows > tokens.capacity / 8) grow(tokens, rows);
  }
}

/** Values of up to this many bytes go through the field's string cache. */
const CachedLength = 8;
/** Fields of fewer rows have no cache: its ~3 KB would outweigh a dictionary's single items (the CCD). */
const CachedMinRows = 64;
/** Slots of the cache, a power of 2 */
const CacheSize = 256;

interface StringCache {
  keys: Int32Array;
  lengths: Uint8Array;
  strings: string[];
}

function StringCache(): StringCache {
  return {
    keys: new Int32Array(CacheSize),
    lengths: new Uint8Array(CacheSize),
    strings: new Array(CacheSize).fill(''),
  };
}

/**
 * Field `column` of byte tokens: `CifField.ofTokens` on the bytes of the file, values decoded as UTF-8.
 *
 * `str` keeps the last string of up to 8 bytes of each of 256 hash slots (up to 4 bytes keyed by the bytes as an
 * int): the repeated values of a column (atom, residue and chain names, element symbols, accessions) are made once,
 * and a reader keying maps by them gets the same string, whose hash V8 has computed. flat-192-cif model stage:
 * allocation 385 MB uncached past 4 bytes, 174 MB cached.
 */
export function byteTokensField(tokens: ByteTokens, column: number): CifField {
  const { data, starts, lengths, long } = tokens;
  const rowCount = Math.floor(tokens.count / tokens.fieldCount);
  const offset = column * tokens.capacity;

  let cache: StringCache | undefined;
  const cached = rowCount >= CachedMinRows;

  const str: CifField['str'] = (row) => {
    const i = offset + row;
    const s = starts[i];
    const l = lengths[i];
    if (l === 1) {
      const c = data[s];
      if (c === 46 /* . */ || c === 63 /* ? */) return '';
    }
    if (l > CachedLength || !cached) return utf8ReadValue(data, s, l < 255 ? l : long!.get(i)! - s);
    if (!cache) cache = StringCache();
    if (l <= 4) {
      const key = data[s] | (data[s + 1] << 8) | (data[s + 2] << 16) | (data[s + 3] << 24);
      // the bytes past the value are its neighbours
      const k = l === 4 ? key : key & ((1 << (8 * l)) - 1);
      const slot = (k ^ (k >>> 13) ^ (k >>> 22) ^ (l << 6)) & (CacheSize - 1);
      if (cache.keys[slot] === k && cache.lengths[slot] === l) return cache.strings[slot];
      const ret = utf8ReadValue(data, s, l);
      cache.keys[slot] = k;
      cache.lengths[slot] = l;
      cache.strings[slot] = ret;
      return ret;
    }
    // 5 to 8 bytes: the key is the index of the cached value, compared by its bytes
    let h = l;
    for (let j = 0; j < l; j++) h = (Math.imul(h, 31) + data[s + j]) | 0;
    const slot = (h ^ (h >>> 8)) & (CacheSize - 1);
    const c = cache.keys[slot];
    if (cache.lengths[slot] === l && lengths[c] === l) {
      const cs = starts[c];
      let j = 0;
      while (j < l && data[cs + j] === data[s + j]) j++;
      if (j === l) return cache.strings[slot];
    }
    const ret = utf8ReadValue(data, s, l);
    cache.keys[slot] = i;
    cache.lengths[slot] = l;
    cache.strings[slot] = ret;
    return ret;
  };

  const int: CifField['int'] = (row) => {
    const i = offset + row;
    const s = starts[i];
    const l = lengths[i];
    return parseIntBytes(data, s, l < 255 ? s + l : long!.get(i)!) || 0;
  };
  const float: CifField['float'] = (row) => {
    const i = offset + row;
    const s = starts[i];
    const l = lengths[i];
    return parseFloatBytes(data, s, l < 255 ? s + l : long!.get(i)!) || 0;
  };

  const valueKind: CifField['valueKind'] = (row) => {
    const i = offset + row;
    const l = lengths[i];
    if (l > 1) return Column.ValueKinds.Present;
    if (l === 0) return Column.ValueKinds.NotPresent;
    const v = data[starts[i]];
    if (v === 46 /* . */) return Column.ValueKinds.NotPresent;
    if (v === 63 /* ? */) return Column.ValueKinds.Unknown;
    return Column.ValueKinds.Present;
  };

  const areValuesEqual: CifField['areValuesEqual'] = (rowA, rowB) => {
    const ia = offset + rowA;
    const ib = offset + rowB;
    const a = starts[ia];
    const b = starts[ib];
    let l = lengths[ia];
    if (l !== lengths[ib]) return false;
    if (l === 255) {
      l = long!.get(ia)! - a;
      if (l !== long!.get(ib)! - b) return false;
    }
    for (let i = 0; i < l; i++) {
      if (data[a + i] !== data[b + i]) return false;
    }
    return true;
  };

  return {
    binaryEncoding: void 0,
    isDefined: true,
    rowCount,
    str,
    int,
    float,
    valueKind,
    areValuesEqual,
    toStringArray: (params) => ColumnHelpers.createAndFillArray(rowCount, str, params),
    toIntArray: (params) => ColumnHelpers.createAndFillArray(rowCount, int, params),
    toFloatArray: (params) => ColumnHelpers.createAndFillArray(rowCount, float, params),
  };
}
