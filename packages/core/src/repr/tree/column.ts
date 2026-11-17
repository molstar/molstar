/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Typed columns (_spec/README.md "Values"; format.md): real data without a per-value cost. A column is a typed array
// of one of eight element types, or a dictionary-coded category. Its identity is its element type and its bytes: two
// columns are equal when both are, and its hash is over them, computed once per column object (a column is never
// mutated once it is in a tree). A category is kept canonical (its dictionary in the order of first use, no unused
// entry), so its identity is the sequence of its values and any encoding that writes the values gives it back.

import { murmurHash3_128, murmurHash3_128_fromBytes } from '../../data/util/hash-functions.js';

export type NumericColumn =
  | Int8Array
  | Uint8Array
  | Int16Array
  | Uint16Array
  | Int32Array
  | Uint32Array
  | Float32Array
  | Float64Array;

export type NumericColumnType = 'int8' | 'uint8' | 'int16' | 'uint16' | 'int32' | 'uint32' | 'float32' | 'float64';

export type ColumnType = NumericColumnType | 'category';

export const ColumnTypes: readonly ColumnType[] = [
  'int8',
  'uint8',
  'int16',
  'uint16',
  'int32',
  'uint32',
  'float32',
  'float64',
  'category',
];

/**
 * A dictionary-coded column of strings: `codes[i]` indexes `dict`, `-1` is a missing value. Made by `categoryColumn`
 * only, so it is canonical: `dict` in the order of first use, every entry used.
 */
export class CategoryColumn {
  /** Use `categoryColumn` (it keeps the column canonical). */
  constructor(
    readonly dict: readonly string[],
    readonly codes: Int32Array,
  ) {}

  get length(): number {
    return this.codes.length;
  }

  /** The value at `i`; `null` when missing. */
  get(i: number): string | null {
    const c = this.codes[i];
    return c < 0 ? null : this.dict[c];
  }
}

export type Column = NumericColumn | CategoryColumn;

const Constructors = {
  int8: Int8Array,
  uint8: Uint8Array,
  int16: Int16Array,
  uint16: Uint16Array,
  int32: Int32Array,
  uint32: Uint32Array,
  float32: Float32Array,
  float64: Float64Array,
} as const;

/** The typed array constructor of a numeric element type. */
export function columnConstructor(type: NumericColumnType) {
  return Constructors[type];
}

/** The integer range of an integer element type; undefined for the float types. */
export function integerRange(type: NumericColumnType): readonly [number, number] | undefined {
  switch (type) {
    case 'int8':
      return [-128, 127];
    case 'uint8':
      return [0, 255];
    case 'int16':
      return [-32768, 32767];
    case 'uint16':
      return [0, 65535];
    case 'int32':
      return [-2147483648, 2147483647];
    case 'uint32':
      return [0, 4294967295];
    default:
      return undefined;
  }
}

export function isFloatType(type: ColumnType): type is 'float32' | 'float64' {
  return type === 'float32' || type === 'float64';
}

/** The element type of a column, or undefined when the value is no column (a DataView, a BigInt64Array, …). */
export function columnType(value: unknown): ColumnType | undefined {
  if (value instanceof CategoryColumn) return 'category';
  if (!ArrayBuffer.isView(value)) return undefined;
  if (value instanceof Float64Array) return 'float64';
  if (value instanceof Float32Array) return 'float32';
  if (value instanceof Int32Array) return 'int32';
  if (value instanceof Uint32Array) return 'uint32';
  if (value instanceof Int16Array) return 'int16';
  if (value instanceof Uint16Array) return 'uint16';
  if (value instanceof Int8Array) return 'int8';
  // a Uint8ClampedArray is no column: its writes clamp, so its bytes are not what was written
  if (value instanceof Uint8Array && !(value instanceof Uint8ClampedArray)) return 'uint8';
  return undefined;
}

export function isColumn(value: unknown): value is Column {
  return columnType(value) !== undefined;
}

/** A canonical category column of `values` (`null`: missing). */
export function categoryColumn(values: ArrayLike<string | null>): CategoryColumn {
  const dict: string[] = [];
  const index = new Map<string, number>();
  const codes = new Int32Array(values.length);
  for (let i = 0; i < values.length; i++) {
    const v = values[i];
    if (v === null) {
      codes[i] = -1;
      continue;
    }
    let c = index.get(v);
    if (c === undefined) {
      c = dict.length;
      dict.push(v);
      index.set(v, c);
    }
    codes[i] = c;
  }
  return new CategoryColumn(dict, codes);
}

/** A canonical category column from a dictionary and codes in any order (`-1`: missing); a code out of range throws. */
export function categoryFromCodes(dict: readonly string[], codes: ArrayLike<number>): CategoryColumn {
  const values: (string | null)[] = new Array(codes.length);
  for (let i = 0; i < codes.length; i++) {
    const c = codes[i];
    if (c === -1) values[i] = null;
    else if (Number.isInteger(c) && c >= 0 && c < dict.length) values[i] = dict[c];
    else throw new Error('A category code is out of its dictionary');
  }
  return categoryColumn(values);
}

/** A column's values as a plain array (a category's as strings and `null`). */
export function columnValues(column: Column): (number | string | null)[] {
  if (column instanceof CategoryColumn) {
    const out: (string | null)[] = new Array(column.length);
    for (let i = 0; i < column.length; i++) out[i] = column.get(i);
    return out;
  }
  return Array.from(column);
}

function bytesOf(view: ArrayBufferView): Uint8Array {
  return new Uint8Array(view.buffer, view.byteOffset, view.byteLength);
}

const Seed = 0x72657072; // 'repr'
const Hashes = new WeakMap<Column, string>();

/**
 * A column's identity as 32 hex digits: its element type and its bytes (a category: its dictionary and its codes),
 * hashed once per column object. Bytes are the platform's, little-endian on every engine that runs this.
 */
export function hashColumn(column: Column): string {
  let hash = Hashes.get(column);
  if (hash !== undefined) return hash;
  if (column instanceof CategoryColumn) {
    const codes = murmurHash3_128_fromBytes(bytesOf(column.codes), Seed);
    hash = murmurHash3_128(`category:${JSON.stringify(column.dict)}:${codes}`, Seed);
  } else {
    // the element type is in the seed, so equal bytes of two types hash apart
    const type = ColumnTypes.indexOf(columnType(column)!);
    hash = murmurHash3_128_fromBytes(bytesOf(column), Seed + type + 1);
  }
  Hashes.set(column, hash);
  return hash;
}

/** Whether two columns are equal: the same element type and the same bytes (so NaN is NaN, and -0 is not 0). */
export function columnsEqual(a: Column, b: Column): boolean {
  if (a === b) return true;
  const type = columnType(a);
  if (type !== columnType(b) || a.length !== b.length) return false;
  if (a instanceof CategoryColumn) {
    const c = b as CategoryColumn;
    if (a.dict.length !== c.dict.length) return false;
    for (let i = 0; i < a.dict.length; i++) if (a.dict[i] !== c.dict[i]) return false;
    return bytesEqual(bytesOf(a.codes), bytesOf(c.codes));
  }
  return bytesEqual(bytesOf(a), bytesOf(b as NumericColumn));
}

function bytesEqual(a: Uint8Array, b: Uint8Array): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

/**
 * A plain array as a column of `type`, or undefined when an item does not fit: an integer type takes integers in its
 * range; a float type takes numbers and `null` (missing, NaN); a category takes strings and, with `missing`, `null`.
 * A column of `type` is returned as it is; a column of another type does not fit.
 */
export function toColumn(type: ColumnType, value: unknown, missing = false): Column | undefined {
  const have = columnType(value);
  if (have !== undefined) return have === type ? (value as Column) : undefined;
  if (!Array.isArray(value)) return undefined;
  if (type === 'category') {
    for (const v of value) if (typeof v !== 'string' && !(missing && v === null)) return undefined;
    return categoryColumn(value as (string | null)[]);
  }
  const out = new Constructors[type](value.length);
  const range = integerRange(type);
  for (let i = 0; i < value.length; i++) {
    const v = value[i];
    if (range) {
      if (typeof v !== 'number' || !Number.isInteger(v) || v < range[0] || v > range[1]) return undefined;
      out[i] = v;
    } else if (v === null) out[i] = Number.NaN;
    else if (typeof v === 'number') out[i] = v;
    else return undefined;
  }
  return out;
}

/** The index of the first item of a plain array that does not fit a column of `type`, or -1 (see `toColumn`). */
export function columnMisfit(type: ColumnType, value: readonly unknown[], missing = false): number {
  const range = type === 'category' ? undefined : integerRange(type);
  for (let i = 0; i < value.length; i++) {
    const v = value[i];
    if (type === 'category') {
      if (typeof v !== 'string' && !(missing && v === null)) return i;
    } else if (range) {
      if (typeof v !== 'number' || !Number.isInteger(v) || v < range[0] || v > range[1]) return i;
    } else if (v !== null && typeof v !== 'number') return i;
  }
  return -1;
}
