/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Column } from '../column.js';

// A dictionary column (`Column.ofDictionary`, `storage().kind === 'dictionary'`) reads as the array column of the same
// values, through views and windows, which keep the dictionary and gather or slice its indices.

const values = ['', 'ALA', 'GLY', 'ala', 'GLY'];
const indices = Uint8Array.of(1, 2, 0, 3, 4, 1, 2, 2);
const valueKinds = [0, 0, 1, 0, 0, 0, 2, 0] as Column.ValueKind[];

function columns(schema: Column.Schema.Str) {
  const valueKind = (row: number) => valueKinds[row];
  return {
    dictionary: Column.ofDictionary({ indices, values, schema, valueKind }),
    array: Column.ofArray({ array: Array.from(indices, (i) => values[i]), schema, valueKind }),
  };
}

function read(c: Column<string>) {
  const rows = Array.from({ length: c.rowCount }, (_, i) => i);
  return {
    rowCount: c.rowCount,
    values: rows.map((i) => c.value(i)),
    kinds: rows.map((i) => c.valueKind(i)),
    equal: rows.map((i) => rows.map((j) => c.areValuesEqual(i, j))),
    array: Array.from(c.toArray()),
    window: Array.from(c.toArray({ start: 2, end: 5 })),
  };
}

describe('dictionary column', () => {
  it('reads as the array column of its values, transforms included', () => {
    for (const schema of [Column.Schema.str, Column.Schema.lstr, Column.Schema.ustr]) {
      const { dictionary, array } = columns(schema);
      expect(read(dictionary)).toEqual(read(array));
      expect(dictionary.storage!()).toEqual({ kind: 'dictionary', indices, values });
    }
  });

  it('views and windows keep the dictionary', () => {
    const { dictionary, array } = columns(Column.Schema.str);
    const map = [7, 0, 3, 3, 6];
    const view = Column.view(dictionary, map);
    expect(read(view)).toEqual(read(Column.view(array, map)));
    const vs = view.storage!();
    expect(vs.kind === 'dictionary' && vs.values).toBe(values);
    expect(vs.kind === 'dictionary' && vs.indices).toBeInstanceOf(Uint8Array);

    const window = Column.window(dictionary, 2, 6);
    expect(read(window)).toEqual(read(Column.window(array, 2, 6)));
    const ws = window.storage!();
    expect(ws.kind === 'dictionary' && (ws.indices as Uint8Array).buffer).toBe(indices.buffer);
  });

  it('compares as its values, repeated values included', () => {
    const { dictionary, array } = columns(Column.Schema.str);
    // rows 1 and 4 are the two GLY entries of the dictionary
    expect(dictionary.areValuesEqual(1, 4)).toBe(true);
    const a = Column.view(dictionary, [1, 2, 3]);
    const b = Column.view(dictionary, [4, 2, 3]);
    expect(Column.areEqual(a, b)).toBe(true);
    expect(Column.areEqual(a, Column.view(dictionary, [1, 2, 4]))).toBe(false);
    expect(Column.areEqual(a, Column.view(array, [1, 2, 3]))).toBe(true);
    expect(Column.asArrayColumn(dictionary)).toBe(dictionary);
  });

  it('a column that reads its dictionary otherwise is viewed by its values', () => {
    const digits = Column.ofDictionary({
      indices: Uint8Array.of(1, 2, 1),
      values: ['', '10', '20'],
      schema: Column.Schema.str,
    });
    // an int column over the strings, as a CIF schema makes over a string field
    const ints: Column<number> = {
      ...digits,
      schema: Column.Schema.int,
      value: (row) => +digits.value(row),
      toArray: (params) => Array.from(digits.toArray(params), (v) => +v),
    } as Column<number>;
    const view = Column.view(ints, [2, 1]);
    expect([view.value(0), view.value(1)]).toEqual([10, 20]);
    expect(Column.window(ints, 1, 3).value(0)).toBe(20);
  });
});

describe('dictionary column edges', () => {
  it('reads alike with Uint16, Int32 and plain array indices', () => {
    const many = Array.from({ length: 300 }, (_, i) => `v${i % 290}`);
    const rows = Array.from({ length: 50 }, (_, i) => (i * 37) % 300);
    for (const make of [
      (xs: number[]) => Uint16Array.from(xs),
      (xs: number[]) => Int32Array.from(xs),
      (xs: number[]) => [...xs],
    ]) {
      const dictionary = Column.ofDictionary({ indices: make(rows), values: many, schema: Column.Schema.str });
      const array = Column.ofArray({ array: rows.map((i) => many[i]), schema: Column.Schema.str });
      expect(read(dictionary)).toEqual(read(array));
      expect(read(Column.view(dictionary, [3, 1, 40]))).toEqual(read(Column.view(array, [3, 1, 40])));
      expect(read(Column.window(dictionary, 10, 20))).toEqual(read(Column.window(array, 10, 20)));
      expect(Column.areEqual(Column.window(dictionary, 0, 10), Column.window(array, 0, 10))).toBe(true);
    }
  });

  it('an empty window has no rows', () => {
    const { dictionary, array } = columns(Column.Schema.str);
    const typed = Column.ofIntArray(Int32Array.of(1, 2, 3));
    for (const c of [dictionary, array, typed] as Column<any>[]) {
      expect(Column.window(c, c.rowCount, c.rowCount).rowCount).toBe(0);
      expect(Column.window(c, 1, 1).rowCount).toBe(0);
    }
  });

  it('reads non-string values of a str dictionary as an array column does', () => {
    const values = [null, 1, 'a', undefined] as any[];
    for (const schema of [Column.Schema.str, Column.Schema.ustr]) {
      const dictionary = Column.ofDictionary({ indices: Uint8Array.of(0, 1, 2, 3), values, schema });
      const array = Column.ofArray({ array: values, schema });
      expect(read(dictionary)).toEqual(read(array));
    }
  });

  it('copies the values of a column, not storage it reads otherwise', () => {
    const strings = Column.ofStringArray(['10', '20']);
    const ints = {
      ...strings,
      schema: Column.Schema.int,
      value: (row: number) => +strings.value(row),
    } as Column<number>;
    const out = [0, 0];
    Column.copyToArray(ints, out);
    expect(out).toEqual([10, 20]);
  });
});

describe('column views and windows', () => {
  it('read toArray from the requested start', () => {
    const c = Column.ofLambda({ value: (i) => i * 10, rowCount: 10, schema: Column.Schema.int });
    expect(Array.from(Column.window(c, 2, 8).toArray({ start: 1, end: 3 }))).toEqual([30, 40]);
    expect(Array.from(Column.view(c, [9, 8, 7, 6]).toArray({ start: 2, end: 4 }))).toEqual([70, 60]);
  });

  it('a window of a typed column asks the column for an array of another kind', () => {
    const source = Column.ofFloatArray(Float64Array.of(0.1, 0.2, 0.3, 0.4));
    // as a BinaryCIF column: Float32 straight from the encoding, not a copy
    const float32 = Float32Array.of(1, 2, 3, 4);
    const column: Column<number> = {
      ...source,
      toArray: (params) =>
        params?.array === Float32Array ? float32.subarray(params.start ?? 0, params.end ?? 4) : source.toArray(params),
    };
    const window = Column.window(column, 1, 3);
    expect(Array.from(window.toArray({ array: Float32Array }))).toEqual([2, 3]);
    expect(Array.from(window.toArray({ array: Float32Array, start: 1, end: 2 }))).toEqual([3]);
    expect(window.toArray()).toBeInstanceOf(Float64Array);
  });
});
