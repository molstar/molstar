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

describe('column views and windows', () => {
  it('read toArray from the requested start', () => {
    const c = Column.ofLambda({ value: (i) => i * 10, rowCount: 10, schema: Column.Schema.int });
    expect(Array.from(Column.window(c, 2, 8).toArray({ start: 1, end: 3 }))).toEqual([30, 40]);
    expect(Array.from(Column.view(c, [9, 8, 7, 6]).toArray({ start: 2, end: 4 }))).toEqual([70, 60]);
  });
});
