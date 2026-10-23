/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Aniruddha Adak <aniruddhaadak80@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Column } from '@molstar/core/data/db';
import {
  ArrayEncoder,
  ArrayEncoding as E,
  type EncodedColumn,
  classifyIntArray,
  classifyFloatArray,
} from '../../../../common/binary-cif.js';
import { Field } from '../field.js';
import { CifCategory } from '../../data-model.js';
import { toTable } from '../../schema.js';
import { decodeSteps } from '../../../../common/binary-cif/decoder.js';
import { parseFloat as fastParseFloat, parseInt as fastParseInt } from '@molstar/core/util/number-parser';

/** The array behind a column or field: its `'array'` storage. */
function arrayOf(c: { storage?(): { kind: string; array?: ArrayLike<any> } | undefined }) {
  const s = c.storage?.();
  return s?.kind === 'array' ? s.array : void 0;
}

/** mirrors `Column.ValueKinds` */
const Present = 0,
  NotPresent = 1,
  Unknown = 2;

function intColumn(data: number[], mask?: number[]): EncodedColumn {
  const array = new Int32Array(data);
  return {
    name: 'test',
    data: classifyIntArray(array).encode(array),
    mask: mask ? ArrayEncoder.by(E.byteArray).encode(new Uint8Array(mask)) : void 0,
  };
}

function floatColumn(data: number[], mask?: number[]): EncodedColumn {
  const array = new Float32Array(data);
  return {
    name: 'test',
    data: classifyFloatArray(array).encode(array),
    mask: mask ? ArrayEncoder.by(E.byteArray).encode(new Uint8Array(mask)) : void 0,
  };
}

const IntValues = [1, -2, 3, -4, 5];
const FloatValues = [1.5, -2.5, 3.5, -4.5, 5.5];
/** rows 1 (`?`) and 3 (`.`) are missing, the rest is present */
const Mask = [Present, Unknown, Present, NotPresent, Present];

describe('binary-cif Field', () => {
  it('masked int values default to 0', () => {
    const f = Field(intColumn(IntValues, Mask));

    expect(f.rowCount).toBe(5);
    expect(f.int(0)).toBe(1);
    expect(f.int(1)).toBe(0);
    expect(f.int(2)).toBe(3);
    expect(f.int(3)).toBe(0);
    expect(f.int(4)).toBe(5);
  });

  it('masked float values default to 0', () => {
    const f = Field(floatColumn(FloatValues, Mask));

    expect(f.float(0)).toBe(1.5);
    expect(f.float(1)).toBe(0);
    expect(f.float(2)).toBe(3.5);
    expect(f.float(3)).toBe(0);
    expect(f.float(4)).toBe(5.5);
  });

  it('toIntArray defaults masked values to 0 and keeps the source data intact', () => {
    const f = Field(intColumn(IntValues, Mask));

    expect(Array.from(f.toIntArray())).toEqual([1, 0, 3, 0, 5]);
    expect(Array.from(f.toIntArray({ start: 1, end: 4 }))).toEqual([0, 3, 0]);

    expect(Array.from(arrayOf(f)!)).toEqual(IntValues);
  });

  it('toFloatArray defaults masked values to 0', () => {
    const f = Field(floatColumn(FloatValues, Mask));

    expect(Array.from(f.toFloatArray())).toEqual([1.5, 0, 3.5, 0, 5.5]);
    expect(Array.from(f.toFloatArray({ start: 1, end: 4 }))).toEqual([0, 3.5, 0]);

    expect(Array.from(arrayOf(f)!)).toEqual(FloatValues);
  });

  it('leaves str and valueKind unchanged', () => {
    const f = Field(intColumn([1, 2, 3], [Present, Unknown, NotPresent]));

    expect(f.str(0)).toBe('1');
    expect(f.str(1)).toBe('');
    expect(f.str(2)).toBe('');

    expect(f.valueKind(0)).toBe(Present);
    expect(f.valueKind(1)).toBe(Unknown);
    expect(f.valueKind(2)).toBe(NotPresent);
  });

  it('returns the encoded values for a trivial mask', () => {
    const f = Field(intColumn([1, 2, 3], [Present, Present, Present]));

    expect(f.int(1)).toBe(2);
    expect(Array.from(f.toIntArray())).toEqual([1, 2, 3]);
    expect(f.valueKind(0)).toBe(Present);
  });

  it('returns the encoded values when there is no mask', () => {
    const f = Field(intColumn([1, 2, 3]));

    expect(f.int(1)).toBe(2);
    expect(Array.from(f.toIntArray())).toEqual([1, 2, 3]);
    expect(f.valueKind(0)).toBe(Present);
  });

  it('decodes on the first access when given the row count, with the same values', () => {
    const column = intColumn(IntValues, Mask);
    // an encoding the decoder rejects: creating the field must not decode it
    const broken: EncodedColumn = {
      ...column,
      data: { ...column.data, encoding: [{ kind: 'ByteArray', type: 99 as any }] },
    };
    const lazy = Field(broken, 5);
    expect(lazy.rowCount).toBe(5);
    expect(lazy.binaryEncoding).toBe(broken.data.encoding);
    expect(() => lazy.int(0)).toThrow();

    const eager = Field(column);
    const f = Field(column, 5);
    for (let i = 0; i < 5; i++) {
      expect([f.str(i), f.int(i), f.float(i), f.valueKind(i)]).toEqual([
        eager.str(i),
        eager.int(i),
        eager.float(i),
        eager.valueKind(i),
      ]);
    }
    expect(f.areValuesEqual(0, 2)).toBe(eager.areValuesEqual(0, 2));
    expect(Array.from(f.toIntArray())).toEqual(Array.from(eager.toIntArray()));
    expect(Array.from(f.toFloatArray({ start: 1, end: 3 }))).toEqual(
      Array.from(eager.toFloatArray({ start: 1, end: 3 })),
    );
    expect(f.toStringArray()).toEqual(eager.toStringArray());
    expect(Array.from(arrayOf(f)!)).toEqual(Array.from(arrayOf(eager)!));
    // decoded once
    expect(arrayOf(f)).toBe(arrayOf(f));
  });

  it('decodes a quantized column read as Float32 straight to Float32, the other accessors to the source type', () => {
    const values = Float64Array.from({ length: 100 }, (_, i) => Math.round(Math.sin(i) * 1e6) / 1000);
    const column: EncodedColumn = {
      name: 'x',
      data: E.by(E.fixedPoint(1000)).and(E.delta).and(E.integerPacking).encode(values),
      mask: ArrayEncoder.by(E.byteArray).encode(Uint8Array.from(values, (_, i) => (i === 3 ? NotPresent : Present))),
    };
    const f = Field(column, values.length);
    const float32 = f.toFloatArray({ array: Float32Array });
    expect(float32).toBeInstanceOf(Float32Array);
    const float64 = f.toFloatArray();
    expect(float64).toBeInstanceOf(Float64Array);
    expect(Array.from(float32)).toEqual(Array.from(float64, (v) => Math.fround(v)));
    expect(float32[3]).toBe(0);
    expect(arrayOf(f)).toBeInstanceOf(Float64Array);
    expect(f.float(5)).toBe(values[5]);
    expect(Array.from(f.toFloatArray({ array: Float32Array, start: 10, end: 20 }))).toEqual(
      Array.from(float32.slice(10, 20)),
    );
  });

  it('keeps a string column as its dictionary, with the values of its array of strings', () => {
    const data = ['ALA', 'GLY', null, '12', 'GLY', '', 'ALA', '3.5'] as string[];
    const encoded = E.by(E.stringArray).encode(data);
    const steps = decodeSteps(encoded) as string[];
    const cases: [string, number[] | undefined, 'dictionary' | 'array'][] = [
      ['no mask', void 0, 'dictionary'],
      // the null row and the empty string masked: both read as ''
      [
        'mask over empty rows',
        [Present, Present, NotPresent, Present, Present, Unknown, Present, Present],
        'dictionary',
      ],
      // a masked row with a string under it keeps the array of strings
      ['mask over a value', [Present, Unknown, NotPresent, Present, Present, Present, Present, Present], 'array'],
    ];
    for (const [name, mask, kind] of cases) {
      const column: EncodedColumn = {
        name: 'x',
        data: encoded,
        mask: mask ? ArrayEncoder.by(E.byteArray).encode(new Uint8Array(mask)) : void 0,
      };
      const f = Field(column, data.length);
      expect(f.storage!()!.kind, name).toBe(kind);
      for (let i = 0; i < data.length; i++) {
        const present = !mask || mask[i] === Present;
        const v = steps[i];
        expect([f.str(i), f.int(i), f.float(i), f.valueKind(i)], `${name} ${i}`).toEqual([
          present ? v : '',
          fastParseInt(v, 0, v.length),
          fastParseFloat(v, 0, v.length),
          mask ? mask[i] : Present,
        ]);
        for (let j = 0; j < data.length; j++) expect(f.areValuesEqual(i, j)).toBe(steps[i] === steps[j]);
      }
    }
  });

  it('reads a dictionary through CIF schemas: transforms, ints over strings, views and windows', () => {
    const names = ['Ala', 'gly', 'ALA', '7', '12', 'Ala'];
    const f = Field({ name: 'x', data: E.by(E.stringArray).encode(names) }, names.length);
    expect(f.storage!()!.kind).toBe('dictionary');
    const category = CifCategory('t', names.length, ['l', 'u', 'n', 's'], { l: f, u: f, n: f, s: f });
    const table = toTable(
      { l: Column.Schema.lstr, u: Column.Schema.ustr, n: Column.Schema.int, s: Column.Schema.str },
      category,
    ) as any;
    const map = [5, 3, 0, 4];
    const expected = {
      l: names.map((v) => v.toLowerCase()),
      u: names.map((v) => v.toUpperCase()),
      n: names.map((v) => fastParseInt(v, 0, v.length)),
      s: names,
    } as Record<string, (string | number)[]>;
    for (const k of ['l', 'u', 'n', 's']) {
      const c = table[k] as Column<any>;
      expect(
        Array.from({ length: c.rowCount }, (_, i) => c.value(i)),
        k,
      ).toEqual(expected[k]);
      const view = Column.view(c, map);
      expect(
        Array.from({ length: view.rowCount }, (_, i) => view.value(i)),
        k,
      ).toEqual(map.map((i) => expected[k][i]));
      const window = Column.window(c, 2, 5);
      expect(
        Array.from({ length: window.rowCount }, (_, i) => window.value(i)),
        k,
      ).toEqual(expected[k].slice(2, 5));
    }
  });
});
