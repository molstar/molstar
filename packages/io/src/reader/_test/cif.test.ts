/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it, test } from 'vitest';
import { CIF } from '../cif.js';
import * as Data from '../cif/data-model.js';
import * as Schema from '../cif/schema.js';
import { Column } from '@molstar/core/data/db';
import { parseCifText } from '../cif/text/parser.js';

const columnData = `123abc d,e,f '4 5 6'`;
// 123abc d,e,f '4 5 6'

const intField = Data.CifField.ofTokens({ data: columnData, indices: [0, 1, 1, 2, 2, 3], count: 3 });
const strField = Data.CifField.ofTokens({ data: columnData, indices: [3, 4, 4, 5, 5, 6], count: 3 });
const strListField = Data.CifField.ofTokens({ data: columnData, indices: [7, 12], count: 1 });
const intListField = Data.CifField.ofTokens({ data: columnData, indices: [14, 19], count: 1 });

const testBlock = Data.CifBlock(
  ['test'],
  {
    test: Data.CifCategory('test', 3, ['int', 'str', 'strList', 'intList'], {
      int: intField,
      str: strField,
      strList: strListField,
      intList: intListField,
    }),
  },
  'test',
);

namespace TestSchema {
  export const test = {
    int: Column.Schema.int,
    str: Column.Schema.str,
    strList: Column.Schema.List(',', (x) => x),
    intList: Column.Schema.List(' ', (x) => parseInt(x, 10)),
  };
  export const schema = { test };
}

test('cif triple quote', async () => {
  const data = `data_test
_test.field1 '''123 " '' 1'''
_test.field2 ''' c glide reflection through the plane (x,1/4,z)
chosen as one of the generators of the space group'''`;

  const result = await parseCifText(data).run();
  if (result.isError) {
    expect(false).toBe(true);
    return;
  }

  const cat = result.result.blocks[0].categories['test'];
  expect(cat.getField('field1')!.str(0)).toBe(`123 " '' 1`);
  expect(cat.getField('field2')!.str(0)).toBe(` c glide reflection through the plane (x,1/4,z)
chosen as one of the generators of the space group`);
});

test('cif case-insensitive field lookup', async () => {
  const data = `data_test
loop_
_atom_site.cartn_x
_atom_site.cartn_y
_atom_site.Cartn_z
1.0 2.0 3.0
4.0 5.0 6.0`;

  const result = await parseCifText(data).run();
  if (result.isError) {
    expect(false).toBe(true);
    return;
  }

  const atom_site = result.result.blocks[0].categories['atom_site'];
  expect(atom_site.getField('Cartn_x')!.float(0)).toBe(1.0);
  expect(atom_site.getField('Cartn_Y')!.float(1)).toBe(5.0);
  expect(atom_site.getField('cartn_z')!.float(1)).toBe(6.0);
});

describe('schema', () => {
  const db = Schema.toDatabase(TestSchema.schema, testBlock);
  it('property access', () => {
    const { int, str, strList, intList } = db.test;
    expect(int.value(0)).toBe(1);
    expect(str.value(1)).toBe('b');
    expect(strList.value(0)).toEqual(['d', 'e', 'f']);
    expect(intList.value(0)).toEqual([4, 5, 6]);
  });

  it('toArray', () => {
    const ret = db.test.int.toArray({ array: Int32Array });
    expect(ret.length).toBe(3);
    expect(ret[0]).toBe(1);
    expect(ret[1]).toBe(2);
    expect(ret[2]).toBe(3);
  });
});

describe('CIF.parse', () => {
  it('reads text CIF from a string or bytes and BinaryCIF from bytes', async () => {
    const text = `data_x\nloop_\n_a.v\n1\n2\n`;
    const binary = fs.readFileSync(path.resolve(import.meta.dirname, '../../../../../data/examples/1cbs_full.bcif'));
    for (const input of [text, new TextEncoder().encode(text), new TextEncoder().encode(`﻿${text}`)]) {
      const r = await CIF.parse(input).run();
      if (r.isError) throw new Error(r.message);
      expect(r.result.blocks[0].categories['a'].getField('v')!.toIntArray()).toEqual([1, 2]);
    }
    // no bytes fail as BinaryCIF, as before bytes were sniffed
    expect((await CIF.parse(new Uint8Array(0)).run()).isError).toBe(true);
    const r = await CIF.parse(new Uint8Array(binary)).run();
    if (r.isError) throw new Error(r.message);
    expect(r.result.blocks[0].categories['atom_site'].getField('Cartn_x')!.binaryEncoding).toBeDefined();
  });
});
