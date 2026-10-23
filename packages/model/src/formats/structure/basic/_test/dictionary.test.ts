/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import type { Column, Table } from '@molstar/core/data/db';
import type { Model } from '../../../../model/model/model.js';
import { models, parseBinary, parseText, toText } from './encode.js';

// A BinaryCIF model, whose string columns are dictionaries (binary/field.ts, perf-audit-2 04/P16) and whose sorted
// `atom_site` gathers their indices and buckets by them (sort.ts), against the same file as text, which has no
// dictionaries: the hierarchy tables equal value by value.

const examples = path.resolve(import.meta.dirname, '../../../../../../../data/examples');

function expectSameTable(a: Table<any>, b: Table<any>, where: string) {
  expect(a._rowCount, where).toBe(b._rowCount);
  for (const k of a._columns) {
    const ca = a[k] as Column<any>;
    const cb = b[k] as Column<any>;
    expect(ca.isDefined, `${where}.${k}`).toBe(cb.isDefined);
    for (let i = 0; i < ca.rowCount; i++) {
      if (ca.value(i) !== cb.value(i) || ca.valueKind(i) !== cb.valueKind(i)) {
        expect([i, ca.value(i), ca.valueKind(i)], `${where}.${k}`).toEqual([i, cb.value(i), cb.valueKind(i)]);
      }
    }
  }
}

function expectSameHierarchy(a: Model, b: Model, where: string) {
  const ha = a.atomicHierarchy;
  const hb = b.atomicHierarchy;
  expectSameTable(ha.atoms, hb.atoms, `${where} atoms`);
  expectSameTable(ha.residues, hb.residues, `${where} residues`);
  expectSameTable(ha.chains, hb.chains, `${where} chains`);
  expect(Array.from(ha.atomSourceIndex.toArray()), where).toEqual(Array.from(hb.atomSourceIndex.toArray()));
  expect(Array.from(ha.residueAtomSegments.offsets), where).toEqual(Array.from(hb.residueAtomSegments.offsets));
  expect(Array.from(ha.chainAtomSegments.offsets), where).toEqual(Array.from(hb.chainAtomSegments.offsets));
}

describe('BinaryCIF string dictionaries', () => {
  it('give the hierarchy of the same file as text', async () => {
    // 21nr sorts; 1tup and 1cbs are in order
    for (const name of ['1tup.bcif', '21nr.bcif', '1cbs_full.bcif']) {
      const binary = await parseBinary(new Uint8Array(fs.readFileSync(path.join(examples, name))));
      const field = binary.blocks[0].categories['atom_site'].getField('label_atom_id')!;
      expect(field.storage!()!.kind).toBe('dictionary');
      const text = await parseText(toText(binary));
      const [mb] = await models(binary);
      const [mt] = await models(text);
      expectSameHierarchy(mb, mt, name);
    }
  });
});
