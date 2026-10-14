/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { SchemaItemFields, isFieldValue, isSchemaItemField, itemWithValues, schemaItemsFromFields } from '../fields.js';

describe('schemaItemsFromFields', () => {
  it('is one row without columns, with the values it has', () => {
    expect(schemaItemsFromFields({ label_asym_id: 'A', label_seq_id: 3 })).toEqual([
      { label_asym_id: 'A', label_seq_id: 3 },
    ]);
    expect(schemaItemsFromFields({})).toEqual([{}]);
    expect(schemaItemsFromFields({ label_asym_id: 'A', label_seq_id: null, auth_seq_id: Number.NaN })).toEqual([
      { label_asym_id: 'A' },
    ]);
  });
  it('is a row for each position of the columns, a field which is not a column is in every row', () => {
    expect(
      schemaItemsFromFields({ label_asym_id: 'A', label_seq_id: [1, null, 3], label_atom_id: ['CA', 'N', undefined] }),
    ).toEqual([
      { label_asym_id: 'A', label_seq_id: 1, label_atom_id: 'CA' },
      { label_asym_id: 'A', label_atom_id: 'N' },
      { label_asym_id: 'A', label_seq_id: 3 },
    ]);
  });
  it('does not make a row from empty columns', () => {
    expect(schemaItemsFromFields({ label_asym_id: 'A', label_seq_id: [] })).toEqual([]);
  });
  it('rejects columns of different lengths and fields which are not known', () => {
    expect(() => schemaItemsFromFields({ label_asym_id: ['A', 'B'], label_seq_id: [1] })).toThrow('different lengths');
    expect(() => schemaItemsFromFields({ label_asym: 'A' })).toThrow("'label_asym' is not a field");
    expect(() => schemaItemsFromFields({ toString: 'A' })).toThrow('is not a field');
  });
});

describe('the fields of an item', () => {
  it('are all known, a missing value is none', () => {
    expect(SchemaItemFields).toContain('beg_auth_ins_code');
    expect(new Set(SchemaItemFields).size).toBe(SchemaItemFields.length);
    expect(SchemaItemFields.every(isSchemaItemField)).toBe(true);
    expect(isSchemaItemField('group_id')).toBe(false);
    expect([0, '', false, 'x'].every(isFieldValue)).toBe(true);
    expect([undefined, null, Number.NaN].some(isFieldValue)).toBe(false);
    expect(itemWithValues({ label_asym_id: 'A', label_seq_id: Number.NaN, group_id: 1 } as any)).toEqual({
      label_asym_id: 'A',
    });
  });
});
