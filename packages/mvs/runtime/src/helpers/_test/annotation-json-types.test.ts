/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { MVSAnnotation } from '../../components/annotation-prop.js';
import { createModel } from './annotation-rows-fixture.js';

/** JSON annotations are not type-checked, numbers can arrive as strings and identifiers as numbers. */
describe('JSON annotation values of the wrong type', async () => {
  const model = await createModel({ nChains: 2, nResidues: 20 });
  const nAtoms = model.atomicHierarchy.atoms._rowCount;
  const rowOfEachAtom = (data: any) => {
    const annot = new MVSAnnotation({ format: 'json', data }, 'all_atomic', {});
    const atoms = (annot as any).getRowForEachAtom(model, '1').atoms as Int32Array | null;
    return Array.from(atoms ?? new Array(nAtoms).fill(-1));
  };
  const selected = (data: any) => rowOfEachAtom(data).flatMap((v, i) => (v >= 0 ? [i] : []));
  const residue10OfA = Array.from({ length: 6 }, (_, i) => 9 * 6 + i); // 6 atoms per residue, residues from 1

  it('a numeric string selects what the number does', () => {
    expect(selected([{ label_asym_id: 'A', label_seq_id: 10 }])).toEqual(residue10OfA);
    expect(selected([{ label_asym_id: 'A', label_seq_id: '10' }])).toEqual(residue10OfA);
    expect(selected([{ label_asym_id: 'A', label_seq_id: ' 10 ' }])).toEqual(residue10OfA);
  });
  it('range bounds given as strings are not concatenated with 1', () => {
    const numeric = selected([{ label_asym_id: 'A', beg_label_seq_id: 5, end_label_seq_id: 7 }]);
    expect(numeric.length).toBe(3 * 6);
    expect(selected([{ label_asym_id: 'A', beg_label_seq_id: '5', end_label_seq_id: '7' }])).toEqual(numeric);
    expect(selected([{ label_asym_id: 'A', end_label_seq_id: '1' }])).toEqual(
      selected([{ label_asym_id: 'A', end_label_seq_id: 1 }]),
    );
    expect(selected([{ label_asym_id: 'A', end_label_seq_id: '1' }]).length).toBe(6);
  });
  it('the object-of-arrays form is converted the same way', () => {
    expect(selected({ label_asym_id: ['A', 'A'], label_seq_id: ['10', 11] }).length).toBe(2 * 6);
  });
  it('a number given for an identifier selects the identifier with that text', () => {
    // entity ids are '1'..'5', chain A is in entity 1
    expect(selected([{ label_entity_id: 1, label_asym_id: 'A' }])).toEqual(
      selected([{ label_entity_id: '1', label_asym_id: 'A' }]),
    );
    expect(selected([{ label_entity_id: 1 }]).length).toBeGreaterThan(0);
  });
  it('values which are not numbers select nothing and non-integers match no residue', () => {
    expect(selected([{ label_asym_id: 'A', label_seq_id: 'ten' }])).toEqual([]);
    expect(selected([{ label_asym_id: 'A', label_seq_id: 10.5 }])).toEqual([]);
    expect(selected([{ label_asym_id: 'A', beg_label_seq_id: 'x' }])).toEqual([]);
  });
});
