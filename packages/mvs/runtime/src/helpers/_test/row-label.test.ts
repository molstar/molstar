/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { describeRow, describeRows } from '../row-label.js';

describe('describeRow', () => {
  it('names chains, residues and ranges in a compact form', () => {
    expect(describeRow({ label_asym_id: 'A' })).toBe('A');
    expect(describeRow({ label_entity_id: '2' })).toBe('Entity 2');
    expect(describeRow({ label_asym_id: 'A', label_seq_id: 432 })).toBe('A432');
    expect(describeRow({ label_asym_id: 'A', beg_label_seq_id: 10, end_label_seq_id: 50 })).toBe('A:10-50');
    expect(describeRow({ label_asym_id: 'A', beg_label_seq_id: 10, end_label_seq_id: 10 })).toBe('A:10');
    expect(describeRow({ label_asym_id: 'A', beg_label_seq_id: 10 })).toBe('A:>=10');
    expect(describeRow({ label_asym_id: 'A', end_label_seq_id: 50 })).toBe('A:<=50');
  });
  it('keeps negative numbers and one-sided ranges apart', () => {
    expect(describeRow({ end_label_seq_id: 20 })).toBe('res <=20');
    expect(describeRow({ label_seq_id: -20 })).toBe('res -20');
    expect(describeRow({ label_asym_id: 'A', label_seq_id: -5 })).toBe('A:-5');
    expect(describeRow({ label_asym_id: 'A', beg_label_seq_id: -5, end_label_seq_id: 3 })).toBe('A:-5..3');
    expect(describeRow({ label_asym_id: 'A', beg_label_seq_id: 3, end_label_seq_id: -5 })).toBe('A:3..-5');
    expect(describeRow({ auth_asym_id: 'B', beg_auth_seq_id: -5, end_auth_seq_id: 3 })).toBe('B:-5..3 (auth)');
    expect(describeRow({ beg_label_seq_id: 10 })).toBe('res >=10');
  });
  it('empty texts are shown, and an empty insertion code is told with the residue', () => {
    expect(describeRow({ label_asym_id: '' })).toBe('""');
    expect(describeRow({ label_asym_id: 'A', label_comp_id: '' })).toBe('A ""');
    expect(describeRow({ auth_asym_id: 'A', auth_seq_id: 10, pdbx_PDB_ins_code: '' })).toBe('A10 no ins code (auth)');
    expect(describeRow({ pdbx_PDB_ins_code: '' })).toBe('no ins code');
    expect(describeRow({ auth_asym_id: 'A', auth_seq_id: 10, pdbx_PDB_ins_code: 'B' })).toBe('A10B (auth)');
  });
  it('marks criteria which use the auth numbering once, or per piece when label and auth are mixed', () => {
    expect(describeRow({ auth_asym_id: 'H' })).toBe('H (auth)');
    expect(describeRow({ auth_asym_id: 'B', auth_seq_id: 5, pdbx_PDB_ins_code: 'A' })).toBe('B5A (auth)');
    expect(describeRow({ auth_asym_id: 'B', beg_auth_seq_id: 3, end_auth_seq_id: 9 })).toBe('B:3-9 (auth)');
    expect(describeRow({ auth_asym_id: 'A', auth_atom_id: 'CA' })).toBe('A CA (auth)');
    expect(describeRow({ label_asym_id: 'A', auth_seq_id: 5 })).toBe('A 5 (auth)');
    expect(describeRow({ label_comp_id: 'ALA', auth_comp_id: 'XALA' })).toBe('ALA XALA (auth)');
    expect(describeRow({ label_comp_id: 'ALA', auth_comp_id: 'ALA' })).toBe('ALA');
  });
  it('names residues without a chain, atoms, components and instances', () => {
    expect(describeRow({ beg_label_seq_id: 3, end_label_seq_id: 9 })).toBe('res 3-9');
    expect(describeRow({ label_seq_id: 3 })).toBe('res 3');
    expect(describeRow({ beg_auth_seq_id: 3, end_auth_seq_id: 9 })).toBe('res 3-9 (auth)');
    expect(describeRow({ pdbx_PDB_ins_code: 'A' })).toBe('ins A');
    expect(describeRow({ residue_index: 7 })).toBe('index 7');
    expect(describeRow({ label_comp_id: 'HEM' })).toBe('HEM');
    expect(describeRow({ label_asym_id: 'A', label_comp_id: 'HEM' })).toBe('A HEM');
    expect(describeRow({ label_asym_id: 'A', label_atom_id: 'CA' })).toBe('A CA');
    expect(describeRow({ label_asym_id: 'A', label_seq_id: 432, label_atom_id: 'CA' })).toBe('A432 CA');
    expect(describeRow({ type_symbol: 'c' })).toBe('C atoms');
    expect(describeRow({ atom_id: 12 })).toBe('atom 12');
    expect(describeRow({ atom_index: 3 })).toBe('atom index 3');
    expect(describeRow({ label_asym_id: 'A', instance_id: '1_555' })).toBe('A instance 1_555');
    expect(describeRow({ label_asym_id: 'A', operator_name: '2_555' })).toBe('A operator 2_555');
  });
  it('a row without criteria selects everything, unrelated fields are ignored', () => {
    expect(describeRow({})).toBe('All');
    expect(describeRow({ group_id: 'g' })).toBe('All');
    expect(describeRow({ label_asym_id: undefined, label_seq_id: undefined })).toBe('All');
  });
});

describe('describeRows', () => {
  it('describes a single row like describeRow, and none as empty', () => {
    expect(describeRows([{ label_asym_id: 'A' }])).toBe('A');
    expect(describeRows([])).toBe('Empty Selection');
  });
  it('lists chains, a few by name and many by count', () => {
    expect(describeRows([{ label_asym_id: 'A' }, { label_asym_id: 'B' }, { label_asym_id: 'C' }])).toBe('A, B, C');
    expect(describeRows([{ auth_asym_id: 'A' }, { auth_asym_id: 'B' }])).toBe('A, B (auth)');
    const many = Array.from({ length: 30 }, (_, i) => ({ label_asym_id: `C${i}` }));
    expect(describeRows(many)).toBe('30 chains');
  });
  it('lists up to three short rows, otherwise counts them', () => {
    expect(describeRows([{ label_asym_id: 'A', label_seq_id: 1 }, { label_comp_id: 'HEM' }])).toBe('A1 + HEM');
    expect(
      describeRows([
        { label_asym_id: 'A', label_seq_id: 432 },
        { label_asym_id: 'A', label_seq_id: 435 },
        { label_asym_id: 'B', beg_label_seq_id: 1, end_label_seq_id: 9 },
      ]),
    ).toBe('A432 + A435 + B:1-9');
    const long = ['A', 'B', 'C'].map((c) => ({
      label_asym_id: c,
      beg_label_seq_id: 100,
      end_label_seq_id: 200,
      label_atom_id: 'CA',
    }));
    expect(describeRows(long)).toBe('3 selections');
  });
  it('counts many residues, of one chain or of several', () => {
    const four = Array.from({ length: 4 }, (_, i) => ({ label_asym_id: 'A', label_seq_id: i + 1 }));
    expect(describeRows(four)).toBe('A: 4 residues');
    const authFour = Array.from({ length: 4 }, (_, i) => ({ auth_asym_id: 'A', auth_seq_id: i + 1 }));
    expect(describeRows(authFour)).toBe('A: 4 residues (auth)');
    const mixed = Array.from({ length: 120 }, (_, i) => ({ label_asym_id: i % 2 ? 'A' : 'B', label_seq_id: i + 1 }));
    expect(describeRows(mixed)).toBe('120 residues');
    const ranges = Array.from({ length: 4 }, (_, i) => ({
      label_asym_id: 'A',
      beg_label_seq_id: i,
      end_label_seq_id: i + 5,
    }));
    expect(describeRows(ranges)).toBe('A: 4 selections');
    expect(describeRows([...ranges, { label_comp_id: 'HEM' }])).toBe('5 selections');
  });
  it('never produces a label which grows with the number of rows', () => {
    for (const n of [1, 2, 3, 10, 1000, 100000]) {
      const rows = Array.from({ length: n }, (_, i) => ({
        label_asym_id: `Chain${i}`,
        label_seq_id: i,
        label_atom_id: 'CA',
      }));
      expect(describeRows(rows).length).toBeLessThan(40);
    }
  });
});
