/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { describeRow, describeRows } from '../row-label.js';

describe('describeRow', () => {
  it('names the usual selections the way a person would', () => {
    expect(describeRow({ label_asym_id: 'A' })).toBe('Chain A');
    expect(describeRow({ label_entity_id: '2' })).toBe('Entity 2');
    expect(describeRow({ auth_asym_id: 'H' })).toBe('Chain H (auth)');
    expect(describeRow({ label_asym_id: 'A', label_seq_id: 42 })).toBe('Chain A, residue 42');
    expect(describeRow({ label_asym_id: 'A', beg_label_seq_id: 10, end_label_seq_id: 50 })).toBe(
      'Chain A, residues 10-50',
    );
    expect(describeRow({ label_asym_id: 'A', beg_label_seq_id: 10, end_label_seq_id: 10 })).toBe(
      'Chain A, residues 10',
    );
    expect(describeRow({ label_asym_id: 'A', beg_label_seq_id: 10 })).toBe('Chain A, residues from 10');
    expect(describeRow({ label_asym_id: 'A', end_label_seq_id: 50 })).toBe('Chain A, residues up to 50');
  });
  it('names residue, atom, and instance criteria', () => {
    expect(describeRow({ auth_asym_id: 'B', auth_seq_id: 5, pdbx_PDB_ins_code: 'A' })).toBe(
      'Chain B (auth), auth residue 5A',
    );
    expect(describeRow({ pdbx_PDB_ins_code: 'A' })).toBe('Insertion code A');
    expect(describeRow({ beg_auth_seq_id: 3, end_auth_seq_id: 9 })).toBe('Auth residues 3-9');
    expect(describeRow({ residue_index: 7 })).toBe('Residue index 7');
    expect(describeRow({ label_comp_id: 'HEM' })).toBe('HEM');
    expect(describeRow({ label_comp_id: 'ALA', auth_comp_id: 'ALA' })).toBe('ALA');
    expect(describeRow({ label_comp_id: 'ALA', auth_comp_id: 'XALA' })).toBe('ALA, XALA (auth)');
    expect(describeRow({ label_asym_id: 'A', label_atom_id: 'CA' })).toBe('Chain A, atom CA');
    expect(describeRow({ type_symbol: 'c' })).toBe('C atoms');
    expect(describeRow({ atom_id: 12 })).toBe('Atom id 12');
    expect(describeRow({ atom_index: 3 })).toBe('Atom index 3');
    expect(describeRow({ label_asym_id: 'A', instance_id: '1_555' })).toBe('Chain A, instance 1_555');
  });
  it('a row without criteria selects everything, unrelated fields are ignored', () => {
    expect(describeRow({})).toBe('All');
    expect(describeRow({ group_id: 'g' })).toBe('All');
    expect(describeRow({ label_asym_id: undefined, label_seq_id: undefined })).toBe('All');
  });
});

describe('describeRows', () => {
  it('describes a single row like describeRow, and none as empty', () => {
    expect(describeRows([{ label_asym_id: 'A' }])).toBe('Chain A');
    expect(describeRows([])).toBe('Empty Selection');
  });
  it('lists chains, a few chains by name and many by count', () => {
    expect(describeRows([{ label_asym_id: 'A' }, { label_asym_id: 'B' }, { label_asym_id: 'C' }])).toBe(
      'Chains A, B, C',
    );
    const many = Array.from({ length: 30 }, (_, i) => ({ label_asym_id: `C${i}` }));
    expect(describeRows(many)).toBe('30 chains');
  });
  it('lists up to three short rows, otherwise counts them', () => {
    expect(describeRows([{ label_asym_id: 'A', label_seq_id: 1 }, { label_comp_id: 'HEM' }])).toBe(
      'Chain A, residue 1 + HEM',
    );
    const long = [
      { label_asym_id: 'A', beg_label_seq_id: 100, end_label_seq_id: 200, label_atom_id: 'CA' },
      { label_asym_id: 'B', beg_label_seq_id: 100, end_label_seq_id: 200, label_atom_id: 'CA' },
      { label_asym_id: 'C', beg_label_seq_id: 100, end_label_seq_id: 200, label_atom_id: 'CA' },
    ];
    expect(describeRows(long)).toBe('3 selections');
    const four = Array.from({ length: 4 }, (_, i) => ({ label_asym_id: 'A', label_seq_id: i + 1 }));
    expect(describeRows(four)).toBe('Chain A: 4 selections');
    const mixed = Array.from({ length: 120 }, (_, i) => ({ label_asym_id: i % 2 ? 'A' : 'B', label_seq_id: i + 1 }));
    expect(describeRows(mixed)).toBe('120 selections');
  });
  it('never produces a label which grows with the number of rows', () => {
    for (const n of [1, 2, 3, 10, 1000, 100000]) {
      const rows = Array.from({ length: n }, (_, i) => ({
        label_asym_id: `Chain${i}`,
        label_seq_id: i,
        label_atom_id: 'CA',
      }));
      expect(describeRows(rows).length).toBeLessThan(80);
    }
  });
});
