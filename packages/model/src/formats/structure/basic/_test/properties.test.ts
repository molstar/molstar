/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Table } from '@molstar/core/data/db';
import { getMissingResidues } from '../properties.js';
import { BasicSchema, createBasic } from '../schema.js';

describe('getMissingResidues', () => {
  it('derives missing residues from entity_poly_seq when pdbx_unobs_or_zero_occ_residues is absent', () => {
    const data = createBasic({
      struct_asym: Table.ofArrays(BasicSchema.struct_asym, {
        id: ['A'],
        entity_id: ['1'],
      }),
      entity_poly_seq: Table.ofArrays(BasicSchema.entity_poly_seq, {
        entity_id: ['1', '1', '1'],
        num: [1, 2, 3],
        mon_id: ['ALA', 'GLY', 'VAL'],
        hetero: ['n', 'n', 'n'],
      }),
      // seq_id 2 has no matching atom_site coordinates, so it should be reported missing
      atom_site: Table.ofArrays(BasicSchema.atom_site, {
        label_asym_id: ['A', 'A'],
        label_seq_id: [1, 3],
        pdbx_PDB_model_num: [1, 1],
      }),
    });

    const missingResidues = getMissingResidues(data);

    expect(missingResidues.size).toBe(1);
    expect(missingResidues.has(1, 'A', 1)).toBe(false);
    expect(missingResidues.has(1, 'A', 2)).toBe(true);
    expect(missingResidues.has(1, 'A', 3)).toBe(false);
    expect(missingResidues.get(1, 'A', 2)).toEqual({ polymer_flag: 'y', occupancy_flag: 1 });
  });

  it('reports no missing residues when all entity_poly_seq positions have coordinates', () => {
    const data = createBasic({
      struct_asym: Table.ofArrays(BasicSchema.struct_asym, {
        id: ['A'],
        entity_id: ['1'],
      }),
      entity_poly_seq: Table.ofArrays(BasicSchema.entity_poly_seq, {
        entity_id: ['1', '1', '1'],
        num: [1, 2, 3],
        mon_id: ['ALA', 'GLY', 'VAL'],
        hetero: ['n', 'n', 'n'],
      }),
      atom_site: Table.ofArrays(BasicSchema.atom_site, {
        label_asym_id: ['A', 'A', 'A'],
        label_seq_id: [1, 2, 3],
        pdbx_PDB_model_num: [1, 1, 1],
      }),
    });

    const missingResidues = getMissingResidues(data);

    expect(missingResidues.size).toBe(0);
  });

  it('prefers pdbx_unobs_or_zero_occ_residues when present, even if entity_poly_seq is also available', () => {
    const data = createBasic({
      struct_asym: Table.ofArrays(BasicSchema.struct_asym, {
        id: ['A'],
        entity_id: ['1'],
      }),
      entity_poly_seq: Table.ofArrays(BasicSchema.entity_poly_seq, {
        entity_id: ['1', '1', '1'],
        num: [1, 2, 3],
        mon_id: ['ALA', 'GLY', 'VAL'],
        hetero: ['n', 'n', 'n'],
      }),
      atom_site: Table.ofArrays(BasicSchema.atom_site, {
        label_asym_id: ['A', 'A'],
        label_seq_id: [1, 3],
        pdbx_PDB_model_num: [1, 1],
      }),
      pdbx_unobs_or_zero_occ_residues: Table.ofArrays(BasicSchema.pdbx_unobs_or_zero_occ_residues, {
        id: [1],
        polymer_flag: ['y'],
        occupancy_flag: [1],
        PDB_model_num: [1],
        label_asym_id: ['A'],
        label_seq_id: [2],
      }),
    });

    const missingResidues = getMissingResidues(data);

    expect(missingResidues.size).toBe(1);
    expect(missingResidues.has(1, 'A', 2)).toBe(true);
  });

  it('derives them from atom rows of several atoms per residue, chains, models and rows without a seq id', () => {
    const data = createBasic({
      struct_asym: Table.ofArrays(BasicSchema.struct_asym, {
        id: ['A', 'B'],
        entity_id: ['1', '1'],
      }),
      entity_poly_seq: Table.ofArrays(BasicSchema.entity_poly_seq, {
        entity_id: ['1', '1', '1', '1'],
        num: [1, 2, 3, 4],
        mon_id: ['ALA', 'GLY', 'VAL', 'SER'],
        hetero: ['n', 'n', 'n', 'n'],
      }),
      // model 1: A 1, 1, 3, (water), A 3 again after B, B 2, 2; model 2: A 4, B 1
      atom_site: Table.ofArrays(BasicSchema.atom_site, {
        label_asym_id: ['A', 'A', 'A', 'W', 'B', 'B', 'A', 'A', 'B'],
        label_seq_id: [1, 1, 3, 0, 2, 2, 3, 4, 1],
        pdbx_PDB_model_num: [1, 1, 1, 1, 1, 1, 1, 2, 2],
      }),
    });
    // the water row has no seq id
    (data.atom_site.label_seq_id as any).valueKind = (i: number) => (i === 3 ? 1 : 0);

    const missingResidues = getMissingResidues(data);

    const missing: string[] = [];
    for (const m of [1, 2]) {
      for (const a of ['A', 'B']) {
        for (const seq of [1, 2, 3, 4]) if (missingResidues.has(m, a, seq)) missing.push(`${m}${a}${seq}`);
      }
    }
    expect(missing).toEqual(['1A2', '1A4', '1B1', '1B3', '1B4', '2A1', '2A2', '2A3', '2B2', '2B3', '2B4']);
    expect(missingResidues.size).toBe(11);
  });
});
