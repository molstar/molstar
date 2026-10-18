/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Structure } from '@molstar/model/model/structure';
import { MVSAnnotation } from '../../components/annotation-prop.js';
import { IndicesAndSortings } from '@molstar/model/model/structure/structure/element/schema-selection/indexing';
import { substructureFromRows } from '@molstar/model/model/structure/structure/element/schema-selection/row-structure';
import type { MVSAnnotationRow } from '../schemas.js';
import { getAtomRangesForRow } from '@molstar/model/model/structure/structure/element/schema-selection/selections';
import { TestOperator, createModel } from './annotation-rows-fixture.js';

/**
 * Explicit expectations (not derived from the implementation or from MolScript) for rows which refer to values which some residues do not have.
 * The structure has 3 polymer chains of 20 residues (6 atoms each; every 7th residue has the insertion code A) and 2 ligands of 4 atoms,
 * which have no `label_seq_id` (and no insertion code).
 */
describe('values which are missing', async () => {
  const model = await createModel({ nChains: 3, nResidues: 20, insertionCodes: true, ligands: 2 });
  const structure = Structure.ofModel(model);
  const indices = IndicesAndSortings.get(model);
  const nLigandAtoms = 2 * 4;
  const nPolymerAtoms = 3 * 20 * 6;

  /** The number of selected atoms by the three ways a row can be resolved, which must agree */
  function counts(row: MVSAnnotationRow) {
    const ranges = getAtomRangesForRow(row, model, TestOperator, indices);
    let viaRanges = 0;
    if (ranges) for (let i = 0; i < ranges.from.length; i++) viaRanges += ranges.to[i] - ranges.from[i];
    const viaSubstructure = substructureFromRows(structure, [row]).elementCount;
    const annot = new MVSAnnotation({ format: 'json', data: [row] as any }, 'all_atomic', {});
    const mapping = (annot as any).getRowForEachAtom(model, TestOperator).atoms as Int32Array | null;
    const viaAnnotation = mapping ? mapping.reduce((n, v) => n + (v >= 0 ? 1 : 0), 0) : 0;
    expect([viaSubstructure, viaAnnotation]).toEqual([viaRanges, viaRanges]);
    return viaRanges;
  }
  it('a residue without label_seq_id never matches a label_seq_id criterion, not even 0', () => {
    expect(counts({ label_seq_id: 0 })).toBe(0);
    expect(counts({ label_asym_id: 'L0', label_seq_id: 1 })).toBe(0);
    expect(counts({ beg_label_seq_id: -5, end_label_seq_id: 5 })).toBe(3 * 5 * 6);
    expect(counts({ end_label_seq_id: 3 })).toBe(3 * 3 * 6);
    expect(counts({ beg_label_seq_id: 1 })).toBe(nPolymerAtoms);
    expect(counts({ beg_label_seq_id: 0, end_label_seq_id: 0 })).toBe(0);
  });
  it('... also when another criterion finds the residues and the label_seq_id is only checked on them', () => {
    const ligandResidueIndex = 3 * 20; // (the residues of the polymer chains come first, the first ligand is the next one)
    expect(counts({ residue_index: ligandResidueIndex, label_seq_id: 0 })).toBe(0);
    expect(counts({ residue_index: ligandResidueIndex, end_label_seq_id: 5 })).toBe(0);
    expect(counts({ residue_index: ligandResidueIndex, beg_label_seq_id: -5, end_label_seq_id: 5 })).toBe(0);
    expect(counts({ auth_seq_id: 500, end_label_seq_id: 5 })).toBe(0);
    expect(counts({ auth_seq_id: 500, beg_label_seq_id: -1 })).toBe(0);
    expect(counts({ pdbx_PDB_ins_code: '', label_comp_id: 'LIG', end_label_seq_id: 5 })).toBe(0);
    // (the residue is found when the criteria are about what it has)
    expect(counts({ residue_index: ligandResidueIndex, auth_seq_id: 500 })).toBe(4);
    expect(counts({ auth_seq_id: 500, auth_asym_id: 'L0' })).toBe(4);
  });
  it('ligands are found by what they have', () => {
    expect(counts({ label_asym_id: 'L0' })).toBe(4);
    expect(counts({ label_comp_id: 'LIG' })).toBe(nLigandAtoms);
    expect(counts({ auth_asym_id: 'L1', auth_seq_id: 501 })).toBe(4);
    expect(counts({ auth_seq_id: 500 })).toBe(4);
    expect(counts({ label_entity_id: '6' })).toBe(nLigandAtoms);
    expect(counts({ label_asym_id: 'L1', label_atom_id: 'X2' })).toBe(1);
    expect(counts({ residue_index: 3 * 20 + 1 })).toBe(4);
  });
  it('an empty insertion code is a residue without one (as in MolScript), a given one is that insertion code', () => {
    expect(counts({ pdbx_PDB_ins_code: 'A' })).toBe(3 * 2 * 6);
    expect(counts({ pdbx_PDB_ins_code: '' })).toBe(3 * 18 * 6 + nLigandAtoms);
    expect(counts({ label_asym_id: 'A', pdbx_PDB_ins_code: '' })).toBe(18 * 6);
    // residue 7 of chain A (auth_seq_id 17) has the insertion code A: asking for it with the code finds it, asking for it without does not
    expect(counts({ auth_asym_id: 'x0', auth_seq_id: 17 })).toBe(6);
    expect(counts({ auth_asym_id: 'x0', auth_seq_id: 17, pdbx_PDB_ins_code: 'A' })).toBe(6);
    expect(counts({ auth_asym_id: 'x0', auth_seq_id: 17, pdbx_PDB_ins_code: '' })).toBe(0);
    // residue 8 (auth_seq_id 18) has none
    expect(counts({ auth_asym_id: 'x0', auth_seq_id: 18, pdbx_PDB_ins_code: '' })).toBe(6);
    expect(counts({ auth_asym_id: 'x0', auth_seq_id: 18, pdbx_PDB_ins_code: 'A' })).toBe(0);
  });
});
