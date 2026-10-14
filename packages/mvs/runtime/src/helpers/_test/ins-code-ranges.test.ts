/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Structure } from '@molstar/model/model/structure';
import type { SchemaItem } from '@molstar/model/model/structure/structure/element/schema';
import { ElementRanges } from '@molstar/model/model/structure/structure/element/schema-selection/element-ranges';
import { IndicesAndSortings } from '@molstar/model/model/structure/structure/element/schema-selection/indexing';
import { substructureFromRows } from '@molstar/model/model/structure/structure/element/schema-selection/row-structure';
import {
  getAtomRangesForRow,
  getAtomRangesForRows,
} from '@molstar/model/model/structure/structure/element/schema-selection/selections';
import { TestOperator, createModel, insCodeRangeCategories, makeRows } from './annotation-rows-fixture.js';
import { createRowOracle } from './row-oracle.js';

/** Ranges of `auth_seq_id` with insertion codes against the definition written on the residues of each chain (`row-oracle.ts`). */
describe.each([
  ['insertion series', { nChains: 12, nResidues: 30, insertionSeries: true }],
  [
    'descending label_seq_id, mixed comp_ids, auth comp_ids',
    {
      nChains: 12,
      nResidues: 30,
      insertionSeries: true,
      descendingSeqIds: true,
      mixedCompIds: true,
      differentAuthComp: true,
    },
  ],
  [
    'ligands without label_seq_id, waters, repeated atom ids',
    { nChains: 12, nResidues: 30, insertionSeries: true, ligands: 3, waters: 14, atomIdWrap: 150 },
  ],
] as const)('ranges with insertion codes vs per-atom brute force (%s)', async (_label, dims) => {
  const model = await createModel(dims);
  const nAtoms = model.atomicHierarchy.atoms._rowCount;
  const indices = IndicesAndSortings.get(model);
  const structure = Structure.ofModel(model);
  const oracle = createRowOracle(model);
  const groups = makeRows<SchemaItem>(insCodeRangeCategories(dims, nAtoms), 7);

  const atomsInRanges = (ranges: ElementRanges | undefined) => {
    const atoms: number[] = [];
    if (ranges)
      ElementRanges.foreach(ranges, (from, to) => {
        for (let a = from; a < to; a++) atoms.push(a);
      });
    return atoms;
  };
  const atomsOfStructure = (s: Structure) => s.units.flatMap((u) => Array.from(u.elements));

  it('the fixture has insertion series (residues which share auth_seq_id)', () => {
    const { auth_seq_id, pdbx_PDB_ins_code } = model.atomicHierarchy.residues;
    // the 6th residue of a chain is 16, the 7th is 16A, the 8th is 16B
    const found = [5, 6, 7].map((i) => `${auth_seq_id.value(i)}${pdbx_PDB_ins_code.value(i)}`);
    expect(found).toEqual(['16', '16A', '16B']);
  });

  it('the waters (when there are some) are in the entity of waters and in the auth chain of the first chain, with numbers of its residues', () => {
    const { chains, residues } = model.atomicHierarchy;
    const iWater = Array.from({ length: chains._rowCount }, (_, i) => i).find(
      (i) => chains.label_asym_id.value(i) === 'W',
    );
    if (!('waters' in dims)) return expect(iWater).toBeUndefined();
    expect(iWater).toBeDefined();
    expect(chains.auth_asym_id.value(iWater!)).toBe(chains.auth_asym_id.value(0));
    expect(model.entities.data.type.value(model.entities.getEntityIndex(chains.label_entity_id.value(iWater!)))).toBe(
      'water',
    );
    const numbers = (iChain: number) => {
      const { residueAtomSegments, chainAtomSegments } = model.atomicHierarchy;
      const from = residueAtomSegments.index[chainAtomSegments.offsets[iChain]];
      const to = residueAtomSegments.index[chainAtomSegments.offsets[iChain + 1] - 1];
      return Array.from({ length: to - from + 1 }, (_, k) => residues.auth_seq_id.value(from + k));
    };
    expect(numbers(iWater!).slice(0, 3)).toEqual([11, 12, 13]);
    expect(numbers(0)).toContain(11);
  });

  for (const g of groups) {
    it(`${g.name}: row by row`, () => {
      let nonEmpty = 0;
      for (const row of g.rows) {
        const expected = oracle.atomsOf([row], TestOperator);
        if (expected.length) nonEmpty++;
        expect(atomsInRanges(getAtomRangesForRow(row, model, TestOperator, indices)), JSON.stringify(row)).toEqual(
          expected,
        );
      }
      // (the rows select something, except for the ones which are about nothing on purpose)
      if (g.name === 'insertion code without its number') expect(nonEmpty).toBe(0);
      else expect(nonEmpty).toBeGreaterThan(0);
    });
  }

  it('all rows together: the union of ranges and the substructure', () => {
    const rows = groups.flatMap((g) => g.rows);
    const expected = oracle.atomsOf(rows, TestOperator);
    expect(expected.length).toBeGreaterThan(0);
    expect(atomsInRanges(getAtomRangesForRows(rows, model, TestOperator, indices))).toEqual(expected);
    // (the instance_id rows have another instance than the structure's unit, which only the structure path tells apart)
    const withoutInstances = rows.filter((r) => r.instance_id === undefined);
    expect(atomsOfStructure(substructureFromRows(structure, withoutInstances))).toEqual(
      oracle.atomsOf(withoutInstances, { instanceId: '1_555', name: '1_555' }),
    );
  });
});
