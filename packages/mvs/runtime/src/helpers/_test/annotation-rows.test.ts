/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { ElementIndex } from '@molstar/model/model/structure';
import { MVSAnnotation } from '../../components/annotation-prop.js';
import { ElementRanges } from '../element-ranges.js';
import { IndicesAndSortings } from '../indexing.js';
import type { MVSAnnotationRow } from '../schemas.js';
import { atomQualifies, getAtomRangesForRow, getAtomRangesForRows } from '../selections.js';
import { TestOperator, combinedCategories, createModel, makeRows, rowCategories } from './annotation-rows-fixture.js';

/** The per-atom definition of "atom satisfies row" (`atomQualifies`, a linear scan) is the oracle for the index-driven implementation. */
describe.each([
  [
    'single comp_id per residue',
    { nChains: 12, nResidues: 30, insertionCodes: true, descendingSeqIds: true, differentAuthComp: true },
  ],
  [
    'mixed comp_ids within residues',
    {
      nChains: 12,
      nResidues: 30,
      insertionCodes: true,
      descendingSeqIds: true,
      differentAuthComp: true,
      mixedCompIds: true,
    },
  ],
  [
    'repeated atom ids',
    {
      nChains: 12,
      nResidues: 30,
      insertionCodes: true,
      descendingSeqIds: true,
      differentAuthComp: true,
      atomIdWrap: 150,
    },
  ],
] as const)('annotation rows vs per-atom brute force (%s)', async (_label, dims) => {
  const model = await createModel(dims);
  const nAtoms = model.atomicHierarchy.atoms._rowCount;
  const indices = IndicesAndSortings.get(model);
  const groups = [
    ...makeRows(rowCategories(dims, nAtoms, 1000), 4242),
    ...makeRows(combinedCategories(dims, nAtoms), 99),
  ];
  const allRows: MVSAnnotationRow[] = groups.flatMap((g) => g.rows);

  /** Last row (by index) which the atom satisfies, -1 for none (later rows win). */
  function bruteForceRowForEachAtom(rows: MVSAnnotationRow[]) {
    const out = new Array<number>(nAtoms).fill(-1);
    for (let iAtom = 0; iAtom < nAtoms; iAtom++) {
      for (let iRow = 0; iRow < rows.length; iRow++) {
        const row = rows[iRow];
        // a row for another instance or operator never applies (`atomQualifies` is only the atom-level part)
        if (row.instance_id !== undefined && row.instance_id !== TestOperator.instanceId) continue;
        if (row.operator_name !== undefined && row.operator_name !== TestOperator.name) continue;
        if (atomQualifies(model, iAtom as ElementIndex, row)) out[iAtom] = iRow;
      }
    }
    return out;
  }

  function atomsInRanges(ranges: ElementRanges | undefined) {
    const atoms: number[] = [];
    if (ranges)
      ElementRanges.foreach(ranges, (from, to) => {
        for (let a = from; a < to; a++) atoms.push(a);
      });
    return atoms;
  }

  it('the fixture has enough variety to be meaningful', () => {
    expect(nAtoms).toBe(dims.nChains * dims.nResidues * 6);
    expect(groups.length).toBeGreaterThan(25);
    expect(indices.atomic!.residuesByLabelCompIdIsPure).toBe(!('mixedCompIds' in dims));
    expect(indices.atomic!.residuesByAuthCompIdIsPure).toBe(!('mixedCompIds' in dims));
  });

  for (const g of groups) {
    it(`getRowForEachAtom: ${g.name}`, () => {
      const annot = new MVSAnnotation({ format: 'json', data: g.rows as any }, 'all_atomic', {});
      const got = (annot as any).getRowForEachAtom(model, TestOperator).atoms as ArrayLike<number> | null;
      const expected = bruteForceRowForEachAtom(g.rows);
      expect(Array.from(got ?? new Array(nAtoms).fill(-1))).toEqual(expected);
    });
    it(`getAtomRangesForRow, one row at a time: ${g.name}`, () => {
      for (const row of g.rows.slice(0, 25)) {
        const expected = bruteForceRowForEachAtom([row]).flatMap((v, i) => (v >= 0 ? [i] : []));
        expect(atomsInRanges(getAtomRangesForRow(row, model, TestOperator, indices))).toEqual(expected);
      }
    });
  }

  it('getRowForEachAtom, all categories mixed', () => {
    const annot = new MVSAnnotation({ format: 'json', data: allRows as any }, 'all_atomic', {});
    const got = (annot as any).getRowForEachAtom(model, TestOperator).atoms as ArrayLike<number>;
    expect(Array.from(got)).toEqual(bruteForceRowForEachAtom(allRows));
  });

  it('getAtomRangesForRows is the union over rows, and its result is not clobbered by later calls', () => {
    const half = Math.floor(allRows.length / 2);
    const first = getAtomRangesForRows(allRows.slice(0, half), model, TestOperator, indices);
    const firstAtoms = atomsInRanges(first);
    const second = getAtomRangesForRows(allRows.slice(half), model, TestOperator, indices);
    getAtomRangesForRow(allRows[0], model, TestOperator, indices); // reuses the shared scratch again
    expect(atomsInRanges(first)).toEqual(firstAtoms); // (a cached result must stay intact)
    const expectedFirst = bruteForceRowForEachAtom(allRows.slice(0, half)).flatMap((v, i) => (v >= 0 ? [i] : []));
    const expectedSecond = bruteForceRowForEachAtom(allRows.slice(half)).flatMap((v, i) => (v >= 0 ? [i] : []));
    expect(firstAtoms).toEqual(expectedFirst);
    expect(atomsInRanges(second)).toEqual(expectedSecond);
  });
});
