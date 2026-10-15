/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Validation's cases, for their heap guard and its quick twin (validate.heap.ts, validate.quick.ts): a table cell checked
// in every column form (a step is one cell, `scanCells` over one position of a column), and the validation of the coded
// genome table (a step is one `validate` of a document of `n` rows; the guard compares a quarter of 3M rows with all
// of them). Columns are built when a case asks, a million cells each.

import { validate } from '../../src/format/index.js';
import { CellKind, CellScan, type Col, ColumnForm, scanCells } from '../../src/format/tables/cells.js';
import { kernel } from './fixtures.js';
import { CHROMOSOMES, codedDoc, freshParams } from './table-fixtures.js';

/** The cells of each column: four times the guard's steps in its shorter run, the index its longer run reaches. */
export const LEN = 1_000_000;

const col = (form: number, data: ArrayLike<unknown>, dict: readonly unknown[] | null = null, name = 'c'): Col => ({
  name,
  form,
  data,
  dict,
  length: data.length,
});
const ints = () => Int32Array.from({ length: LEN }, (_, i) => i);
const doubles = () => Array.from({ length: LEN }, (_, i) => (i % 100) / 100 + 0.001);

export interface CellColumn {
  name: string;
  /** The column and the kind of cell its check reads. */
  make: () => { col: Col; kind: number };
}

export const CELL_COLUMNS: CellColumn[] = [
  { name: 'typed indices', make: () => ({ col: col(ColumnForm.typed, ints()), kind: CellKind.index }) },
  {
    name: 'typed values',
    make: () => ({ col: col(ColumnForm.typed, Float64Array.from(doubles())), kind: CellKind.value }),
  },
  {
    name: 'coded refs',
    make: () => ({
      col: col(
        ColumnForm.dict,
        Uint8Array.from(ints(), (i) => i % 24),
        CHROMOSOMES,
      ),
      kind: CellKind.ref,
    }),
  },
  {
    name: 'string refs',
    make: () => ({
      col: col(
        ColumnForm.array,
        Array.from(ints(), (i) => CHROMOSOMES[i % 24]),
      ),
      kind: CellKind.ref,
    }),
  },
  {
    name: 'label strings',
    make: () => ({
      col: col(
        ColumnForm.array,
        Array.from(ints(), (i) => `${i % 900}A`),
      ),
      kind: CellKind.label,
    }),
  },
  { name: 'an array of doubles', make: () => ({ col: col(ColumnForm.array, doubles()), kind: CellKind.value }) },
  {
    name: 'rows',
    make: () => ({
      col: col(
        ColumnForm.rows,
        Array.from(ints(), (i) => ({ start: i, label: 'x' })),
        null,
        'start',
      ),
      kind: CellKind.index,
    }),
  },
];

/** A column's step, one cell checked, and the scan that counts the invalid cells. */
export function cellStep({ col, kind }: { col: Col; kind: number }) {
  const s = new CellScan();
  s.reset(col);
  return { step: (i: number) => scanCells(col, kind, i % LEN, (i % LEN) + 1, s), scan: s };
}

/**
 * The coded genome table of `n` rows, validated once so its first-run compiles are behind it. A step is one validation
 * of the document with its table's checks run again (`freshParams`).
 */
export function genome(n: number) {
  const reg = kernel();
  const d = codedDoc(n);
  validate(d, reg);
  return {
    step: () => {
      validate(freshParams(d), reg);
    },
  };
}
