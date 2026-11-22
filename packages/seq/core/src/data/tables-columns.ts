/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  CategoryColumn,
  type Column,
  categoryColumn,
  categoryFromCodes,
  columnConstructor,
  columnType,
  columnValues,
  type NumericColumn,
  type NumericColumnType,
} from '@molstar/core/repr/tree';

// Columns as the tables' cuts and joins handle them (_spec/README.md "Data", Tables): rows gathered, parts
// concatenated, the checks' small predicates.

export function isInt(x: unknown): boolean {
  return typeof x === 'number' && Number.isInteger(x);
}

export function sameLength(n: number, ...cols: (Column | undefined)[]): boolean {
  for (const c of cols) if (c !== undefined && c.length !== n) return false;
  return true;
}

export function isObject(x: unknown): x is Record<string, unknown> {
  return typeof x === 'object' && x !== null;
}

/** A category's value at `i` (`null` missing). */
export function categoryOf(c: CategoryColumn, i: number): string | null {
  return c.get(i);
}

/** The rows of a column at the given indices, in order. */
export function gather<C extends Column>(c: C, rows: readonly number[]): C {
  if (c instanceof CategoryColumn) {
    const codes = new Int32Array(rows.length);
    for (let i = 0; i < rows.length; i++) codes[i] = c.codes[rows[i]];
    return categoryFromCodes(c.dict, codes) as C;
  }
  const out = new (columnConstructor(columnType(c) as NumericColumnType))(rows.length) as NumericColumn;
  for (let i = 0; i < rows.length; i++) out[i] = (c as NumericColumn)[rows[i]];
  return out as C;
}

/** Columns of one type joined in order. */
export function concatColumns<C extends Column>(cols: readonly C[]): C {
  if (cols.length === 1) return cols[0];
  if (cols[0] instanceof CategoryColumn) {
    const values: (string | null)[] = [];
    for (const c of cols) for (const v of columnValues(c)) values.push(v as string | null);
    return categoryColumn(values) as C;
  }
  let n = 0;
  for (const c of cols) n += c.length;
  const out = new (columnConstructor(columnType(cols[0]) as NumericColumnType))(n) as NumericColumn;
  let at = 0;
  for (const c of cols) {
    out.set(c as NumericColumn, at);
    at += c.length;
  }
  return out as C;
}
