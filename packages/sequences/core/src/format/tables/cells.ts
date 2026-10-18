/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The per-cell part of a table check (seqspec:format#annotation-tables, seqspec:format#validation): a column read cell
// by cell with nothing allocated per cell. Failures are counted per column, to be reported once with their count; a ref
// column's values are collected per distinct value, each to become one `RefUse` with its count; a dictionary column
// counts rows per entry, so each entry is checked once. Not exported from the package: `table.ts` reads it, and so do
// the heap guards.

import type { Geometry } from '@molstar/seq-core/model';
import { isPlainObject } from '../schema/field-schema.js';

/**
 * How a column is held (seqspec:format#positions): one key of each row object, an array, a typed array, coded, or JSON.
 */
export const ColumnForm = { rows: 0, array: 1, typed: 2, dict: 3, json: 4 } as const;

export interface Col {
  readonly name: string;
  readonly form: number;
  /** The row objects (rows form), the values, the codes (dict) or the JSON values. */
  readonly data: ArrayLike<unknown>;
  /** The dictionary of a coded column; null otherwise. */
  readonly dict: readonly unknown[] | null;
  readonly length: number;
}

/** What a cell of a column the geometry reads must be. */
export const CellKind = { index: 0, label: 1, ref: 2, value: 3, states: 4, direction: 5, str: 6, cells: 7 } as const;
const INDEX = new Set(['start', 'stop', 'boundary', 'index', 'b_start', 'b_stop', 'row']);
const LABEL = new Set(['beg', 'end', 'pos', 'b_beg', 'b_end']);

export function cellKind(name: string, geometry: Geometry): number {
  if (INDEX.has(name)) return CellKind.index;
  if (LABEL.has(name)) return CellKind.label;
  if (name === 'target' || name === 'b_target') return CellKind.ref;
  if (name === 'value')
    return geometry === 'matrix' ? CellKind.states : geometry === 'grid' ? CellKind.cells : CellKind.value;
  return name === 'direction' ? CellKind.direction : CellKind.str;
}

/** A cell's value; null for a missing one (a key a row lacks, NaN in a typed column, a null or invalid code). */
export function cellAt(c: Col, i: number): unknown {
  const d = c.data;
  if (c.form === ColumnForm.rows) {
    const row = d[i];
    return isPlainObject(row) ? (row[c.name] ?? null) : null;
  }
  if (c.form === ColumnForm.dict) {
    const code = d[i];
    return typeof code === 'number' && Number.isInteger(code) && code >= 0 ? (c.dict?.[code] ?? null) : null;
  }
  const x = d[i];
  return x === undefined || (c.form === ColumnForm.typed && Number.isNaN(x)) ? null : x;
}

const NO_SPACE = /^\S+$/;

/** Whether a non-null value fits a cell kind (the shared rule; the scan inlines its own for numbers). */
export function cellOk(kind: number, x: unknown): boolean {
  switch (kind) {
    case CellKind.index:
      return Number.isInteger(x) && (x as number) >= 0;
    case CellKind.label:
      return Number.isInteger(x) || (typeof x === 'string' && NO_SPACE.test(x));
    case CellKind.value:
      return typeof x === 'number' && Number.isFinite(x);
    case CellKind.states:
    case CellKind.cells: {
      // A grid's cells may be a Float32Array or a Float64Array, NaN missing (seqspec:format#annotation-tables).
      if (kind === CellKind.cells && (x instanceof Float32Array || x instanceof Float64Array)) return true;
      if (!Array.isArray(x)) return false;
      for (let k = 0; k < x.length; k++) {
        const y = x[k];
        if (y !== null && !(typeof y === 'number' && Number.isFinite(y))) return false;
      }
      return true;
    }
    case CellKind.direction:
      return x === 'ab' || x === 'ba';
    default:
      return typeof x === 'string';
  }
}

/** One column's scan, reset per column and reused (validation is synchronous): nothing is allocated per cell. */
export class CellScan {
  /** Cells that fail their kind, and the first one's row. */
  bad = 0;
  firstBad = -1;
  /** Codes neither null nor an index into the dictionary, and the first one's row. */
  badCodes = 0;
  firstBadCode = -1;
  /** A ref column's distinct values (not coded), each to its slot in `counts` and `firsts`. */
  readonly distinct = new Map<string, number>();
  readonly counts: number[] = [];
  readonly firsts: number[] = [];
  /** A coded column's rows per dictionary entry. */
  perEntry: Int32Array = new Int32Array(0);

  reset(c: Col): void {
    this.bad = this.badCodes = 0;
    this.firstBad = this.firstBadCode = -1;
    this.distinct.clear();
    this.counts.length = this.firsts.length = 0;
    const n = c.dict?.length ?? 0;
    if (this.perEntry.length < n) this.perEntry = new Int32Array(n);
    else this.perEntry.fill(0, 0, n);
  }
}

/** Scans rows `from..to` of a column for a kind into `s`. */
export function scanCells(c: Col, kind: number, from: number, to: number, s: CellScan): void {
  const d = c.data;
  if (c.form === ColumnForm.dict) {
    const n = c.dict?.length ?? 0;
    for (let i = from; i < to; i++) {
      const code = d[i] as number | null;
      if (code === null || Number.isNaN(code)) continue;
      if (Number.isInteger(code) && code >= 0 && code < n) s.perEntry[code]++;
      else if (s.badCodes++ === 0) s.firstBadCode = i;
    }
    return;
  }
  if (c.form === ColumnForm.typed) {
    // Numbers only, read as doubles and never passed on, so none is boxed.
    const a = d as Float64Array;
    for (let i = from; i < to; i++) {
      const x = a[i];
      if (Number.isNaN(x)) continue;
      const ok =
        kind === CellKind.value
          ? Number.isFinite(x)
          : (kind === CellKind.index && Number.isInteger(x) && x >= 0) ||
            (kind === CellKind.label && Number.isInteger(x));
      if (!ok && s.bad++ === 0) s.firstBad = i;
    }
    return;
  }
  for (let i = from; i < to; i++) {
    const x = cellAt(c, i);
    if (x === null) continue;
    if (kind === CellKind.ref && typeof x === 'string') {
      const slot = s.distinct.get(x);
      if (slot !== undefined) s.counts[slot]++;
      else {
        s.distinct.set(x, s.counts.length);
        s.counts.push(1);
        s.firsts.push(i);
      }
    } else if (!cellOk(kind, x) && s.bad++ === 0) s.firstBad = i;
  }
}
