/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, type Column, columnType, isColumn } from '@molstar/core/repr/tree';
import type { RangesField } from './protocol.js';
import { concatColumns, gather, isObject } from './tables-columns.js';

// The row tables' forms as one schema (_spec/README.md "Data", Tables): a form whose rows are columns of one length
// (a ranges list, a density, pairs and bin pairs, a matrix's items, a sparse grid's cells, keyed rows) is its columns'
// names and types, those that come in twos, its fields and its cells; its rows, bytes, checks, gathers and joins are
// read from that, one routine for every form, so a form is declared once. Forms by offset (letters, values, cells
// codes and majority, a dense grid, a matrix's summary) stay written out.

/** A column's type in a form. */
export type ColKind = 'i32' | 'f32' | 'i8' | 'u8' | 'cat';

/** A form's column: its name, its type, whether it is required. */
export type ColSpec = readonly [name: string, kind: ColKind, required?: 1];

/** A row table's form. */
export interface RowForm {
  /** Its row columns, the first required: the rows are its length. */
  readonly cols: readonly ColSpec[];
  /** Columns that come together or not at all (an end's system and the other's). */
  readonly both?: readonly (readonly [string, string])[];
  /** Named columns besides (`RangesField`s). */
  readonly fields?: 1;
  /** Cells: columns of `states` values a row, exactly one present (a matrix's values or codes). */
  readonly cells?: readonly ColSpec[];
  /** Numbers carried as they are (equal across joined parts). */
  readonly scalars?: readonly string[];
}

type Row = Record<string, unknown>;

const Is: Record<ColKind, (c: unknown) => boolean> = {
  i32: (c) => c instanceof Int32Array,
  f32: (c) => c instanceof Float32Array,
  i8: (c) => c instanceof Int8Array,
  u8: (c) => c instanceof Uint8Array,
  cat: (c) => c instanceof CategoryColumn,
};

const KindName: Record<ColKind, string> = { i32: 'int32', f32: 'float32', i8: 'int8', u8: 'uint8', cat: 'a category' };

export function columnBytes(c: Column | undefined): number {
  if (c === undefined) return 0;
  if (c instanceof CategoryColumn) {
    let n = c.codes.byteLength;
    for (const s of c.dict) n += 2 * s.length;
    return n;
  }
  return c.byteLength;
}

/** A row table's rows: its first column's length. */
export function formRows(t: object, f: RowForm): number {
  return ((t as Row)[f.cols[0][0]] as Column).length;
}

function fieldsOf(t: object): readonly RangesField[] | undefined {
  return (t as { fields?: readonly RangesField[] }).fields;
}

function cellsOf(t: object, f: RowForm): Column | undefined {
  for (const [name] of f.cells ?? []) {
    const c = (t as Row)[name] as Column | undefined;
    if (c) return c;
  }
  return undefined;
}

/** The columns a row table holds (its cells a column a state): what a block's values are counted by. */
export function formColumnCount(t: object, f: RowForm): number {
  let n = 0;
  for (const [name] of f.cols) if ((t as Row)[name] !== undefined) n++;
  if (f.cells) n += Math.max(1, (t as { states: number }).states);
  return n + (fieldsOf(t)?.length ?? 0);
}

/** What a row table costs the byte bound: its columns' bytes. */
export function formBytes(t: object, f: RowForm): number {
  let n = columnBytes(cellsOf(t, f));
  for (const [name] of f.cols) n += columnBytes((t as Row)[name] as Column | undefined);
  for (const x of fieldsOf(t) ?? []) n += columnBytes(x.column);
  return n;
}

/** What is wrong with a row table's columns, as a message, or `undefined`. */
export function checkForm(t: object, f: RowForm): string | undefined {
  const r = t as Row;
  for (const [name, kind, required] of f.cols) {
    const c = r[name];
    if (c === undefined ? required : !Is[kind](c)) return `\`${name}\` is not ${KindName[kind]}`;
  }
  for (const [a, b] of f.both ?? [])
    if ((r[a] === undefined) !== (r[b] === undefined)) return `\`${a}\` without \`${b}\``;
  const n = formRows(t, f);
  for (const [name] of f.cols) {
    const c = r[name] as Column | undefined;
    if (c !== undefined && c.length !== n) return 'columns of unequal lengths';
  }
  if (f.cells) {
    const states = r.states;
    if (typeof states !== 'number' || !Number.isInteger(states) || states < 0)
      return '`states` is not an integer of at least 0';
    const present = f.cells.filter(([name]) => r[name] !== undefined);
    if (present.length !== 1) return `not one of ${f.cells.map(([name]) => `\`${name}\``).join(' and ')}`;
    const [name, kind] = present[0];
    if (!Is[kind](r[name])) return `\`${name}\` is not ${KindName[kind]}`;
    if ((r[name] as Column).length !== n * states) return 'cells are not items × states';
  }
  const fields = r.fields;
  if (fields === undefined) return undefined;
  if (!f.fields || !Array.isArray(fields)) return '`fields` is not a list';
  for (const x of fields) {
    if (!isObject(x) || typeof x.name !== 'string' || !isColumn(x.column as Column))
      return 'a field is not a named column';
    if ((x.column as Column).length !== n) return 'columns of unequal lengths';
  }
  return undefined;
}

/** Blocks of `size` values at the given block indices, in order (a matrix's states per item). */
function gatherBlocks<A extends Float32Array | Uint8Array>(a: A, rows: readonly number[], size: number): A {
  const out = new (a.constructor as new (n: number) => A)(rows.length * size);
  for (let i = 0; i < rows.length; i++) out.set(a.subarray(rows[i] * size, (rows[i] + 1) * size), i * size);
  return out;
}

/** A row table's rows at the given indices, in order: every column present gathered, its scalars kept. */
export function gatherForm<T extends object>(t: T, f: RowForm, rows: readonly number[]): T {
  const r = t as Row;
  const out: Row = { form: r.form };
  for (const s of f.scalars ?? []) if (r[s] !== undefined) out[s] = r[s];
  for (const [name] of f.cols) if (r[name] !== undefined) out[name] = gather(r[name] as Column, rows);
  for (const [name] of f.cells ?? [])
    if (r[name] !== undefined) out[name] = gatherBlocks(r[name] as Float32Array, rows, r.states as number);
  const fields = fieldsOf(t);
  if (fields) out.fields = fields.map((x) => ({ name: x.name, column: gather(x.column, rows) }));
  return out as T;
}

/** What columns a row table carries: equal across parts that join. */
function shapeOf(t: object, f: RowForm): string {
  const r = t as Row;
  let s = '';
  for (const [name] of f.cols) s += r[name] === undefined ? '0' : '1';
  for (const [name] of f.cells ?? []) s += r[name] === undefined ? '0' : '1';
  for (const x of f.scalars ?? []) s += `|${String(r[x])}`;
  for (const x of fieldsOf(t) ?? []) s += `|${x.name}:${columnType(x.column)}`;
  return s;
}

/** Row tables of one form joined in order (each part's rows after the last's), or a message when columns differ. */
export function joinForm<T extends object>(parts: readonly T[], f: RowForm): T | string {
  const first = parts[0] as Row;
  const shape = shapeOf(first, f);
  if (parts.some((p) => shapeOf(p, f) !== shape)) return 'the parts carry different columns';
  const out: Row = { form: first.form };
  for (const s of f.scalars ?? []) if (first[s] !== undefined) out[s] = first[s];
  for (const [name] of [...f.cols, ...(f.cells ?? [])])
    if (first[name] !== undefined) out[name] = concatColumns(parts.map((p) => (p as Row)[name] as Column));
  const fields = fieldsOf(first);
  if (fields)
    out.fields = fields.map((x, j) => ({
      name: x.name,
      column: concatColumns(parts.map((p) => fieldsOf(p)![j].column)),
    }));
  return out as T;
}

/** Rows `[a, b)` by offset a table starts at `start` holds outside a region at `step` positions a row, or a message. */
export function offsetOutside(
  start: number,
  n: number,
  step: number,
  r: { start: number; end: number },
): string | undefined {
  if (step > 1 && start % step !== 0) return 'not aligned to the level';
  if (n > 0 && (start < Math.floor(r.start / step) * step || start + n * step > Math.ceil(r.end / step) * step))
    return 'outside the region';
  return undefined;
}

/** A list's rows `[s, e)` that end before they start or miss the region, as a message (`what` names a row). */
export function rangesOutside(
  starts: Int32Array,
  ends: Int32Array,
  r: { start: number; end: number },
  maxRows: number,
  what: string,
): string | undefined {
  const n = starts.length;
  if (n > maxRows) return `more than ${maxRows} rows`;
  for (let i = 0; i < n; i++) {
    const s = starts[i];
    const e = ends[i];
    if (s > e) return `${what} ends before it starts`;
    if (s >= r.end || Math.max(e, s + 1) <= r.start) return `${what} outside the region`;
  }
  return undefined;
}

/** Strands outside −1, 0, 1, as a message. */
export function strandsOutside(...cols: (Int8Array | undefined)[]): string | undefined {
  for (const s of cols)
    if (s) for (let i = 0; i < s.length; i++) if (s[i] < -1 || s[i] > 1) return 'a strand not −1, 0, 1';
  return undefined;
}

// ---- the forms ------------------------------------------------------------------------------------------------------

const Fields = 1 as const;

export const RangesListForm: RowForm = {
  cols: [
    ['start', 'i32', 1],
    ['end', 'i32', 1],
    ['strand', 'i8'],
    ['id', 'cat'],
    ['row', 'i32'],
    ['itemStart', 'i32'],
    ['itemEnd', 'i32'],
  ],
  both: [['itemStart', 'itemEnd']],
  fields: Fields,
};

export const DensityForm: RowForm = {
  cols: [
    ['bin', 'i32', 1],
    ['count', 'f32', 1],
    ['category', 'cat'],
  ],
};

export const PairsListForm: RowForm = {
  cols: [
    ['aStart', 'i32', 1],
    ['aEnd', 'i32', 1],
    ['bStart', 'i32', 1],
    ['bEnd', 'i32', 1],
    ['aSystem', 'cat'],
    ['bSystem', 'cat'],
    ['aStrand', 'i8'],
    ['bStrand', 'i8'],
    ['type', 'cat'],
    ['id', 'cat'],
    ['row', 'i32'],
  ],
  both: [['aSystem', 'bSystem']],
  fields: Fields,
};

export const PairsBinsForm: RowForm = {
  cols: [
    ['aBin', 'i32', 1],
    ['bBin', 'i32', 1],
    ['count', 'f32', 1],
    ['aSystem', 'cat'],
    ['bSystem', 'cat'],
    ['strand', 'i8'],
    ['category', 'cat'],
    ['sum', 'f32'],
  ],
  both: [['aSystem', 'bSystem']],
};

export const MatrixItemsForm: RowForm = {
  cols: [
    ['start', 'i32', 1],
    ['end', 'i32', 1],
    ['id', 'cat'],
    ['row', 'i32'],
  ],
  cells: [
    ['value', 'f32'],
    ['codes', 'u8'],
  ],
  scalars: ['states'],
};

export const SparseForm: RowForm = {
  cols: [
    ['xBin', 'i32', 1],
    ['yBin', 'i32', 1],
    ['value', 'f32', 1],
  ],
  scalars: ['fill'],
};

export const KeyedForm: RowForm = { cols: [['key', 'cat', 1]], fields: Fields };
