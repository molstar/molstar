/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  CategoryColumn,
  type Column,
  categoryFromCodes,
  columnConstructor,
  columnType,
  isColumn,
  type NumericColumn,
  type NumericColumnType,
} from '@molstar/core/repr/tree';
import type { DynamicNode, Shape, Table, Tables } from './protocol.js';

// Tables of the protocol's shapes (_spec/README.md "Data", Tables): their rows, their bytes, their check, and the
// slice of a table to a region (a `seq` source answers a tile from a larger region).

/** The rows (positions, bins or ranges) a table holds. */
export function tableRows(t: Table): number {
  switch (t.form) {
    case 'letters':
      return t.letters.length;
    case 'rows':
      return t.start.length;
    case 'density':
      return t.bin.length;
    case 'positions':
      return t.value.length;
    case 'bins':
      return t.count.length;
  }
}

/** Whether an answer holds nothing: no rows in any table and no nodes (an `empty` entry). */
export function isEmptyAnswer(tables: Tables, nodes: readonly DynamicNode[] | undefined): boolean {
  if (nodes?.length) return false;
  for (const t of [tables.letters, tables.ranges, tables.values]) if (t && tableRows(t) > 0) return false;
  return true;
}

function columnBytes(c: Column): number {
  if (c instanceof CategoryColumn) {
    let n = c.codes.byteLength;
    for (const s of c.dict) n += 2 * s.length;
    return n;
  }
  return c.byteLength;
}

/** A held node's cost in the byte bound: nodes are small and few against tables; a flat estimate. */
export const NodeBytes = 256;

/** What an answer costs the store's byte bound: its columns' bytes, a string's two bytes a letter, nodes flat. */
export function answerBytes(tables: Tables, nodes: readonly DynamicNode[] | undefined): number {
  let n = (nodes?.length ?? 0) * NodeBytes;
  const { letters: l, ranges: r, values: v } = tables;
  if (l) n += typeof l.letters === 'string' ? 2 * l.letters.length : l.letters.byteLength;
  if (r?.form === 'rows') {
    n += r.start.byteLength + r.end.byteLength + (r.strand?.byteLength ?? 0) + (r.id ? columnBytes(r.id) : 0);
    if (r.fields) for (const c of Object.values(r.fields)) n += columnBytes(c);
  } else if (r) n += r.bin.byteLength + r.count.byteLength + (r.category ? columnBytes(r.category) : 0);
  if (v?.form === 'positions') n += v.value.byteLength;
  else if (v) n += v.min.byteLength + v.max.byteLength + v.mean.byteLength + v.count.byteLength;
  return n;
}

const Forms: Record<Shape, readonly string[]> = {
  letters: ['letters'],
  ranges: ['rows', 'density'],
  values: ['positions', 'bins'],
};

function isInt(x: unknown): boolean {
  return typeof x === 'number' && Number.isInteger(x);
}

function sameLength(n: number, ...cols: (Column | undefined)[]): boolean {
  for (const c of cols) if (c !== undefined && c.length !== n) return false;
  return true;
}

/**
 * What is wrong with an answer's tables, as a message, or `undefined` when each is of its shape's forms with the
 * protocol's column types and equal lengths.
 */
export function checkTables(tables: Tables): string | undefined {
  if (typeof tables !== 'object' || tables === null) return 'the tables are not a record';
  for (const shape of Object.keys(tables)) {
    const t = (tables as Record<string, Table | undefined>)[shape];
    if (t === undefined) continue;
    const forms = Object.hasOwn(Forms, shape) ? Forms[shape as Shape] : undefined;
    if (!forms) return 'a table of a shape other than letters, ranges and values';
    if (typeof t !== 'object' || t === null || !forms.includes(t.form))
      return `the ${shape} table has no form of its shape`;
    const bad = checkTable(t);
    if (bad) return `the ${shape} table: ${bad}`;
  }
  return undefined;
}

function checkTable(t: Table): string | undefined {
  switch (t.form) {
    case 'letters':
      if (!isInt(t.start)) return '`start` is not an integer';
      if (typeof t.letters === 'string') return undefined;
      if (!(t.letters instanceof Uint8Array)) return '`letters` is neither a string nor uint8 codes';
      return typeof t.alphabet === 'string' ? undefined : 'codes without an `alphabet`';
    case 'rows': {
      if (!(t.start instanceof Int32Array) || !(t.end instanceof Int32Array)) return '`start` and `end` are not int32';
      if (t.strand !== undefined && !(t.strand instanceof Int8Array)) return '`strand` is not int8';
      if (t.id !== undefined && !(t.id instanceof CategoryColumn)) return '`id` is not a category';
      const n = t.start.length;
      if (!sameLength(n, t.end, t.strand, t.id)) return 'columns of unequal lengths';
      for (const [name, c] of Object.entries(t.fields ?? {})) {
        if (!isColumn(c)) return `the field ${name} is not a column`;
        if (c.length !== n) return 'columns of unequal lengths';
      }
      return undefined;
    }
    case 'density':
      if (!(t.bin instanceof Int32Array) || !(t.count instanceof Float32Array))
        return '`bin` and `count` are not int32 and float32';
      if (t.category !== undefined && !(t.category instanceof CategoryColumn)) return '`category` is not a category';
      return sameLength(t.bin.length, t.count, t.category) ? undefined : 'columns of unequal lengths';
    case 'positions':
      if (!isInt(t.start)) return '`start` is not an integer';
      return t.value instanceof Float32Array ? undefined : '`value` is not float32';
    case 'bins':
      if (!isInt(t.start)) return '`start` is not an integer';
      for (const c of [t.min, t.max, t.mean, t.count])
        if (!(c instanceof Float32Array)) return 'a summary is not float32';
      return sameLength(t.count.length, t.min, t.max, t.mean) ? undefined : 'columns of unequal lengths';
  }
}

function gather<C extends Column>(c: C, rows: readonly number[]): C {
  if (c instanceof CategoryColumn) {
    const codes = new Int32Array(rows.length);
    for (let i = 0; i < rows.length; i++) codes[i] = c.codes[rows[i]];
    return categoryFromCodes(c.dict, codes) as C;
  }
  const out = new (columnConstructor(columnType(c) as NumericColumnType))(rows.length) as NumericColumn;
  for (let i = 0; i < rows.length; i++) out[i] = (c as NumericColumn)[rows[i]];
  return out as C;
}

function clip(lo: number, hi: number, n: number): [number, number] {
  return [Math.max(0, Math.min(n, lo)), Math.max(0, Math.min(n, hi))];
}

/**
 * `t` cut to the positions `[start, end)` at `level` (the table's own level): letters, positions and bins by offset,
 * ranges rows by overlap (a range over a tile's edge is in each tile it touches), density bins by index. The table
 * itself when nothing is cut.
 */
export function sliceTable<T extends Table>(t: T, start: number, end: number, level: number): T {
  switch (t.form) {
    case 'letters':
    case 'positions': {
      const n = t.form === 'letters' ? t.letters.length : t.value.length;
      const [a, b] = clip(start - t.start, end - t.start, n);
      if (a === 0 && b === n) return t;
      if (t.form === 'letters') return { ...t, start: t.start + a, letters: t.letters.slice(a, b) } as T;
      return { ...t, start: t.start + a, value: t.value.slice(a, b) } as T;
    }
    case 'bins': {
      const n = t.count.length;
      const [a, b] = clip(Math.floor((start - t.start) / level), Math.ceil((end - t.start) / level), n);
      if (a === 0 && b === n) return t;
      return {
        ...t,
        start: t.start + a * level,
        min: t.min.slice(a, b),
        max: t.max.slice(a, b),
        mean: t.mean.slice(a, b),
        count: t.count.slice(a, b),
      } as T;
    }
    case 'rows': {
      const rows: number[] = [];
      for (let i = 0; i < t.start.length; i++) if (t.start[i] < end && t.end[i] > start) rows.push(i);
      if (rows.length === t.start.length) return t;
      const fields: Record<string, Column> | undefined = t.fields ? {} : undefined;
      if (t.fields) for (const [k, c] of Object.entries(t.fields)) fields![k] = gather(c, rows);
      return {
        form: 'rows',
        start: gather(t.start, rows),
        end: gather(t.end, rows),
        ...(t.strand && { strand: gather(t.strand, rows) }),
        ...(t.id && { id: gather(t.id, rows) }),
        ...(fields && { fields }),
      } as T;
    }
    case 'density': {
      const lo = Math.floor(start / level);
      const hi = Math.ceil(end / level);
      const rows: number[] = [];
      for (let i = 0; i < t.bin.length; i++) if (t.bin[i] >= lo && t.bin[i] < hi) rows.push(i);
      if (rows.length === t.bin.length) return t;
      return {
        form: 'density',
        bin: gather(t.bin, rows),
        count: gather(t.count, rows),
        ...(t.category && { category: gather(t.category, rows) }),
      } as T;
    }
  }
}
