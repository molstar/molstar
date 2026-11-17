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
  isColumn,
  isNode,
  type NumericColumn,
  type NumericColumnType,
  type Value,
} from '@molstar/core/repr/tree';
import {
  type DynamicNode,
  FailureKinds,
  type LettersTable,
  type RangesDensity,
  type RangesField,
  type RangesList,
  type RegionRequest,
  type Shape,
  type Table,
  type Tables,
  type ValuesBins,
  type ValuesPositions,
} from './protocol.js';

// Tables of the protocol's shapes (_spec/README.md "Data", Tables): their rows, their bytes, their checks (shape, and
// an answer against its request), and cutting and joining tables by region (a `seq` source answers a tile from its
// pairs).

/** The rows (positions, bins or ranges) a table holds. */
export function tableRows(t: Table): number {
  switch (t.form) {
    case 'letters':
      return t.letters.length;
    case 'list':
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

function columnBytes(c: Column | undefined): number {
  if (c === undefined) return 0;
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
  if (r?.form === 'list') {
    n += columnBytes(r.start) + columnBytes(r.end) + columnBytes(r.strand) + columnBytes(r.id) + columnBytes(r.row);
    for (const f of r.fields ?? []) n += columnBytes(f.column);
  } else if (r) n += columnBytes(r.bin) + columnBytes(r.count) + columnBytes(r.category);
  if (v?.form === 'positions') n += v.value.byteLength;
  else if (v) n += v.min.byteLength + v.max.byteLength + v.mean.byteLength + v.count.byteLength;
  return n;
}

const Forms: Record<Shape, readonly string[]> = {
  letters: ['letters'],
  ranges: ['list', 'density'],
  values: ['positions', 'bins'],
};

function isInt(x: unknown): boolean {
  return typeof x === 'number' && Number.isInteger(x);
}

function sameLength(n: number, ...cols: (Column | undefined)[]): boolean {
  for (const c of cols) if (c !== undefined && c.length !== n) return false;
  return true;
}

function isObject(x: unknown): x is Record<string, unknown> {
  return typeof x === 'object' && x !== null;
}

/**
 * What is wrong with tables, as a message, or `undefined` when each is of its shape's forms with the protocol's column
 * types and equal lengths.
 */
export function checkTables(tables: Tables): string | undefined {
  if (!isObject(tables)) return 'the tables are not a record';
  for (const shape of Object.keys(tables)) {
    const t = (tables as Record<string, Table | undefined>)[shape];
    if (t === undefined) continue;
    const forms = Object.hasOwn(Forms, shape) ? Forms[shape as Shape] : undefined;
    if (!forms) return 'a table of a shape other than letters, ranges and values';
    if (!isObject(t) || !forms.includes(t.form)) return `the ${shape} table has no form of its shape`;
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
    case 'list': {
      if (!(t.start instanceof Int32Array) || !(t.end instanceof Int32Array)) return '`start` and `end` are not int32';
      if (t.strand !== undefined && !(t.strand instanceof Int8Array)) return '`strand` is not int8';
      if (t.id !== undefined && !(t.id instanceof CategoryColumn)) return '`id` is not a category';
      if (t.row !== undefined && !(t.row instanceof Int32Array)) return '`row` is not int32';
      const n = t.start.length;
      if (!sameLength(n, t.end, t.strand, t.id, t.row)) return 'columns of unequal lengths';
      if (t.fields !== undefined && !Array.isArray(t.fields)) return '`fields` is not a list';
      for (const f of t.fields ?? []) {
        if (!isObject(f) || typeof f.name !== 'string' || !isColumn(f.column)) return 'a field is not a named column';
        if (f.column.length !== n) return 'columns of unequal lengths';
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

/** Rows a list or density table may hold in one answer, by default. */
export const DefaultMaxRows = 1_000_000;

/**
 * What is wrong with an answer to `request`, as a message, or `undefined`. Beyond the tables' shapes: a record, a
 * failure of a known kind with an issue, nodes with string keys, only the shapes asked, each table within the region
 * (bins aligned to the level, positions at level 1 only, bins above it), ranges with `start ≤ end`, strands −1, 0, 1,
 * codes within their alphabet, at most `maxRows` rows.
 */
export function checkAnswer(answer: unknown, request: RegionRequest, maxRows = DefaultMaxRows): string | undefined {
  if (!isObject(answer)) return 'the answer is not a record';
  if (answer.ok === false) {
    if (!FailureKinds.includes(answer.kind as never)) return 'a failure of an unknown kind';
    const issue = answer.issue;
    if (!isObject(issue) || typeof issue.code !== 'string' || typeof issue.message !== 'string')
      return 'a failure without an issue';
    return undefined;
  }
  if (answer.ok !== true) return '`ok` is neither true nor false';
  if (answer.nodes !== undefined) {
    if (!Array.isArray(answer.nodes)) return '`nodes` is not a list';
    for (const n of answer.nodes)
      if (!isObject(n) || typeof n.key !== 'string' || !isNode(n.node as Value))
        return 'a node is not a key and a node';
  }
  const tables = answer.tables as Tables;
  const bad = checkTables(tables);
  if (bad) return bad;
  for (const shape of Object.keys(tables)) {
    if ((tables as Record<string, unknown>)[shape] === undefined) continue;
    if (!request.shapes.includes(shape as Shape)) return `a ${shape} table, not asked for`;
    const t = tables[shape as Shape]!;
    const out = checkExtent(t, request, maxRows);
    if (out) return `the ${shape} table: ${out}`;
  }
  return undefined;
}

function checkExtent(t: Table, r: RegionRequest, maxRows: number): string | undefined {
  const L = r.level;
  switch (t.form) {
    case 'letters':
    case 'positions': {
      if (t.form === 'positions' && L !== 1) return 'values per position above level 1';
      const n = tableRows(t);
      if (n > 0 && (t.start < r.start || t.start + n > r.end)) return 'outside the region';
      if (t.form === 'letters' && typeof t.letters !== 'string') {
        const size = t.alphabet!.length;
        for (let i = 0; i < n; i++) if (t.letters[i] >= size) return 'a code outside its alphabet';
      }
      return undefined;
    }
    case 'bins': {
      if (L === 1) return 'bins at level 1';
      const n = t.count.length;
      if (t.start % L !== 0) return 'bins not aligned to the level';
      if (n > 0 && (t.start < Math.floor(r.start / L) * L || t.start + n * L > Math.ceil(r.end / L) * L))
        return 'outside the region';
      return undefined;
    }
    case 'list': {
      const n = t.start.length;
      if (n > maxRows) return `more than ${maxRows} rows`;
      for (let i = 0; i < n; i++) {
        const s = t.start[i];
        const e = t.end[i];
        if (s > e) return 'a range ends before it starts';
        if (s >= r.end || Math.max(e, s + 1) <= r.start) return 'a range outside the region';
      }
      if (t.strand) for (let i = 0; i < n; i++) if (t.strand[i] < -1 || t.strand[i] > 1) return 'a strand not −1, 0, 1';
      return undefined;
    }
    case 'density': {
      const n = t.bin.length;
      if (n > maxRows) return `more than ${maxRows} rows`;
      const lo = Math.floor(r.start / L);
      const hi = Math.ceil(r.end / L);
      for (let i = 0; i < n; i++) if (t.bin[i] < lo || t.bin[i] >= hi) return 'a bin outside the region';
      return undefined;
    }
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

function concatColumns<C extends Column>(cols: readonly C[]): C {
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

/** Where a cut keeps a row: the region `[start, end)`, and the part of it `[from, to)` this cut answers for. */
export interface Cut {
  readonly start: number;
  readonly end: number;
  readonly from?: number;
  readonly to?: number;
  /** A list without `id` or `row` gets `row` = `rowBase` + its index, a stable identity across cuts. */
  readonly rowBase?: number;
}

/**
 * `t` cut to a region at `level` (the table's own). Letters, positions and bins by offset; a range or a density bin
 * when it overlaps the region and its first position in the region (its anchor) is in `[from, to)` (the whole
 * region by default), so joined cuts of one region hold each row once. The table itself when nothing is cut.
 */
export function sliceTable<T extends Table>(t: T, level: number, cut: Cut): T {
  const { start, end } = cut;
  const from = cut.from ?? start;
  const to = cut.to ?? end;
  const anchored = (s: number) => {
    const a = Math.max(s, start);
    return a >= from && a < to;
  };
  switch (t.form) {
    case 'letters':
    case 'positions': {
      const n = tableRows(t);
      const a = Math.max(0, Math.min(n, from - t.start));
      const b = Math.max(a, Math.min(n, to - t.start));
      if (a === 0 && b === n) return t;
      if (t.form === 'letters') return { ...t, start: t.start + a, letters: t.letters.slice(a, b) } as T;
      return { ...t, start: t.start + a, value: t.value.slice(a, b) } as T;
    }
    case 'bins': {
      const n = t.count.length;
      let a = n;
      let b = 0;
      for (let i = 0; i < n; i++) {
        const s = t.start + i * level;
        if (s < end && s + level > start && anchored(s)) {
          a = Math.min(a, i);
          b = i + 1;
        }
      }
      if (b <= a) a = b = 0;
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
    case 'list': {
      const rows: number[] = [];
      for (let i = 0; i < t.start.length; i++) {
        const s = t.start[i];
        if (s < end && Math.max(t.end[i], s + 1) > start && anchored(s)) rows.push(i);
      }
      const derive = !t.id && !t.row && cut.rowBase !== undefined;
      if (rows.length === t.start.length && !derive) return t;
      const row = derive ? Int32Array.from(rows, (i) => cut.rowBase! + i) : t.row && gather(t.row, rows);
      return {
        form: 'list',
        start: gather(t.start, rows),
        end: gather(t.end, rows),
        ...(t.strand && { strand: gather(t.strand, rows) }),
        ...(t.id && { id: gather(t.id, rows) }),
        ...(row && { row }),
        ...(t.fields && { fields: t.fields.map((f) => ({ name: f.name, column: gather(f.column, rows) })) }),
      } as T;
    }
    case 'density': {
      const rows: number[] = [];
      for (let i = 0; i < t.bin.length; i++) {
        const s = t.bin[i] * level;
        if (s < end && s + level > start && anchored(s)) rows.push(i);
      }
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

/**
 * Tables of one shape from consecutive parts of a region joined in order, or a message when they cannot be: forms
 * differ, letters or values leave a gap, lists carry different columns.
 */
export function joinTables(parts: readonly Table[], level: number): Table | string {
  if (parts.length === 1) return parts[0];
  const form = parts[0].form;
  if (parts.some((p) => p.form !== form)) return 'the parts are of different forms';
  const nonEmpty = parts.filter((p) => tableRows(p) > 0);
  if (nonEmpty.length <= 1) return nonEmpty[0] ?? parts[0];
  switch (form) {
    case 'letters':
    case 'positions':
    case 'bins': {
      const step = form === 'bins' ? level : 1;
      for (let i = 1; i < nonEmpty.length; i++) {
        const p = nonEmpty[i - 1] as { start: number };
        if ((nonEmpty[i] as { start: number }).start !== p.start + tableRows(nonEmpty[i - 1]) * step)
          return 'the parts leave a gap';
      }
      if (form === 'letters') {
        const ls = nonEmpty as LettersTable[];
        if (ls.every((l) => typeof l.letters === 'string'))
          return { form, start: ls[0].start, letters: ls.map((l) => l.letters).join('') };
        if (ls.some((l) => typeof l.letters === 'string' || l.alphabet !== ls[0].alphabet))
          return 'the parts have different alphabets';
        return {
          form,
          start: ls[0].start,
          alphabet: ls[0].alphabet,
          letters: concatColumns(ls.map((l) => l.letters as Uint8Array)),
        };
      }
      if (form === 'positions') {
        const ps = nonEmpty as ValuesPositions[];
        return { form, start: ps[0].start, value: concatColumns(ps.map((p) => p.value)) };
      }
      const bs = nonEmpty as ValuesBins[];
      return {
        form,
        start: bs[0].start,
        min: concatColumns(bs.map((b) => b.min)),
        max: concatColumns(bs.map((b) => b.max)),
        mean: concatColumns(bs.map((b) => b.mean)),
        count: concatColumns(bs.map((b) => b.count)),
      };
    }
    case 'list': {
      const ls = nonEmpty as RangesList[];
      const sig = (l: (typeof ls)[number]) =>
        [!!l.strand, !!l.id, !!l.row, ...(l.fields ?? []).map((f) => `${f.name}:${columnType(f.column)}`)].join('|');
      if (ls.some((l) => sig(l) !== sig(ls[0]))) return 'the parts carry different columns';
      const fields: RangesField[] | undefined = ls[0].fields?.map((f, j) => ({
        name: f.name,
        column: concatColumns(ls.map((l) => l.fields![j].column)),
      }));
      return {
        form,
        start: concatColumns(ls.map((l) => l.start)),
        end: concatColumns(ls.map((l) => l.end)),
        ...(ls[0].strand && { strand: concatColumns(ls.map((l) => l.strand!)) }),
        ...(ls[0].id && { id: concatColumns(ls.map((l) => l.id!)) }),
        ...(ls[0].row && { row: concatColumns(ls.map((l) => l.row!)) }),
        ...(fields && { fields }),
      };
    }
    case 'density': {
      const ds = nonEmpty as RangesDensity[];
      if (ds.some((d) => !!d.category !== !!ds[0].category)) return 'the parts carry different columns';
      return {
        form,
        bin: concatColumns(ds.map((d) => d.bin)),
        count: concatColumns(ds.map((d) => d.count)),
        ...(ds[0].category && { category: concatColumns(ds.map((d) => d.category!)) }),
      };
    }
  }
}
