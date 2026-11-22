/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, type Column, columnType, isColumn } from '@molstar/core/repr/tree';
import type {
  GridDense,
  GridSparse,
  GridTable,
  MatrixBins,
  MatrixItems,
  MatrixTable,
  PairsBins,
  PairsList,
  PairsTable,
  RangesField,
  RegionRequest,
  RowsKeyed,
  Table,
} from './protocol.js';
import { categoryOf, concatColumns, gather, isInt, isObject, sameLength } from './tables-columns.js';

// The region shapes' tables (_spec/README.md "Data", Tables): `pairs`, `grid`, `matrix` and `rows`, their rows, bytes,
// checks, cuts and joins. Their forms share names with the 1D shapes' (`list`, `bins`), so a table's shape is told by
// its columns (`regionShapeOf`).

export type RegionShape = 'pairs' | 'grid' | 'matrix' | 'rows';
export type RegionTable = PairsTable | GridTable | MatrixTable | RowsKeyed;

/** A region shape's table's shape by its columns; undefined for the 1D shapes' tables. */
export function regionShapeOf(t: Table): RegionShape | undefined {
  switch (t.form) {
    case 'list':
      return 'aStart' in t ? 'pairs' : undefined;
    case 'bins':
      return 'aBin' in t ? 'pairs' : 'states' in t ? 'matrix' : undefined;
    case 'dense':
    case 'sparse':
      return 'grid';
    case 'items':
      return 'matrix';
    case 'keyed':
      return 'rows';
    default:
      return undefined;
  }
}

/** The forms each region shape has. */
export const RegionForms: Record<RegionShape, readonly string[]> = {
  pairs: ['list', 'bins'],
  grid: ['dense', 'sparse'],
  matrix: ['items', 'bins'],
  rows: ['keyed'],
};

function matrixColumns(t: MatrixBins): number {
  const c = t.mean ?? t.count ?? t.min ?? t.max ?? t.classes ?? t.shares;
  return c && t.states > 0 ? Math.floor(c.length / t.states) : 0;
}

/** The rows (pairs, bin pairs, cells, items, bins, records) a region table holds. */
export function regionRows(t: RegionTable): number {
  switch (t.form) {
    case 'list':
      return t.aStart.length;
    case 'bins':
      return 'aBin' in t ? t.count.length : matrixColumns(t);
    case 'dense':
      return t.nx * t.ny;
    case 'sparse':
      return t.value.length;
    case 'items':
      return t.start.length;
    case 'keyed':
      return t.key.length;
  }
}

/** The columns a region table holds besides its rows' extent: what a block's values are counted by. */
export function regionColumnCount(t: RegionTable): number {
  switch (t.form) {
    case 'list':
      return 4 + [t.aSystem, t.bSystem, t.aStrand, t.bStrand, t.type, t.id, t.row].filter(Boolean).length;
    case 'bins':
      if ('aBin' in t) return 3 + [t.aSystem, t.bSystem, t.strand, t.category, t.sum].filter(Boolean).length;
      return t.states;
    case 'dense':
      return 1;
    case 'sparse':
      return 3;
    case 'items':
      return 2 + Math.max(1, t.states);
    case 'keyed':
      return 1 + t.fields.length;
  }
}

export function columnBytes(c: Column | undefined): number {
  if (c === undefined) return 0;
  if (c instanceof CategoryColumn) {
    let n = c.codes.byteLength;
    for (const s of c.dict) n += 2 * s.length;
    return n;
  }
  return c.byteLength;
}

function fieldsBytes(fields: readonly RangesField[] | undefined): number {
  let n = 0;
  for (const f of fields ?? []) n += columnBytes(f.column);
  return n;
}

/** What a region table costs the store's byte bound: its columns (a dense grid tile `nx × ny × 4`, sparse 12 a cell). */
export function regionBytes(t: RegionTable): number {
  switch (t.form) {
    case 'list':
      return (
        [t.aStart, t.aEnd, t.bStart, t.bEnd, t.aSystem, t.bSystem, t.aStrand, t.bStrand, t.type, t.id, t.row]
          .map(columnBytes)
          .reduce((a, b) => a + b, 0) + fieldsBytes(t.fields)
      );
    case 'bins':
      if ('aBin' in t)
        return [t.aBin, t.bBin, t.count, t.aSystem, t.bSystem, t.strand, t.category, t.sum]
          .map(columnBytes)
          .reduce((a, b) => a + b, 0);
      return [t.mean, t.min, t.max, t.count, t.classes, t.shares].map(columnBytes).reduce((a, b) => a + b, 0);
    case 'dense':
      return t.value.byteLength;
    case 'sparse':
      return t.xBin.byteLength + t.yBin.byteLength + t.value.byteLength;
    case 'items':
      return [t.start, t.end, t.value, t.codes, t.id, t.row].map(columnBytes).reduce((a, b) => a + b, 0);
    case 'keyed':
      return columnBytes(t.key) + fieldsBytes(t.fields);
  }
}

function checkFields(fields: unknown, n: number): string | undefined {
  if (fields === undefined) return undefined;
  if (!Array.isArray(fields)) return '`fields` is not a list';
  for (const f of fields) {
    if (!isObject(f) || typeof f.name !== 'string' || !isColumn(f.column as Column))
      return 'a field is not a named column';
    if ((f.column as Column).length !== n) return 'columns of unequal lengths';
  }
  return undefined;
}

const isCat = (c: unknown) => c === undefined || c instanceof CategoryColumn;
const isI32 = (c: unknown) => c instanceof Int32Array;
const isF32 = (c: unknown) => c instanceof Float32Array;
const isI8 = (c: unknown) => c === undefined || c instanceof Int8Array;
const isU8 = (c: unknown) => c === undefined || c instanceof Uint8Array;

/** What is wrong with a region table's columns, as a message, or `undefined`. */
export function checkRegionTable(t: RegionTable): string | undefined {
  switch (t.form) {
    case 'list': {
      if (![t.aStart, t.aEnd, t.bStart, t.bEnd].every(isI32)) return 'the ends are not int32';
      if (!isI8(t.aStrand) || !isI8(t.bStrand)) return 'a strand is not int8';
      if (![t.aSystem, t.bSystem, t.type, t.id].every(isCat)) return 'a system, `type` or `id` is not a category';
      if (t.row !== undefined && !isI32(t.row)) return '`row` is not int32';
      const n = t.aStart.length;
      if (!sameLength(n, t.aEnd, t.bStart, t.bEnd, t.aStrand, t.bStrand, t.aSystem, t.bSystem, t.type, t.id, t.row))
        return 'columns of unequal lengths';
      if (!!t.aSystem !== !!t.bSystem) return 'one end’s system without the other’s';
      return checkFields(t.fields, n);
    }
    case 'bins':
      if ('aBin' in t) {
        if (!isI32(t.aBin) || !isI32(t.bBin) || !isF32(t.count))
          return '`aBin`, `bBin`, `count` are not int32 and float32';
        if (!isI8(t.strand) || !isCat(t.category) || !isCat(t.aSystem) || !isCat(t.bSystem))
          return 'a strand, system or category of another type';
        if (t.sum !== undefined && !isF32(t.sum)) return '`sum` is not float32';
        if (!!t.aSystem !== !!t.bSystem) return 'one end’s system without the other’s';
        return sameLength(t.count.length, t.aBin, t.bBin, t.strand, t.category, t.sum, t.aSystem, t.bSystem)
          ? undefined
          : 'columns of unequal lengths';
      }
      return checkMatrixBins(t);
    case 'dense':
      if (![t.x, t.y, t.nx, t.ny].every(isInt) || t.x < 0 || t.y < 0 || t.nx < 0 || t.ny < 0)
        return '`x`, `y`, `nx`, `ny` are not integers of at least 0';
      if (!isF32(t.value)) return '`value` is not float32';
      return t.value.length === t.nx * t.ny ? undefined : '`value` is not `nx × ny` cells';
    case 'sparse':
      if (!isInt(t.x) || !isInt(t.y)) return '`x` and `y` are not integers';
      if (!isI32(t.xBin) || !isI32(t.yBin) || !isF32(t.value)) return 'cells are not int32 bins and float32 values';
      if (typeof t.fill !== 'number') return '`fill` is not a number';
      return sameLength(t.value.length, t.xBin, t.yBin) ? undefined : 'columns of unequal lengths';
    case 'items': {
      if (!isI32(t.start) || !isI32(t.end)) return '`start` and `end` are not int32';
      if (!isInt(t.states) || t.states < 0) return '`states` is not an integer of at least 0';
      if ((t.value === undefined) === (t.codes === undefined)) return 'not one of `value` and `codes`';
      if (t.value !== undefined && !isF32(t.value)) return '`value` is not float32';
      if (!isU8(t.codes)) return '`codes` are not uint8';
      const n = t.start.length;
      if ((t.value ?? t.codes)!.length !== n * t.states) return 'cells are not items × states';
      if (!isCat(t.id) || (t.row !== undefined && !isI32(t.row))) return '`id` or `row` of another type';
      return sameLength(n, t.end, t.id, t.row) ? undefined : 'columns of unequal lengths';
    }
    case 'keyed':
      if (!isInt(t.start)) return '`start` is not an integer';
      if (!(t.key instanceof CategoryColumn)) return '`key` is not a category';
      return checkFields(t.fields ?? [], t.key.length);
  }
}

function checkMatrixBins(t: MatrixBins): string | undefined {
  if (!isInt(t.start) || !isInt(t.states) || t.states < 0) return '`start` and `states` are not integers';
  for (const c of [t.mean, t.min, t.max, t.count]) if (c !== undefined && !isF32(c)) return 'a summary is not float32';
  if (!isU8(t.classes) || !isU8(t.shares)) return '`classes` and `shares` are not uint8';
  if (!!t.classes !== !!t.shares) return '`classes` without `shares`';
  const cols = [t.mean, t.min, t.max, t.count, t.classes, t.shares].filter((c) => c !== undefined);
  if (!cols.length) return 'no summary';
  const n = cols[0]!.length;
  if (cols.some((c) => c!.length !== n)) return 'columns of unequal lengths';
  return t.states === 0 || n % t.states === 0 ? undefined : 'bins are not whole columns';
}

/** Whether an end `[s, e)` meets `[start, end)` (a point end, `s = e`, at its position). */
function meets(s: number, e: number, start: number, end: number): boolean {
  return s < end && Math.max(e, s + 1) > start;
}

/** Whether end `a` (else `b`) of pair `i` is on `system` (no systems column: both are). */
function onSystem(sys: CategoryColumn | undefined, i: number, system: string | undefined): boolean {
  return !sys || system === undefined || categoryOf(sys, i) === system;
}

/** The ends a request reads: a sided one its side's alone (1 `a`, 2 `b`), an unsided one either (3). */
export function sideEnds(side: 'a' | 'b' | undefined): number {
  return side === 'a' ? 1 : side === 'b' ? 2 : 3;
}

/** Each end of pair `i` meeting the region on its system: 1 for `a`, 2 for `b`, 3 both, 0 none (of `ends`). */
export function pairEndsIn(
  t: PairsList,
  i: number,
  system: string | undefined,
  start: number,
  end: number,
  ends = 3,
): number {
  let m = 0;
  if (ends & 1 && onSystem(t.aSystem, i, system) && meets(t.aStart[i], t.aEnd[i], start, end)) m |= 1;
  if (ends & 2 && onSystem(t.bSystem, i, system) && meets(t.bStart[i], t.bEnd[i], start, end)) m |= 2;
  return m;
}

/** Each end of bin pair `i` whose bin meets the region at `level`: as `pairEndsIn`. */
export function binEndsIn(
  t: PairsBins,
  i: number,
  level: number,
  system: string | undefined,
  start: number,
  end: number,
  ends = 3,
) {
  let m = 0;
  if (ends & 1 && onSystem(t.aSystem, i, system) && meets(t.aBin[i] * level, (t.aBin[i] + 1) * level, start, end))
    m |= 1;
  if (ends & 2 && onSystem(t.bSystem, i, system) && meets(t.bBin[i] * level, (t.bBin[i] + 1) * level, start, end))
    m |= 2;
  return m;
}

/** The region of bins `[lo, hi)` a request covers at its level. */
function binSpan(start: number, end: number, level: number): [number, number] {
  return [Math.floor(start / level), Math.ceil(end / level)];
}

/** What a region table holds outside its request, as a message, or `undefined` (`checkAnswer`'s extent). */
export function checkRegionExtent(t: RegionTable, r: RegionRequest, maxRows: number): string | undefined {
  const L = r.level;
  const n = regionRows(t);
  switch (t.form) {
    case 'list':
      if (n > maxRows) return `more than ${maxRows} rows`;
      for (let i = 0; i < n; i++) {
        if (t.aStart[i] > t.aEnd[i] || t.bStart[i] > t.bEnd[i]) return 'an end ends before it starts';
        if (!pairEndsIn(t, i, r.system, r.start, r.end, sideEnds(r.side))) return 'a pair with no end in the region';
      }
      for (const s of [t.aStrand, t.bStrand])
        if (s) for (let i = 0; i < n; i++) if (s[i] < -1 || s[i] > 1) return 'a strand not −1, 0, 1';
      return undefined;
    case 'bins':
      if (L === 1) return 'bins at level 1';
      if (r.items) return 'bins asked as items';
      if ('aBin' in t) {
        if (n > maxRows) return `more than ${maxRows} rows`;
        for (let i = 0; i < n; i++)
          if (!binEndsIn(t, i, L, r.system, r.start, r.end, sideEnds(r.side))) return 'a bin pair outside the region';
        if (t.strand)
          for (let i = 0; i < n; i++) if (t.strand[i] < -1 || t.strand[i] > 1) return 'a strand not −1, 0, 1';
        return undefined;
      }
      if (t.start % L !== 0) return 'bins not aligned to the level';
      if (n > 0 && (t.start < Math.floor(r.start / L) * L || t.start + n * L > Math.ceil(r.end / L) * L))
        return 'outside the region';
      return undefined;
    case 'dense':
    case 'sparse': {
      if (!r.y) return 'a grid without a second axis asked';
      const [xlo, xhi] = binSpan(r.start, r.end, L);
      const [ylo, yhi] = binSpan(r.y.start, r.y.end, L);
      if (t.form === 'dense') {
        if (t.nx * t.ny > 0 && (t.x < xlo || t.x + t.nx > xhi || t.y < ylo || t.y + t.ny > yhi))
          return 'outside the region';
        return undefined;
      }
      for (let i = 0; i < n; i++)
        if (t.xBin[i] < xlo || t.xBin[i] >= xhi || t.yBin[i] < ylo || t.yBin[i] >= yhi)
          return 'a cell outside the region';
      return undefined;
    }
    case 'items':
      if (n > maxRows) return `more than ${maxRows} rows`;
      for (let i = 0; i < n; i++) {
        if (t.start[i] > t.end[i]) return 'an item ends before it starts';
        if (!meets(t.start[i], t.end[i], r.start, r.end)) return 'an item outside the region';
      }
      return undefined;
    case 'keyed':
      if (n > 0 && (t.start < r.start || t.start + n > r.end)) return 'outside the region';
      return undefined;
  }
}

/** Where a region cut keeps a row (`Cut`, with the system its ends are read on and the grid's second axis). */
export interface RegionCut {
  readonly start: number;
  readonly end: number;
  readonly from?: number;
  readonly to?: number;
  readonly rowBase?: number;
  /** The system the region is on: a pair's end on another system never anchors it. */
  readonly system?: string;
  /** A sided request's side: only that end anchors a pair. */
  readonly side?: 'a' | 'b';
  /** A grid's second axis, cut to `[start, end)` (no `from`: grid parts join along x only). */
  readonly y?: { readonly start: number; readonly end: number };
}

/**
 * The anchor of a pair in a region: its ends meeting the region, the first position of the earlier one in it. A pair
 * is kept by the one cut whose `[from, to)` holds it, so joined cuts hold each pair once though both its ends meet.
 */
function pairAnchor(m: number, as: number, bs: number, start: number): number {
  const a = m & 1 ? Math.max(as, start) : Number.POSITIVE_INFINITY;
  const b = m & 2 ? Math.max(bs, start) : Number.POSITIVE_INFINITY;
  return Math.min(a, b);
}

/** A region table cut to a region at `level` (as `sliceTable`); the table itself when nothing is cut. */
export function sliceRegion<T extends RegionTable>(t: T, level: number, cut: RegionCut): T {
  const { start, end } = cut;
  const from = cut.from ?? start;
  const to = cut.to ?? end;
  const inPart = (a: number) => a >= from && a < to;
  switch (t.form) {
    case 'list': {
      const rows: number[] = [];
      for (let i = 0; i < t.aStart.length; i++) {
        const m = pairEndsIn(t, i, cut.system, start, end, sideEnds(cut.side));
        if (m && inPart(pairAnchor(m, t.aStart[i], t.bStart[i], start))) rows.push(i);
      }
      const derive = !t.id && !t.row && cut.rowBase !== undefined;
      if (rows.length === t.aStart.length && !derive) return t;
      const row = derive ? Int32Array.from(rows, (i) => cut.rowBase! + i) : t.row && gather(t.row, rows);
      const g = <C extends Column>(c: C | undefined) => c && gather(c, rows);
      return {
        form: 'list',
        ...(t.aSystem && { aSystem: g(t.aSystem), bSystem: g(t.bSystem) }),
        aStart: g(t.aStart)!,
        aEnd: g(t.aEnd)!,
        bStart: g(t.bStart)!,
        bEnd: g(t.bEnd)!,
        ...(t.aStrand && { aStrand: g(t.aStrand) }),
        ...(t.bStrand && { bStrand: g(t.bStrand) }),
        ...(t.type && { type: g(t.type) }),
        ...(t.id && { id: g(t.id) }),
        ...(row && { row }),
        ...(t.fields && { fields: t.fields.map((f) => ({ name: f.name, column: gather(f.column, rows) })) }),
      } as T;
    }
    case 'bins': {
      if (!('aBin' in t)) return sliceMatrixBins(t, level, from, to) as T;
      const rows: number[] = [];
      for (let i = 0; i < t.count.length; i++) {
        const m = binEndsIn(t, i, level, cut.system, start, end, sideEnds(cut.side));
        if (m && inPart(pairAnchor(m, t.aBin[i] * level, t.bBin[i] * level, start))) rows.push(i);
      }
      if (rows.length === t.count.length) return t;
      const g = <C extends Column>(c: C | undefined) => c && gather(c, rows);
      return {
        form: 'bins',
        ...(t.aSystem && { aSystem: g(t.aSystem), bSystem: g(t.bSystem) }),
        aBin: g(t.aBin)!,
        bBin: g(t.bBin)!,
        ...(t.strand && { strand: g(t.strand) }),
        ...(t.category && { category: g(t.category) }),
        count: g(t.count)!,
        ...(t.sum && { sum: g(t.sum) }),
      } as T;
    }
    case 'dense':
    case 'sparse':
      return sliceGrid(t, level, start, from, to, cut.y) as T;
    case 'items': {
      const rows: number[] = [];
      for (let i = 0; i < t.start.length; i++)
        if (meets(t.start[i], t.end[i], start, end) && inPart(Math.max(t.start[i], start))) rows.push(i);
      const derive = !t.id && !t.row && cut.rowBase !== undefined;
      if (rows.length === t.start.length && !derive) return t;
      const row = derive ? Int32Array.from(rows, (i) => cut.rowBase! + i) : t.row && gather(t.row, rows);
      return {
        form: 'items',
        start: gather(t.start, rows),
        end: gather(t.end, rows),
        states: t.states,
        ...(t.value && { value: gatherBlocks(t.value, rows, t.states) }),
        ...(t.codes && { codes: gatherBlocks(t.codes, rows, t.states) }),
        ...(t.id && { id: gather(t.id, rows) }),
        ...(row && { row }),
      } as T;
    }
    case 'keyed': {
      const n = t.key.length;
      const a = Math.max(0, Math.min(n, from - t.start));
      const b = Math.max(a, Math.min(n, to - t.start));
      if (a === 0 && b === n) return t;
      const rows = Array.from({ length: b - a }, (_, i) => a + i);
      return {
        form: 'keyed',
        start: t.start + a,
        key: gather(t.key, rows),
        fields: t.fields.map((f) => ({ name: f.name, column: gather(f.column, rows) })),
      } as RowsKeyed as T;
    }
  }
}

/** Blocks of `size` values at the given block indices, in order (a matrix's states per item). */
function gatherBlocks<A extends Float32Array | Uint8Array>(a: A, rows: readonly number[], size: number): A {
  const out = new (a.constructor as new (n: number) => A)(rows.length * size);
  for (let i = 0; i < rows.length; i++) out.set(a.subarray(rows[i] * size, (rows[i] + 1) * size), i * size);
  return out;
}

function sliceMatrixBins(t: MatrixBins, level: number, from: number, to: number): MatrixBins {
  const n = matrixColumns(t);
  const a = Math.max(0, Math.min(n, Math.floor((from - t.start) / level)));
  const b = Math.max(a, Math.min(n, Math.ceil((to - t.start) / level)));
  if (a === 0 && b === n) return t;
  const s = t.states;
  const cut = <A extends Float32Array | Uint8Array>(c: A | undefined) => c && (c.slice(a * s, b * s) as A);
  return {
    form: 'bins',
    start: t.start + a * level,
    states: s,
    ...(t.mean && { mean: cut(t.mean) }),
    ...(t.min && { min: cut(t.min) }),
    ...(t.max && { max: cut(t.max) }),
    ...(t.count && { count: cut(t.count) }),
    ...(t.classes && { classes: cut(t.classes), shares: cut(t.shares) }),
  };
}

/** A grid tile cut to the bins anchored in `[from, to)` on x and meeting `y` on y. */
function sliceGrid(
  t: GridTable,
  level: number,
  start: number,
  from: number,
  to: number,
  y: { readonly start: number; readonly end: number } | undefined,
): GridTable {
  // a bin is anchored at its first position in the region: the region's first bin at `start`, the others at their own
  const xlo = from === start ? Math.floor(start / level) : Math.ceil(from / level);
  const xhi = Math.ceil(to / level);
  const [ylo, yhi] = y ? [Math.floor(y.start / level), Math.ceil(y.end / level)] : [-Infinity, Infinity];
  if (t.form === 'sparse') {
    const rows: number[] = [];
    for (let i = 0; i < t.value.length; i++)
      if (t.xBin[i] >= xlo && t.xBin[i] < xhi && t.yBin[i] >= ylo && t.yBin[i] < yhi) rows.push(i);
    if (rows.length === t.value.length) return t;
    return { ...t, xBin: gather(t.xBin, rows), yBin: gather(t.yBin, rows), value: gather(t.value, rows) };
  }
  const x0 = Math.max(t.x, Math.min(t.x + t.nx, xlo));
  const x1 = Math.max(x0, Math.min(t.x + t.nx, xhi));
  const y0 = Math.max(t.y, Math.min(t.y + t.ny, ylo));
  const y1 = Math.max(y0, Math.min(t.y + t.ny, yhi));
  if (x0 === t.x && x1 === t.x + t.nx && y0 === t.y && y1 === t.y + t.ny) return t;
  const nx = x1 - x0;
  const ny = y1 - y0;
  const value = new Float32Array(nx * ny);
  for (let i = 0; i < nx; i++) {
    const o = (x0 - t.x + i) * t.ny + (y0 - t.y);
    value.set(t.value.subarray(o, o + ny), i * ny);
  }
  return { form: 'dense', x: x0, y: y0, nx, ny, value };
}

function sig(fields: readonly RangesField[] | undefined): string {
  return (fields ?? []).map((f) => `${f.name}:${columnType(f.column)}`).join('|');
}

/**
 * Region tables of one shape from consecutive parts of a region joined in order (a grid's along x), or a message when
 * they cannot be: forms differ, columns differ, a gap.
 */
export function joinRegion(parts: readonly RegionTable[], level: number): RegionTable | string {
  if (parts.length === 1) return parts[0];
  const form = parts[0].form;
  if (parts.some((p) => p.form !== form || regionShapeOf(p) !== regionShapeOf(parts[0])))
    return 'the parts are of different forms';
  const nonEmpty = parts.filter((p) => regionRows(p) > 0);
  if (nonEmpty.length <= 1) return nonEmpty[0] ?? parts[0];
  const cat = <C extends Column>(cols: (C | undefined)[]) => (cols[0] ? concatColumns(cols as C[]) : undefined);
  switch (form) {
    case 'list': {
      const ls = nonEmpty as PairsList[];
      const shape = (l: PairsList) =>
        [l.aSystem, l.aStrand, l.bStrand, l.type, l.id, l.row].map((c) => !!c).join() + sig(l.fields);
      if (ls.some((l) => shape(l) !== shape(ls[0]))) return 'the parts carry different columns';
      const pick = <K extends keyof PairsList>(k: K) => cat(ls.map((l) => l[k] as Column | undefined));
      const out: Record<string, unknown> = { form };
      for (const k of [
        'aSystem',
        'bSystem',
        'aStart',
        'aEnd',
        'bStart',
        'bEnd',
        'aStrand',
        'bStrand',
        'type',
        'id',
        'row',
      ] as const) {
        const c = pick(k);
        if (c) out[k] = c;
      }
      if (ls[0].fields)
        out.fields = ls[0].fields.map((f, j) => ({
          name: f.name,
          column: concatColumns(ls.map((l) => l.fields![j].column)),
        }));
      return out as unknown as PairsList;
    }
    case 'bins': {
      if (!('aBin' in nonEmpty[0])) return joinMatrixBins(nonEmpty as MatrixBins[], level);
      const bs = nonEmpty as PairsBins[];
      const shape = (b: PairsBins) => [b.aSystem, b.strand, b.category, b.sum].map((c) => !!c).join();
      if (bs.some((b) => shape(b) !== shape(bs[0]))) return 'the parts carry different columns';
      const out: Record<string, unknown> = { form };
      for (const k of ['aSystem', 'bSystem', 'aBin', 'bBin', 'strand', 'category', 'count', 'sum'] as const) {
        const c = cat(bs.map((b) => b[k] as Column | undefined));
        if (c) out[k] = c;
      }
      return out as unknown as PairsBins;
    }
    case 'dense': {
      const ds = nonEmpty as GridDense[];
      for (let i = 1; i < ds.length; i++)
        if (ds[i].y !== ds[0].y || ds[i].ny !== ds[0].ny || ds[i].x !== ds[i - 1].x + ds[i - 1].nx)
          return 'the parts leave a gap';
      // x-major: consecutive columns of x are consecutive runs of values
      const nx = ds.reduce((a, d) => a + d.nx, 0);
      return { form, x: ds[0].x, y: ds[0].y, nx, ny: ds[0].ny, value: concatColumns(ds.map((d) => d.value)) };
    }
    case 'sparse': {
      const ss = nonEmpty as GridSparse[];
      if (ss.some((s) => !Object.is(s.fill, ss[0].fill))) return 'the parts have different fills';
      return {
        form,
        x: Math.min(...ss.map((s) => s.x)),
        y: Math.min(...ss.map((s) => s.y)),
        xBin: concatColumns(ss.map((s) => s.xBin)),
        yBin: concatColumns(ss.map((s) => s.yBin)),
        value: concatColumns(ss.map((s) => s.value)),
        fill: ss[0].fill,
      };
    }
    case 'items': {
      const is = nonEmpty as MatrixItems[];
      const shape = (m: MatrixItems) => [m.states, !!m.value, !!m.codes, !!m.id, !!m.row].join();
      if (is.some((m) => shape(m) !== shape(is[0]))) return 'the parts carry different columns';
      const out: Record<string, unknown> = { form, states: is[0].states };
      for (const k of ['start', 'end', 'value', 'codes', 'id', 'row'] as const) {
        const c = cat(is.map((m) => m[k] as Column | undefined));
        if (c) out[k] = c;
      }
      return out as unknown as MatrixItems;
    }
    case 'keyed': {
      const ks = nonEmpty as RowsKeyed[];
      for (let i = 1; i < ks.length; i++)
        if (ks[i].start !== ks[i - 1].start + ks[i - 1].key.length || sig(ks[i].fields) !== sig(ks[0].fields))
          return 'the parts leave a gap';
      return {
        form,
        start: ks[0].start,
        key: concatColumns(ks.map((k) => k.key)),
        fields: ks[0].fields.map((f, j) => ({
          name: f.name,
          column: concatColumns(ks.map((k) => k.fields[j].column)),
        })),
      };
    }
  }
}

function joinMatrixBins(bs: readonly MatrixBins[], level: number): MatrixBins | string {
  for (let i = 1; i < bs.length; i++)
    if (bs[i].states !== bs[0].states || bs[i].start !== bs[i - 1].start + matrixColumns(bs[i - 1]) * level)
      return 'the parts leave a gap';
  const shape = (b: MatrixBins) => [b.mean, b.min, b.max, b.count, b.classes].map((c) => !!c).join();
  if (bs.some((b) => shape(b) !== shape(bs[0]))) return 'the parts carry different columns';
  const out: Record<string, unknown> = { form: 'bins', start: bs[0].start, states: bs[0].states };
  for (const k of ['mean', 'min', 'max', 'count', 'classes', 'shares'] as const)
    if (bs[0][k]) out[k] = concatColumns(bs.map((b) => b[k] as Float32Array | Uint8Array));
  return out as unknown as MatrixBins;
}
