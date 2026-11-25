/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CategoryColumn } from '@molstar/core/repr/tree';
import type {
  GridDense,
  GridSparse,
  GridTable,
  MatrixSummary,
  PairsBins,
  PairsList,
  PairsTable,
  MatrixItems,
  MatrixTable,
  RegionRequest,
  RowsKeyed,
  Table,
} from './protocol.js';
import { categoryOf, concatColumns, isInt } from './tables-columns.js';
import {
  checkForm,
  formBytes,
  formColumnCount,
  formRows,
  gatherForm,
  joinForm,
  KeyedForm,
  MatrixItemsForm,
  offsetOutside,
  PairsBinsForm,
  PairsListForm,
  type RowForm,
  rangesOutside,
  SparseForm,
  strandsOutside,
} from './tables-forms.js';

// The region shapes' tables (_spec/README.md "Data", Tables): `pairs`, `grid`, `matrix` and `rows`, their rows, bytes,
// checks, cuts and joins. The row forms are declared once (`tables-forms.ts`); a dense grid and a matrix's summary are
// by offset. Pairs share the 1D forms' names (`list`, `bins`), so a table's shape is told by its columns
// (`regionShapeOf`).

export type RegionShape = 'pairs' | 'grid' | 'matrix' | 'rows';
export type RegionTable = PairsTable | GridTable | MatrixTable | RowsKeyed;

/** A region shape's table's shape by its columns; undefined for the 1D shapes' tables. */
export function regionShapeOf(t: Table): RegionShape | undefined {
  switch (t.form) {
    case 'list':
    case 'bins':
      return 'aStart' in t || 'aBin' in t ? 'pairs' : undefined;
    case 'dense':
    case 'sparse':
      return 'grid';
    case 'items':
    case 'summary':
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
  matrix: ['items', 'summary'],
  rows: ['keyed'],
};

/** A row table's form; undefined for the forms by offset (a dense grid, a matrix's summary). */
function rowForm(t: RegionTable): RowForm | undefined {
  switch (t.form) {
    case 'list':
      return PairsListForm;
    case 'bins':
      return PairsBinsForm;
    case 'items':
      return MatrixItemsForm;
    case 'sparse':
      return SparseForm;
    case 'keyed':
      return KeyedForm;
    default:
      return undefined;
  }
}

/** A matrix summary's columns: its values a column a state. */
const SummaryCols = ['mean', 'min', 'max', 'count', 'classes', 'shares'] as const;

function summaryCols(t: MatrixSummary): (Float32Array | Uint8Array)[] {
  const out: (Float32Array | Uint8Array)[] = [];
  for (const k of SummaryCols) if (t[k]) out.push(t[k]!);
  return out;
}

function summaryColumns(t: MatrixSummary): number {
  const c = summaryCols(t)[0];
  return c && t.states > 0 ? Math.floor(c.length / t.states) : 0;
}

/** The rows (pairs, bin pairs, cells, items, bins, records) a region table holds. */
export function regionRows(t: RegionTable): number {
  if (t.form === 'dense') return t.nx * t.ny;
  if (t.form === 'summary') return summaryColumns(t);
  return formRows(t, rowForm(t)!);
}

/** The columns a region table holds besides its rows' extent: what a block's values are counted by. */
export function regionColumnCount(t: RegionTable): number {
  if (t.form === 'dense') return 1;
  if (t.form === 'summary') return t.states;
  return formColumnCount(t, rowForm(t)!);
}

/** What a region table costs the store's byte bound: its columns (a dense grid tile `nx × ny × 4`, sparse 12 a cell). */
export function regionBytes(t: RegionTable): number {
  if (t.form === 'dense') return t.value.byteLength;
  if (t.form === 'summary') return summaryCols(t).reduce((a, c) => a + c.byteLength, 0);
  return formBytes(t, rowForm(t)!);
}

/** What is wrong with a region table's columns, as a message, or `undefined`. */
export function checkRegionTable(t: RegionTable): string | undefined {
  switch (t.form) {
    case 'dense':
      if (![t.x, t.y, t.nx, t.ny].every(isInt) || t.x < 0 || t.y < 0 || t.nx < 0 || t.ny < 0)
        return '`x`, `y`, `nx`, `ny` are not integers of at least 0';
      if (!(t.value instanceof Float32Array)) return '`value` is not float32';
      return t.value.length === t.nx * t.ny ? undefined : '`value` is not `nx × ny` cells';
    case 'summary':
      return checkSummary(t);
    case 'sparse':
      if (!isInt(t.x) || !isInt(t.y)) return '`x` and `y` are not integers';
      if (typeof t.fill !== 'number') return '`fill` is not a number';
      break;
    case 'keyed':
      if (!isInt(t.start)) return '`start` is not an integer';
      if (!Array.isArray(t.fields)) return '`fields` is not a list';
      break;
  }
  return checkForm(t, rowForm(t)!);
}

function checkSummary(t: MatrixSummary): string | undefined {
  if (!isInt(t.start) || !isInt(t.states) || t.states < 0) return '`start` and `states` are not integers';
  for (const k of ['mean', 'min', 'max', 'count'] as const)
    if (t[k] !== undefined && !(t[k] instanceof Float32Array)) return 'a summary is not float32';
  for (const c of [t.classes, t.shares])
    if (c !== undefined && !(c instanceof Uint8Array)) return '`classes` and `shares` are not uint8';
  if (!!t.classes !== !!t.shares) return '`classes` without `shares`';
  const cols = summaryCols(t);
  if (!cols.length) return 'no summary';
  const n = cols[0].length;
  if (cols.some((c) => c.length !== n)) return 'columns of unequal lengths';
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
      return strandsOutside(t.aStrand, t.bStrand);
    case 'bins':
    case 'summary':
      if (L === 1) return 'bins at level 1';
      if (r.items) return 'bins asked as items';
      if (t.form === 'summary') return offsetOutside(t.start, n, L, r);
      if (n > maxRows) return `more than ${maxRows} rows`;
      for (let i = 0; i < n; i++)
        if (!binEndsIn(t, i, L, r.system, r.start, r.end, sideEnds(r.side))) return 'a bin pair outside the region';
      return strandsOutside(t.strand);
    case 'dense':
    case 'sparse': {
      if (!r.y) return 'a grid without a second axis asked';
      const xlo = Math.floor(r.start / L);
      const xhi = Math.ceil(r.end / L);
      const ylo = Math.floor(r.y.start / L);
      const yhi = Math.ceil(r.y.end / L);
      if (t.form === 'dense')
        return t.nx * t.ny > 0 && (t.x < xlo || t.x + t.nx > xhi || t.y < ylo || t.y + t.ny > yhi)
          ? 'outside the region'
          : undefined;
      for (let i = 0; i < n; i++)
        if (t.xBin[i] < xlo || t.xBin[i] >= xhi || t.yBin[i] < ylo || t.yBin[i] >= yhi)
          return 'a cell outside the region';
      return undefined;
    }
    case 'items':
      return rangesOutside(t.start, t.end, r, maxRows, 'an item');
    case 'keyed':
      return offsetOutside(t.start, n, 1, r);
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
  const rows: number[] = [];
  switch (t.form) {
    case 'list': {
      for (let i = 0; i < t.aStart.length; i++) {
        const m = pairEndsIn(t, i, cut.system, start, end, sideEnds(cut.side));
        if (m && inPart(pairAnchor(m, t.aStart[i], t.bStart[i], start))) rows.push(i);
      }
      return rowsCut(t, rows, cut) as T;
    }
    case 'bins': {
      for (let i = 0; i < t.count.length; i++) {
        const m = binEndsIn(t, i, level, cut.system, start, end, sideEnds(cut.side));
        if (m && inPart(pairAnchor(m, t.aBin[i] * level, t.bBin[i] * level, start))) rows.push(i);
      }
      return rows.length === t.count.length ? t : gatherForm(t, PairsBinsForm, rows);
    }
    case 'summary':
      return sliceSummary(t, level, from, to) as T;
    case 'dense':
    case 'sparse':
      return sliceGrid(t, level, start, from, to, cut.y) as T;
    case 'items':
      for (let i = 0; i < t.start.length; i++)
        if (meets(t.start[i], t.end[i], start, end) && inPart(Math.max(t.start[i], start))) rows.push(i);
      return rowsCut(t, rows, cut) as T;
    case 'keyed': {
      const n = t.key.length;
      const a = Math.max(0, Math.min(n, from - t.start));
      const b = Math.max(a, Math.min(n, to - t.start));
      if (a === 0 && b === n) return t;
      for (let i = a; i < b; i++) rows.push(i);
      return { ...gatherForm(t, KeyedForm, rows), start: t.start + a };
    }
  }
}

/** A list's or items' rows kept: gathered, a `row` derived from `rowBase` when it has no identity of its own. */
function rowsCut(t: PairsList | MatrixItems, rows: number[], cut: RegionCut): RegionTable {
  const f = t.form === 'list' ? PairsListForm : MatrixItemsForm;
  const own = t as { id?: unknown; row?: unknown };
  const derive = !own.id && !own.row && cut.rowBase !== undefined;
  if (rows.length === formRows(t, f) && !derive) return t;
  const out = gatherForm(t, f, rows);
  if (derive) (out as { row?: Int32Array }).row = Int32Array.from(rows, (i) => cut.rowBase! + i);
  return out;
}

function sliceSummary(t: MatrixSummary, level: number, from: number, to: number): MatrixSummary {
  const n = summaryColumns(t);
  const a = Math.max(0, Math.min(n, Math.floor((from - t.start) / level)));
  const b = Math.max(a, Math.min(n, Math.ceil((to - t.start) / level)));
  if (a === 0 && b === n) return t;
  const s = t.states;
  const out: Record<string, unknown> = { form: 'summary', start: t.start + a * level, states: s };
  for (const k of SummaryCols) if (t[k]) out[k] = t[k]!.slice(a * s, b * s);
  return out as unknown as MatrixSummary;
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
    return { ...gatherForm(t, SparseForm, rows), x: t.x, y: t.y };
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
  switch (form) {
    case 'dense': {
      const ds = nonEmpty as GridDense[];
      for (let i = 1; i < ds.length; i++)
        if (ds[i].y !== ds[0].y || ds[i].ny !== ds[0].ny || ds[i].x !== ds[i - 1].x + ds[i - 1].nx)
          return 'the parts leave a gap';
      // x-major: consecutive columns of x are consecutive runs of values
      const nx = ds.reduce((a, d) => a + d.nx, 0);
      return { form, x: ds[0].x, y: ds[0].y, nx, ny: ds[0].ny, value: concatColumns(ds.map((d) => d.value)) };
    }
    case 'summary':
      return joinSummary(nonEmpty as MatrixSummary[], level);
    case 'sparse': {
      const ss = nonEmpty as GridSparse[];
      if (ss.some((s) => !Object.is(s.fill, ss[0].fill))) return 'the parts have different fills';
      const out = joinForm(ss, SparseForm);
      return typeof out === 'string'
        ? out
        : { ...out, x: Math.min(...ss.map((s) => s.x)), y: Math.min(...ss.map((s) => s.y)) };
    }
    case 'keyed': {
      const ks = nonEmpty as RowsKeyed[];
      for (let i = 1; i < ks.length; i++)
        if (ks[i].start !== ks[i - 1].start + ks[i - 1].key.length) return 'the parts leave a gap';
      const out = joinForm(ks, KeyedForm);
      return typeof out === 'string' ? 'the parts leave a gap' : { ...out, start: ks[0].start };
    }
    default:
      return joinForm(nonEmpty, rowForm(nonEmpty[0])!);
  }
}

function joinSummary(bs: readonly MatrixSummary[], level: number): MatrixSummary | string {
  for (let i = 1; i < bs.length; i++)
    if (bs[i].states !== bs[0].states || bs[i].start !== bs[i - 1].start + summaryColumns(bs[i - 1]) * level)
      return 'the parts leave a gap';
  const shape = (b: MatrixSummary) => SummaryCols.map((k) => !!b[k]).join();
  if (bs.some((b) => shape(b) !== shape(bs[0]))) return 'the parts carry different columns';
  const out: Record<string, unknown> = { form: 'summary', start: bs[0].start, states: bs[0].states };
  for (const k of SummaryCols) if (bs[0][k]) out[k] = concatColumns(bs.map((b) => b[k] as Float32Array | Uint8Array));
  return out as unknown as MatrixSummary;
}
