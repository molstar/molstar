/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Column, type ColumnType, categoryColumn } from '@molstar/core/repr/tree';
import type { FieldInfo, RangesField, RangesList, ValuesBins, ValuesPositions } from '../data/protocol.js';
import type { BedCauses } from './bed-issues.js';

// The BED family's one line parser (_spec/README.md "Data", BED family): a line's cells into a system's items (BED's
// pieces), a bedGraph's intervals or a chrom.sizes' lengths, the same whether the lines come from a file read whole or
// from a tabix region, so a region read equals the whole read clipped. Positions as the format writes them: base 0,
// ends exclusive (EC-195), never a param.

/** The largest position a column holds (int32). */
export const MaxPosition = 0x7fffffff;

/** A line's text and its place for an issue (`line 12`). */
export interface BedLine {
  readonly text: string;
  readonly place: string;
}

/** BED's columns by place (UCSC): the names a column past the third has when no mapping names it. */
export const BedColumns = [
  'chrom',
  'start',
  'end',
  'name',
  'score',
  'strand',
  'thick_start',
  'thick_end',
  'item_rgb',
  'block_count',
  'block_sizes',
  'block_starts',
] as const;

const BedTypes: Readonly<Record<string, ColumnType>> = {
  score: 'float32',
  thick_start: 'int32',
  thick_end: 'int32',
  block_count: 'int32',
};

/** A field's declared type (the source's `types`); `string` is held as a category. */
export type DeclaredType = 'int32' | 'float32' | 'float64' | 'category' | 'string';

/** Whether a line is read past: `#` comments, `track` and `browser` lines. */
export function isMetaLine(text: string): boolean {
  return text.startsWith('#') || /^(track|browser)(\s|$)/.test(text);
}

/** A line's cells: tab-separated, or split on whitespace when it has no tab (said once by the caller). */
export function cellsOf(text: string): { cells: string[]; spaces: boolean } {
  if (text.includes('\t')) return { cells: text.split('\t'), spaces: false };
  const t = text.trim();
  return { cells: t ? t.split(/\s+/) : [], spaces: true };
}

/** A whole number from 0 to int32's largest, else undefined. */
export function positionOf(cell: string | undefined): number | undefined {
  if (cell === undefined || cell.length === 0 || cell.length > 10 || !/^\d+$/.test(cell)) return undefined;
  const v = Number(cell);
  return v <= MaxPosition ? v : undefined;
}

/** `+` 1, `-` −1, `.` (or empty) 0; else undefined. */
export function strandOf(cell: string | undefined): number | undefined {
  if (cell === undefined || cell === '.' || cell === '' || cell === '?') return 0;
  if (cell === '+') return 1;
  if (cell === '-') return -1;
  return undefined;
}

/** A field a schema reads: its name, its 0-based column and its type. */
export interface FieldSpec {
  readonly name: string;
  readonly column: number;
  readonly type: ColumnType;
}

/** What a format's lines mean: the columns its first line has, the fields read, the name's and strand's columns. */
export interface BedSchema {
  readonly columns: number;
  readonly fields: readonly FieldSpec[];
  /** The column of the item's name (its identity): a mapping's `name`, else BED's fourth column when there is one. */
  readonly name?: number;
  readonly strand?: number;
  /** BED12: blocks as pieces (columns 9 to 11, none of them mapped). */
  readonly blocks: boolean;
}

const typeOf = (declared: DeclaredType | undefined, fallback: ColumnType): ColumnType =>
  declared === undefined ? fallback : declared === 'string' ? 'category' : declared;

/** A mapping's spelling as a column: an all-digit segment is its 0-based number (the records' rule), else a name. */
function columnOf(spelling: string, names: readonly string[]): number | undefined {
  if (/^\d+$/.test(spelling)) return Number(spelling);
  const i = names.indexOf(spelling);
  return i < 0 ? undefined : i;
}

/**
 * A BED file's schema from its first data line's column count (BED `N`) and the source's `fields` (a field by its
 * column's number or BED name; `name` and `strand` the roles) and `types`. A column a mapping reads is not read under
 * BED's name (EC-190); a mapping naming no column, or a position column, is the error returned.
 */
export function bedSchema(
  columns: number,
  mapping: Readonly<Record<string, string>> = {},
  types: Readonly<Record<string, DeclaredType>> = {},
): BedSchema | { readonly error: string } {
  const mapped = new Map<number, string>();
  let name: number | undefined;
  let strand: number | undefined;
  const fields: FieldSpec[] = [];
  for (const [field, spelling] of Object.entries(mapping)) {
    const c = columnOf(spelling, BedColumns);
    if (c === undefined) return { error: `The field “${field}” reads “${spelling}”, which is no column of BED` };
    if (c < 3) return { error: `The field “${field}” reads column ${c}, a position` };
    if (c >= columns) return { error: `The field “${field}” reads column ${c}; the first line has ${columns}` };
    mapped.set(c, field);
    if (field === 'strand') strand = c;
    else {
      // a mapping's `name` is the identity's column and the label
      if (field === 'name') name = c;
      const label = field === 'name' ? 'label' : field;
      fields.push({ name: label, column: c, type: typeOf(types[label], 'category') });
    }
  }
  // a `bed N+M` file: BED's names stop at the first column a mapping reads; the unmapped ones past it are `col<k>`
  const plus = mapped.size ? Math.min(...mapped.keys()) : columns;
  const blocks = columns >= 12 && plus > 11;
  // a field a mapping names under BED's name (`label`, `score`): the mapping's column holds it, BED's is not read
  const own = new Set([...fields.map((f) => f.name), ...(strand !== undefined ? ['strand'] : [])]);
  for (let c = 3; c < columns; c++) {
    if (mapped.has(c) || (blocks && c >= 9 && c <= 11)) continue;
    // BED's name is the label (the items' label field), its identity too (below)
    const bed = c >= 12 || c >= plus ? `col${c}` : c === 3 ? 'label' : BedColumns[c];
    if (own.has(bed)) continue;
    if (bed === 'strand') {
      strand = c;
      continue;
    }
    if (c === 3) name = c;
    fields.push({ name: bed, column: c, type: typeOf(types[bed], BedTypes[bed] ?? 'category') });
  }
  // identity is BED's name column (EC-196) when a mapping reads it under another field (`label`) too
  name ??= columns > 3 ? 3 : undefined;
  return { columns, fields, ...(name !== undefined && { name }), ...(strand !== undefined && { strand }), blocks };
}

/** The coverage's fields of a schema's tables. */
export function fieldInfos(schema: BedSchema): FieldInfo[] {
  const out: FieldInfo[] = schema.fields.map((f) => ({ name: f.name, type: f.type }));
  if (schema.strand !== undefined) out.push({ name: 'strand', type: 'int8' });
  return out;
}

// ---- BED: a system's items and their pieces ----

/** A system's items as lines are read, then finished into columns sorted by piece start. */
export class ItemsBuilder {
  readonly pStart: number[] = [];
  readonly pEnd: number[] = [];
  readonly pItem: number[] = [];
  readonly iStart: number[] = [];
  readonly iEnd: number[] = [];
  readonly strand: number[] = [];
  readonly key: string[] = [];
  readonly values: (number | string | null)[][];
  private readonly seen = new Map<string, number>();

  constructor(readonly schema: BedSchema) {
    this.values = schema.fields.map(() => []);
  }

  get pieces() {
    return this.pStart.length;
  }

  /**
   * Reads one data line's cells (its chrom already told); false when it is left out (said). Identity (EC-196):
   * `name@start-end`, `#2`, `#3` on a later identical key in read order (identical items are read together).
   */
  add(
    cells: readonly string[],
    place: string,
    causes: BedCauses,
    room = Number.POSITIVE_INFINITY,
  ): boolean | 'too-large' {
    const s = this.schema;
    const start = positionOf(cells[1]);
    const end = positionOf(cells[2]);
    if (start === undefined || end === undefined) {
      causes.add('seq.source.bed.line', 'A position is not a whole number from 0 to 2^31 − 1', place);
      return false;
    }
    if (end < start) {
      causes.add('seq.source.bed.line', 'An end before its start', place);
      return false;
    }
    let blockStarts: number[] | undefined;
    let blockSizes: number[] | undefined;
    if (s.blocks) {
      const got = blocksOf(cells, end - start, room);
      if (got === TooMany) return 'too-large';
      if (typeof got === 'string') {
        causes.add('seq.source.bed.line', got, place);
        return false;
      }
      [blockStarts, blockSizes] = got;
    }
    if (cells.length > s.columns)
      causes.add('seq.source.bed.cells', 'More columns than the first line: not read', place);
    let strand = 0;
    if (s.strand !== undefined) {
      const v = strandOf(cells[s.strand]);
      if (v === undefined) causes.add('seq.source.bed.cells', 'A strand not +, - or .: read as .', place);
      strand = v ?? 0;
    }
    const item = this.iStart.length;
    for (let f = 0; f < s.fields.length; f++) {
      const spec = s.fields[f];
      this.values[f].push(cellValue(cells[spec.column], spec.type, place, causes));
    }
    const n = s.name !== undefined ? cells[s.name] : undefined;
    const base = `${n === undefined || n === '.' ? '' : n}@${start}-${end}`;
    const k = (this.seen.get(base) ?? 0) + 1;
    this.seen.set(base, k);
    this.key.push(k === 1 ? base : `${base}#${k}`);
    this.iStart.push(start);
    this.iEnd.push(end);
    this.strand.push(strand);
    if (blockStarts && blockSizes) {
      for (let b = 0; b < blockStarts.length; b++) {
        this.pStart.push(start + blockStarts[b]);
        this.pEnd.push(start + blockStarts[b] + blockSizes[b]);
        this.pItem.push(item);
      }
    } else {
      this.pStart.push(start);
      this.pEnd.push(end);
      this.pItem.push(item);
    }
    return true;
  }

  finish(): Items {
    const n = this.pStart.length;
    const order = sortedOrder(this.pStart);
    const pStart = new Int32Array(n);
    const pEnd = new Int32Array(n);
    const pItem = new Int32Array(n);
    let maxSpan = 1;
    let maxEnd = 0;
    for (let i = 0; i < n; i++) {
      const k = order ? order[i] : i;
      pStart[i] = this.pStart[k];
      pEnd[i] = this.pEnd[k];
      pItem[i] = this.pItem[k];
      maxSpan = Math.max(maxSpan, pEnd[i] - pStart[i]);
      maxEnd = Math.max(maxEnd, pEnd[i]);
    }
    const s = this.schema;
    return {
      schema: s,
      pStart,
      pEnd,
      pItem,
      maxSpan,
      maxEnd,
      iStart: Int32Array.from(this.iStart),
      iEnd: Int32Array.from(this.iEnd),
      strand: s.strand !== undefined ? Int8Array.from(this.strand) : undefined,
      key: this.key,
      values: this.values,
    };
  }
}

/** A system's items finished: pieces sorted by start (file order among equal starts), each with its item. */
export interface Items {
  readonly schema: BedSchema;
  readonly pStart: Int32Array;
  readonly pEnd: Int32Array;
  readonly pItem: Int32Array;
  /** The longest piece (at least 1): a piece meeting `start` starts at `start − maxSpan` or later. */
  readonly maxSpan: number;
  readonly maxEnd: number;
  readonly iStart: Int32Array;
  readonly iEnd: Int32Array;
  readonly strand?: Int8Array;
  readonly key: readonly string[];
  readonly values: readonly (readonly (number | string | null)[])[];
}

/** The order sorting `starts` (stable), or undefined when they are sorted already. */
function sortedOrder(starts: readonly number[]): Int32Array | undefined {
  let sorted = true;
  for (let i = 1; i < starts.length && sorted; i++) if (starts[i] < starts[i - 1]) sorted = false;
  if (sorted) return undefined;
  const order = Int32Array.from({ length: starts.length }, (_, i) => i);
  order.sort((a, b) => starts[a] - starts[b] || a - b);
  return order;
}

/** BED12's blocks: counts equal to `block_count`, each inside the item, sorted and apart; else why not. */
/** A line's blocks past the rows left: refused before its lists are split. */
const TooMany = 'too-many';

/**
 * BED12's blocks by UCSC's rule: counts equal to `block_count`, each at least a position long, the first from the item's
 * start, the last to its end, sorted and apart; else why not. More blocks than `room` (the rows left) is `TooMany`, told
 * from the count before the lists are split.
 */
function blocksOf(cells: readonly string[], length: number, room: number): [number[], number[]] | string {
  const count = positionOf(cells[9]);
  if (count !== undefined && count > room) return TooMany;
  const list = (cell: string | undefined) => {
    const parts = (cell ?? '').split(',');
    if (parts[parts.length - 1] === '') parts.pop();
    return parts.map(positionOf);
  };
  const sizes = list(cells[10]);
  const starts = list(cells[11]);
  if (count === undefined || count < 1 || sizes.length !== count || starts.length !== count)
    return 'Blocks whose counts are not `block_count`';
  let last = 0;
  for (let b = 0; b < count; b++) {
    const s = starts[b];
    const z = sizes[b];
    if (s === undefined || z === undefined) return 'A block position that is not a whole number';
    if (z < 1) return 'A block of no length';
    if ((b === 0 && s !== 0) || s < last || s + z > length)
      return 'Blocks out of order, overlapping or outside the item';
    last = s + z;
  }
  if (last !== length) return 'Blocks that do not reach the item’s end';
  return [starts as number[], sizes as number[]];
}

/** A cell as its field's type: a number (NaN missing, `.` and empty silently), else its text (null missing). */
function cellValue(
  cell: string | undefined,
  type: ColumnType,
  place: string,
  causes: BedCauses,
): number | string | null {
  if (type === 'category') return cell === undefined || cell === '' || cell === '.' ? null : cell;
  if (cell === undefined || cell === '' || cell === '.') return Number.NaN;
  const v = Number(cell);
  if (Number.isNaN(v) || (type === 'int32' && !Number.isInteger(v))) {
    causes.add('seq.source.bed.cells', 'A number field’s cell is no number: read as missing', place);
    return Number.NaN;
  }
  return v;
}

/** A field's values over the rows `items` as its column. */
function columnOf2(values: readonly (number | string | null)[], type: ColumnType, items: Int32Array): Column {
  if (type === 'category') return categoryColumn(Array.from(items, (i) => values[i] as string | null));
  const out =
    type === 'int32'
      ? new Int32Array(items.length)
      : type === 'float64'
        ? new Float64Array(items.length)
        : new Float32Array(items.length);
  for (let r = 0; r < items.length; r++) {
    const v = values[items[r]] as number;
    out[r] = type === 'int32' && Number.isNaN(v) ? 0 : v;
  }
  return out as Column;
}

/** The pieces meeting `[start, end)` (a zero-length one at its start; none in an empty region), in the items' order. */
export function piecesIn(items: Items, start: number, end: number): number[] {
  const { pStart, pEnd } = items;
  if (end <= start) return [];
  const from = start - items.maxSpan;
  let lo = 0;
  let hi = pStart.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (pStart[mid] < from) lo = mid + 1;
    else hi = mid;
  }
  const out: number[] = [];
  for (let i = lo; i < pStart.length && pStart[i] < end; i++) if (Math.max(pEnd[i], pStart[i] + 1) > start) out.push(i);
  return out;
}

/** The list of the pieces `rows` (of `piecesIn`): identity, strand, fields, and the item's extent on BED12's pieces. */
export function rangesList(items: Items, rows: readonly number[]): RangesList {
  const s = items.schema;
  const item = Int32Array.from(rows, (r) => items.pItem[r]);
  const fields: RangesField[] = s.fields.map((f, k) => ({
    name: f.name,
    column: columnOf2(items.values[k], f.type, item),
  }));
  const strand = items.strand && Int8Array.from(item, (i) => items.strand![i]);
  // the strand a field too, so a layer's data filter and a theme read it (`(item :strand -1)`)
  if (strand) fields.push({ name: 'strand', column: strand });
  return {
    form: 'list',
    start: Int32Array.from(rows, (r) => items.pStart[r]),
    end: Int32Array.from(rows, (r) => items.pEnd[r]),
    ...(strand && { strand }),
    id: categoryColumn(Array.from(item, (i) => items.key[i])),
    ...(fields.length && { fields }),
    // every row its item's extent (a BED12 block's, its line's): a reader marks an item's ends from any tile
    itemStart: Int32Array.from(item, (i) => items.iStart[i]),
    itemEnd: Int32Array.from(item, (i) => items.iEnd[i]),
  };
}

// ---- bedGraph: a system's intervals as values ----

/** A bedGraph system's intervals: sorted by start (file order among equal starts). */
export interface Intervals {
  readonly start: Int32Array;
  readonly end: Int32Array;
  readonly value: Float32Array;
  readonly maxSpan: number;
}

export class IntervalsBuilder {
  readonly start: number[] = [];
  readonly end: number[] = [];
  readonly value: number[] = [];

  /** `chrom start end value`; false when left out (said). */
  add(cells: readonly string[], place: string, causes: BedCauses): boolean {
    const start = positionOf(cells[1]);
    const end = positionOf(cells[2]);
    if (start === undefined || end === undefined) {
      causes.add('seq.source.bed.line', 'A position is not a whole number from 0 to 2^31 − 1', place);
      return false;
    }
    if (end < start) {
      causes.add('seq.source.bed.line', 'An end before its start', place);
      return false;
    }
    const v = cells[3] === undefined ? Number.NaN : Number(cells[3]);
    if (cells.length < 4 || Number.isNaN(v)) {
      causes.add('seq.source.bed.line', 'A bedGraph line without a number in its fourth column', place);
      return false;
    }
    if (cells.length > 4) causes.add('seq.source.bed.cells', 'More columns than the first line: not read', place);
    this.start.push(start);
    this.end.push(end);
    this.value.push(v);
    return true;
  }

  finish(): Intervals {
    const order = sortedOrder(this.start);
    const n = this.start.length;
    const start = new Int32Array(n);
    const end = new Int32Array(n);
    const value = new Float32Array(n);
    let maxSpan = 1;
    for (let i = 0; i < n; i++) {
      const k = order ? order[i] : i;
      start[i] = this.start[k];
      end[i] = this.end[k];
      value[i] = this.value[k];
      maxSpan = Math.max(maxSpan, end[i] - start[i]);
    }
    return { start, end, value, maxSpan };
  }
}

/**
 * The values of `[start, end)` binned at `level` (bins of its width from 0, those meeting the region): `min`, `max`,
 * count-weighted `mean` and `count` of the positions with a value, exactly as `binValues` over `valuesIn` (a later
 * interval over an earlier), without a value a position: the intervals meeting the region painted into segments (each
 * later one cutting what it covers of the earlier), each segment added to its bins whole. Bounded by the bins, not
 * the positions.
 */
export function binsIn(iv: Intervals, start: number, end: number, level: number): ValuesBins {
  const lo = Math.floor(start / level);
  const n = Math.max(0, Math.ceil(end / level) - lo);
  const min = new Float32Array(n).fill(Number.NaN);
  const max = new Float32Array(n).fill(Number.NaN);
  const mean = new Float32Array(n).fill(Number.NaN);
  const count = new Float32Array(n);
  const out: ValuesBins = { form: 'bins', start: lo * level, min, max, mean, count };
  if (end <= start) return out;
  // the effective segments, sorted and apart: `[s, e)` of value `v`, three numbers each
  const seg: number[] = [];
  const tail: number[] = [];
  const from = start - iv.maxSpan;
  let a = 0;
  let b = iv.start.length;
  while (a < b) {
    const mid = (a + b) >> 1;
    if (iv.start[mid] < from) a = mid + 1;
    else b = mid;
  }
  for (let i = a; i < iv.start.length && iv.start[i] < end; i++) {
    if (iv.end[i] <= start) continue;
    const s = Math.max(iv.start[i], start);
    const e = Math.min(iv.end[i], end);
    if (e <= s) continue;
    // what the earlier segments hold past `s`: kept before `s` and after `e`, the new one between
    tail.length = 0;
    while (seg.length && seg[seg.length - 2] > s) {
      const v = seg.pop()!;
      const se = seg.pop()!;
      const ss = seg.pop()!;
      if (se > e) tail.push(Math.max(ss, e), se, v);
      if (ss < s) seg.push(ss, s, v);
    }
    seg.push(s, e, iv.value[i]);
    for (let k = tail.length - 3; k >= 0; k -= 3) seg.push(tail[k], tail[k + 1], tail[k + 2]);
  }
  const sum = new Float64Array(n);
  for (let k = 0; k < seg.length; k += 3) {
    const s = seg[k];
    const e = seg[k + 1];
    const v = seg[k + 2];
    for (let x = Math.floor(s / level); x * level < e; x++) {
      const c = Math.min(e, (x + 1) * level) - Math.max(s, x * level);
      const j = x - lo;
      if (count[j] === 0 || v < min[j]) min[j] = v;
      if (count[j] === 0 || v > max[j]) max[j] = v;
      count[j] += c;
      sum[j] += v * c;
    }
  }
  for (let j = 0; j < n; j++) if (count[j] > 0) mean[j] = sum[j] / count[j];
  return out;
}

/**
 * Values a position over `[start, end)`: from `start` to the last position an interval there covers (NaN between), a
 * later interval over an earlier one; empty when none meets it. The caller bounds `end − start` first.
 */
export function valuesIn(iv: Intervals, start: number, end: number): ValuesPositions {
  if (end <= start) return { form: 'positions', start, value: new Float32Array(0) };
  const from = start - iv.maxSpan;
  let lo = 0;
  let hi = iv.start.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (iv.start[mid] < from) lo = mid + 1;
    else hi = mid;
  }
  let last = start;
  const hits: number[] = [];
  for (let i = lo; i < iv.start.length && iv.start[i] < end; i++) {
    if (iv.end[i] <= start) continue;
    hits.push(i);
    last = Math.max(last, Math.min(iv.end[i], end));
  }
  const value = new Float32Array(last - start).fill(Number.NaN);
  for (const i of hits) value.fill(iv.value[i], Math.max(iv.start[i], start) - start, Math.min(iv.end[i], end) - start);
  return { form: 'positions', start, value };
}
