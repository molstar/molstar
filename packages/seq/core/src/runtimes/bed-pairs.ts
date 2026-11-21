/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Column, type ColumnType, categoryColumn } from '@molstar/core/repr/tree';
import type { PairsList, RangesField } from '../data/protocol.js';
import type { BedCauses } from './bed-issues.js';
import { type DeclaredType, positionOf, strandOf } from './bed-parse.js';

// BEDPE (_spec/README.md "Data", BED family): read whole, then indexed by both ends once (two `Int32Array`
// permutations, each by system then start), so a region finds every pair with either end in it, each once (EC-208).
// What it holds is counted in bytes against the runtime's `maxPairsBytes` before the index is made.

/** BEDPE's columns by place. */
export const BedpeColumns = [
  'chrom1',
  'start1',
  'end1',
  'chrom2',
  'start2',
  'end2',
  'name',
  'score',
  'strand1',
  'strand2',
] as const;

/** A field of a BEDPE: by a mapping (any column past the positions, by number) or `score` by place. */
interface PairsFieldSpec {
  readonly name: string;
  readonly column: number;
  readonly type: ColumnType;
}

export interface PairsSchema {
  readonly columns: number;
  readonly name?: number;
  readonly strand1?: number;
  readonly strand2?: number;
  readonly fields: readonly PairsFieldSpec[];
}

/** A BEDPE's schema from its first line's columns and the source's `fields` and `types` (a column by its number). */
export function pairsSchema(
  columns: number,
  mapping: Readonly<Record<string, string>> = {},
  types: Readonly<Record<string, DeclaredType>> = {},
): PairsSchema | { readonly error: string } {
  if (columns < 6) return { error: `The first line has ${columns} columns; BEDPE has 6 at least` };
  const mapped = new Set<number>();
  const fields: PairsFieldSpec[] = [];
  const type = (name: string, fallback: ColumnType): ColumnType => {
    const d = types[name];
    return d === undefined ? fallback : d === 'string' ? 'category' : d;
  };
  for (const [field, spelling] of Object.entries(mapping)) {
    const c = /^\d+$/.test(spelling) ? Number(spelling) : BedpeColumns.indexOf(spelling as never);
    if (c < 0) return { error: `The field “${field}” reads “${spelling}”, which is no column of BEDPE` };
    if (c < 6) return { error: `The field “${field}” reads column ${c}, a position` };
    if (c >= columns) return { error: `The field “${field}” reads column ${c}; the first line has ${columns}` };
    mapped.add(c);
    fields.push({ name: field, column: c, type: type(field, 'category') });
  }
  const own = new Set(fields.map((f) => f.name));
  const at = (c: number) => (c < columns && !mapped.has(c) && !own.has(BedpeColumns[c]) ? c : undefined);
  const score = at(7);
  if (score !== undefined) fields.push({ name: 'score', column: 7, type: type('score', 'float32') });
  const name = at(6);
  const strand1 = at(8);
  const strand2 = at(9);
  return {
    columns,
    fields,
    ...(name !== undefined && { name }),
    ...(strand1 !== undefined && { strand1 }),
    ...(strand2 !== undefined && { strand2 }),
  };
}

/** A BEDPE read whole: each end's system (an index into `names`), positions and strand; the name; the fields. */
export interface PairsRows {
  readonly names: readonly string[];
  readonly aSys: Uint8Array | Uint16Array | Int32Array;
  readonly bSys: Uint8Array | Uint16Array | Int32Array;
  readonly aStart: Int32Array;
  readonly aEnd: Int32Array;
  readonly bStart: Int32Array;
  readonly bEnd: Int32Array;
  readonly aStrand?: Int8Array;
  readonly bStrand?: Int8Array;
  readonly name?: readonly (string | null)[];
  readonly fields: readonly {
    readonly name: string;
    readonly type: ColumnType;
    readonly values: readonly (number | string | null)[];
  }[];
}

/** Bytes a row holds with its two orders: its columns' and the orders' (an estimate's floor before it is read). */
export function rowBytes(schema: PairsSchema, systems: number): number {
  const sys = systems <= 256 ? 1 : systems <= 65536 ? 2 : 4;
  let n = 2 * sys + 4 * 4 + 2 * 4;
  if (schema.strand1 !== undefined) n++;
  if (schema.strand2 !== undefined) n++;
  if (schema.name !== undefined) n += 4;
  return n + 4 * schema.fields.length;
}

export class PairsBuilder {
  private readonly sysOf = new Map<string, number>();
  readonly names: string[] = [];
  private readonly a: number[] = [];
  private readonly b: number[] = [];
  private readonly pos: number[] = [];
  private readonly strands: number[] = [];
  private readonly name: (string | null)[] = [];
  private readonly values: (number | string | null)[][];

  constructor(readonly schema: PairsSchema) {
    this.values = schema.fields.map(() => []);
  }

  get rows() {
    return this.a.length;
  }

  private system(name: string) {
    let s = this.sysOf.get(name);
    if (s === undefined) {
      s = this.names.length;
      this.sysOf.set(name, s);
      this.names.push(name);
    }
    return s;
  }

  /** Reads one line's cells; false when it is left out (said). */
  add(cells: readonly string[], place: string, causes: BedCauses): boolean {
    const s = this.schema;
    const p = [positionOf(cells[1]), positionOf(cells[2]), positionOf(cells[4]), positionOf(cells[5])];
    if (!cells[0] || !cells[3] || p.some((v) => v === undefined)) {
      causes.add('seq.source.bed.line', 'A BEDPE line without two ends of whole numbers from 0 to 2^31 − 1', place);
      return false;
    }
    if (p[1]! < p[0]! || p[3]! < p[2]!) {
      causes.add('seq.source.bed.line', 'An end before its start', place);
      return false;
    }
    if (cells.length > s.columns)
      causes.add('seq.source.bed.cells', 'More columns than the first line: not read', place);
    const strand = (c: number | undefined) => {
      if (c === undefined) return 0;
      const v = strandOf(cells[c]);
      if (v === undefined) causes.add('seq.source.bed.cells', 'A strand not +, - or .: read as .', place);
      return v ?? 0;
    };
    this.a.push(this.system(cells[0]));
    this.b.push(this.system(cells[3]));
    this.pos.push(p[0]!, p[1]!, p[2]!, p[3]!);
    this.strands.push(strand(s.strand1), strand(s.strand2));
    if (s.name !== undefined) {
      const n = cells[s.name];
      // `.` is no identity (BEDPE writes it for every unnamed pair)
      this.name.push(n === undefined || n === '' || n === '.' ? null : n);
    }
    for (let f = 0; f < s.fields.length; f++) {
      const { column, type } = s.fields[f];
      const cell = cells[column];
      if (type === 'category') this.values[f].push(cell === undefined || cell === '' || cell === '.' ? null : cell);
      else {
        const v = cell === undefined || cell === '' || cell === '.' ? Number.NaN : Number(cell);
        if (Number.isNaN(v) && cell !== undefined && cell !== '' && cell !== '.')
          causes.add('seq.source.bed.cells', 'A number field’s cell is no number: read as missing', place);
        this.values[f].push(v);
      }
    }
    return true;
  }

  finish(): PairsRows {
    const n = this.a.length;
    const sys = (v: number[]) =>
      this.names.length <= 256
        ? Uint8Array.from(v)
        : this.names.length <= 65536
          ? Uint16Array.from(v)
          : Int32Array.from(v);
    const col = (k: number) => Int32Array.from({ length: n }, (_, i) => this.pos[4 * i + k]);
    const s = this.schema;
    return {
      names: this.names,
      aSys: sys(this.a),
      bSys: sys(this.b),
      aStart: col(0),
      aEnd: col(1),
      bStart: col(2),
      bEnd: col(3),
      ...(s.strand1 !== undefined && { aStrand: Int8Array.from({ length: n }, (_, i) => this.strands[2 * i]) }),
      ...(s.strand2 !== undefined && { bStrand: Int8Array.from({ length: n }, (_, i) => this.strands[2 * i + 1]) }),
      ...(s.name !== undefined && { name: this.name }),
      fields: s.fields.map((f, k) => ({ name: f.name, type: f.type, values: this.values[k] })),
    };
  }
}

/** One end's order: rows by system then start, each system's run (`first[s]` to `first[s + 1]`) and longest end. */
export interface EndOrder {
  readonly order: Int32Array;
  readonly first: Int32Array;
  readonly span: Int32Array;
}

/** Both ends' orders of a BEDPE. */
export interface PairsIndex {
  readonly a: EndOrder;
  readonly b: EndOrder;
}

/** Keys exact in a float64: a start (< 2^31) × 2^21 + a row (< 2^21). */
const RowBits = 2 ** 21;

/**
 * Rows by system (a counting sort) then start (each system's starts × 2^21 + the row, sorted as float64: exact, and
 * native; a file past 2^21 rows by a comparison sort), row order among equal starts.
 */
function endOrder(sys: ArrayLike<number>, start: Int32Array, end: Int32Array, systems: number): EndOrder {
  const n = start.length;
  const first = new Int32Array(systems + 1);
  for (let i = 0; i < n; i++) first[sys[i] + 1]++;
  for (let s = 0; s < systems; s++) first[s + 1] += first[s];
  const order = new Int32Array(n);
  const at = first.slice(0, systems);
  const span = new Int32Array(systems).fill(1);
  for (let i = 0; i < n; i++) {
    const s = sys[i];
    order[at[s]++] = i;
    span[s] = Math.max(span[s], end[i] - start[i]);
  }
  for (let s = 0; s < systems; s++) {
    const lo = first[s];
    const hi = first[s + 1];
    if (hi - lo < 2) continue;
    if (n < RowBits) {
      const keys = new Float64Array(hi - lo);
      for (let k = lo; k < hi; k++) keys[k - lo] = start[order[k]] * RowBits + order[k];
      keys.sort();
      for (let k = lo; k < hi; k++) order[k] = keys[k - lo] % RowBits;
    } else {
      order.subarray(lo, hi).sort((x, y) => start[x] - start[y] || x - y);
    }
  }
  return { order, first, span };
}

/** The two orders of a BEDPE's rows (`maxPairsBytes` checked by the caller first). */
export function indexPairs(rows: PairsRows): PairsIndex {
  const systems = rows.names.length;
  return {
    a: endOrder(rows.aSys, rows.aStart, rows.aEnd, systems),
    b: endOrder(rows.bSys, rows.bStart, rows.bEnd, systems),
  };
}

/** Whether `[s, e)` (a zero-length one at its start) meets `[start, end)`. */
const meets = (s: number, e: number, start: number, end: number) => s < end && Math.max(e, s + 1) > start;

function scan(
  o: EndOrder,
  sys: number,
  starts: Int32Array,
  ends: Int32Array,
  start: number,
  end: number,
  visit: (row: number) => void,
) {
  let lo = o.first[sys];
  let hi = o.first[sys + 1];
  const from = start - o.span[sys];
  const stop = hi;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (starts[o.order[mid]] < from) lo = mid + 1;
    else hi = mid;
  }
  for (let k = lo; k < stop; k++) {
    const r = o.order[k];
    if (starts[r] >= end) break;
    if (meets(starts[r], ends[r], start, end)) visit(r);
  }
}

/**
 * The rows with an end meeting `[start, end)` of system `sys`, each once (a pair with both ends there is found from its
 * first), in row order (EC-208).
 */
export function pairsIn(rows: PairsRows, index: PairsIndex, sys: number, start: number, end: number): Int32Array {
  const found: number[] = [];
  if (end <= start) return new Int32Array(0);
  scan(index.a, sys, rows.aStart, rows.aEnd, start, end, (r) => found.push(r));
  scan(index.b, sys, rows.bStart, rows.bEnd, start, end, (r) => {
    if (!(rows.aSys[r] === sys && meets(rows.aStart[r], rows.aEnd[r], start, end))) found.push(r);
  });
  return Int32Array.from(found).sort();
}

function gathered(values: readonly (number | string | null)[], type: ColumnType, rows: Int32Array): Column {
  if (type === 'category') return categoryColumn(Array.from(rows, (r) => values[r] as string | null));
  const out =
    type === 'int32'
      ? new Int32Array(rows.length)
      : type === 'float64'
        ? new Float64Array(rows.length)
        : new Float32Array(rows.length);
  for (let i = 0; i < rows.length; i++) {
    const v = values[rows[i]] as number;
    out[i] = type === 'int32' && Number.isNaN(v) ? 0 : v;
  }
  return out as Column;
}

/** The pairs list of `found` rows: both ends' systems, positions and strands, `id` the name (`.` none), `row`. */
export function pairsList(rows: PairsRows, found: Int32Array): PairsList {
  const fields: RangesField[] = rows.fields.map((f) => ({ name: f.name, column: gathered(f.values, f.type, found) }));
  const pick = (c: Int32Array) => Int32Array.from(found, (r) => c[r]);
  return {
    form: 'list',
    aSystem: categoryColumn(Array.from(found, (r) => rows.names[rows.aSys[r]])),
    bSystem: categoryColumn(Array.from(found, (r) => rows.names[rows.bSys[r]])),
    aStart: pick(rows.aStart),
    aEnd: pick(rows.aEnd),
    bStart: pick(rows.bStart),
    bEnd: pick(rows.bEnd),
    ...(rows.aStrand && { aStrand: Int8Array.from(found, (r) => rows.aStrand![r]) }),
    ...(rows.bStrand && { bStrand: Int8Array.from(found, (r) => rows.bStrand![r]) }),
    ...(rows.name && { id: categoryColumn(Array.from(found, (r) => rows.name![r])) }),
    row: found.slice(),
    ...(fields.length && { fields }),
  };
}
