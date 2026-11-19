/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Column, categoryColumn, issueOf } from '@molstar/core/repr/tree';
import { SourceIssues } from './issues.js';
import {
  abortError,
  type Coverage,
  type FormatRuntime,
  type OpenedSource,
  type RangesField,
  type RegionAnswer,
  type SourceEntry,
} from './protocol.js';
import { SourceOpenError } from './seq-runtime.js';

// An alignment held in memory as a source of the protocol (_spec/README.md "Data", Cells): its rows system (`rows`,
// a position a row in data order) answered as `cells` tiles, a tile's rows × every column; its columns system
// (`columns`) declared for its length. Inline rows read through it (`alignment`); the file readers open their text into
// the same (`openAlignment`).

/** The names an alignment source gives its two systems. */
export const AlignmentSystems = { rows: 'rows', columns: 'columns' } as const;

/** An alignment's rows as plain data: a label and aligned letters a row, every row as wide. */
export interface AlignmentRows {
  readonly label: readonly string[];
  readonly letters: readonly string[];
  /** Each row's first residue's number (EC-39); from its label (`NAME/27-137`) when absent. */
  readonly first?: ArrayLike<number>;
  /**
   * The rows' own fields (Stockholm's `#=GS` tags, an inline table's other columns): a value a row each, `null` (or a
   * missing value) where a row has none; numbers when every value is one, else text.
   */
  readonly fields?: Readonly<Record<string, readonly unknown[]>>;
}

const RangeLabel = /\/(-?\d+)-(-?\d+)$/;

/** The first residue's number a label says (`NAME/27-137`), else 1. */
export function firstNumberOf(label: string): number {
  const m = RangeLabel.exec(label);
  return m ? Number(m[1]) : 1;
}

function malformed(source: SourceEntry, why: string): SourceOpenError {
  return new SourceOpenError([
    issueOf(SourceIssues, 'seq.source.malformed', `The alignment cannot be read: ${why}`, { address: source.address }),
  ]);
}

/**
 * Rows as an opened source: the letters as bytes a cell, row-major, made once; a tile answers its rows' bytes (a view
 * of them, no copy) and its rows table (labels and first numbers). Rows of another width are refused.
 */
export function openAlignment(rows: AlignmentRows, source: SourceEntry = { address: '', format: 'alignment' }) {
  const n = rows.letters.length;
  if (!Array.isArray(rows.label) || rows.label.length !== n) throw malformed(source, 'not one label a row');
  const columns = n ? rows.letters[0].length : 0;
  const codes = new Uint8Array(n * columns);
  for (let i = 0; i < n; i++) {
    const l = rows.letters[i];
    if (typeof l !== 'string' || l.length !== columns) throw malformed(source, `the row ${i} is not ${columns} wide`);
    for (let j = 0; j < columns; j++) {
      const c = l.charCodeAt(j);
      if (c > 127) throw malformed(source, `the row ${i} holds a letter that is not ASCII`);
      codes[i * columns + j] = c;
    }
  }
  const first = new Int32Array(n);
  for (let i = 0; i < n; i++) {
    const f = rows.first?.[i];
    first[i] = f !== undefined ? f : firstNumberOf(rows.label[i]);
  }
  const fields = fieldColumns(rows.fields, n, source);
  return cellsSource(n, columns, (start, end) => ({
    codes: codes.subarray(start * columns, end * columns),
    label: rows.label.slice(start, end),
    first: first.subarray(start, end),
    ...(fields.length
      ? { fields: fields.map((f) => ({ name: f.name, column: sliceColumn(f.column, start, end) })) }
      : {}),
  }));
}

/** The rows' own fields as columns: numbers (NaN missing) when every value is one, else text (`null` missing). */
function fieldColumns(
  fields: AlignmentRows['fields'],
  n: number,
  source: SourceEntry,
): { name: string; column: Column }[] {
  const out: { name: string; column: Column }[] = [];
  if (!fields || typeof fields !== 'object') return out;
  for (const name of Object.keys(fields)) {
    const v = fields[name];
    if (!Array.isArray(v) || v.length !== n) throw malformed(source, `the field ${name} is not one value a row`);
    const numeric = v.every((x) => x === null || x === undefined || typeof x === 'number');
    out.push({
      name,
      column: numeric
        ? Float32Array.from(v, (x) => (typeof x === 'number' ? x : Number.NaN))
        : categoryColumn(v.map((x) => (x === null || x === undefined ? null : String(x)))),
    });
  }
  return out;
}

function sliceColumn(c: Column, start: number, end: number): Column {
  if (c instanceof Float32Array) return c.subarray(start, end);
  const cat = c as ReturnType<typeof categoryColumn>;
  return categoryColumn(Array.from(cat.codes.subarray(start, end), (k) => (k < 0 ? null : cat.dict[k])));
}

/**
 * A source of `cells` over `rows` × `columns` whose tiles `rowsOf(start, end)` makes (a generated alignment makes them
 * on demand: no row is held before it is asked for, EC-41). Level 1 only: coarse row bins are a source's own.
 */
export function cellsSource(
  rows: number,
  columns: number,
  rowsOf: (
    start: number,
    end: number,
  ) => {
    readonly codes: Uint8Array;
    readonly label: readonly string[];
    readonly first: Int32Array;
    readonly fields?: readonly RangesField[];
  },
): OpenedSource {
  const coverage: Coverage = {
    systems: [
      { name: AlignmentSystems.rows, length: rows },
      { name: AlignmentSystems.columns, length: columns },
    ],
    shapes: [{ system: AlignmentSystems.rows, shape: 'cells' }],
    levels: [1],
  };
  return {
    coverage,
    read(request, signal): Promise<RegionAnswer> {
      if (signal.aborted) return Promise.reject(abortError());
      if (request.system !== AlignmentSystems.rows || request.level !== 1 || !request.shapes.includes('cells'))
        return Promise.resolve({
          ok: false,
          kind: 'not-served',
          issue: issueOf(SourceIssues, 'seq.source.not_served', 'An alignment serves its rows’ cells at level 1'),
        });
      const start = Math.max(0, Math.min(request.start, rows));
      const end = Math.max(start, Math.min(request.end, rows));
      const r = rowsOf(start, end);
      return Promise.resolve({
        ok: true,
        tables: {
          cells: {
            form: 'codes',
            start,
            columns,
            codes: r.codes,
            label: categoryColumn(r.label),
            first: r.first,
            ...(r.fields ? { fields: r.fields } : {}),
          },
        },
      });
    },
    dispose() {},
  };
}

/** Inline rows (`value: { label, letters }`) as an alignment source (`format: 'alignment'`). */
export const AlignmentRuntime: FormatRuntime = {
  name: 'alignment',
  formats: ['alignment'],
  async open(source) {
    const v = source.value as AlignmentRows | undefined;
    if (!v || typeof v !== 'object') throw malformed(source, 'no inline rows');
    return openAlignment(v, source);
  },
};
