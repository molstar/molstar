/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { decodeMsgPack } from '@molstar/core/data/msgpack';
import { readDelimited } from '@molstar/core/data/text/delimited';
import { msgPackColumnDecoder } from '@molstar/core/repr/codecs/msgpack';
import {
  CategoryColumn,
  type Column,
  type ColumnType,
  categoryColumn,
  columnType,
  defineIssues,
  type Issue,
  issueOf,
  type Value,
} from '@molstar/core/repr/tree';
import { SourceIssues } from './issues.js';
import type {
  FieldInfo,
  FormatRuntime,
  LettersTable,
  OpenContext,
  OpenedSource,
  OpenLimits,
  RangesField,
  RangesList,
  Shape,
  ShapeInfo,
  SourceEntry,
  SystemInfo,
  Tables,
  ValuesPositions,
} from './protocol.js';
import { isAbortError } from './protocol.js';
import type { SeqData, SeqPair } from './seq-format.js';
import { openSeq, SourceOpenError } from './seq-runtime.js';
import { DefaultMaxBlockValues, DefaultMaxColumnLength } from './tables.js';

// The records runtime (_spec/README.md "Data", Records): a file of records (TSV, CSV, JSON, MessagePack) or an inline
// value, read whole at open, its conventions undone by the source's transform (`fields`, `records`, `letters`,
// `system`, `base`, `ends`, `types`) into `seq` data: a pair a system over the whole system at level 1, answered by
// region by the `seq` runtime's join. An entry of its own (`@molstar/seqx-core/data/records`): a reader of a spec never
// pays for it.

/** The records runtime's codes: one issue a cause, with its count and its first places. */
export const RecordsIssues = defineIssues('seq.source.records', {
  'seq.source.records.malformed': {
    severity: 'error',
    description:
      'The file cannot be read as its format (no header, a document that does not parse, records of neither form, columns of unequal lengths, a column named twice, nested too deep); the source is unavailable',
  },
  'seq.source.records.field': {
    severity: 'error',
    description: 'A field’s column or path is found in no record; the source is unavailable',
    hint: 'Name one of the columns there are in the source’s `fields`',
  },
  'seq.source.records.value': {
    severity: 'warning',
    description: 'A position or value is not a number, or an end is past int32; the row is left out',
  },
  'seq.source.records.cell': {
    severity: 'info',
    description:
      'A cell its field cannot hold, read as missing (an object or a list, a number that is none, a strand read as 0), or a position named twice (the last kept); the row is kept',
  },
  'seq.source.records.missing': {
    severity: 'warning',
    description: 'A row’s position is missing; the row is left out',
  },
  'seq.source.records.order': {
    severity: 'warning',
    description: 'A row ends before it starts; the row is left out',
  },
  'seq.source.records.range': {
    severity: 'warning',
    description:
      'A row starts outside its system (left out), or runs past the letters’ end (clipped); both counted here',
  },
  'seq.source.records.cells': {
    severity: 'warning',
    description: 'A row has another count of cells than the header: missing cells empty, extra ones dropped',
  },
  'seq.source.records.size': {
    severity: 'error',
    description:
      'A matrix is not the size its axes’ total says, or its rows are of unequal lengths; the source is unavailable',
    hint: 'A matrix has one row a position of the first axis and one column a position of the second',
  },
  'seq.source.records.encoding': {
    severity: 'info',
    description: 'The text holds replacement characters (U+FFFD): it may not be UTF-8',
    hint: 'Save the file as UTF-8: Latin-1 accents read as replacement characters',
  },
});

type Code = keyof typeof RecordsIssues.codes;

/** The formats the runtime reads. */
export const RecordsFormats = ['tsv', 'csv', 'json', 'msgpack', 'gff3'] as const;
export type RecordsFormat = (typeof RecordsFormats)[number];

/** A field's declared type; `string` is held as a category. */
export type RecordsFieldType = 'int32' | 'float32' | 'float64' | 'category' | 'string';

/** The transform: the `source` node's params that say what the file's columns mean (language 0.4.0). */
export interface RecordsTransform {
  readonly fields?: Readonly<Record<string, string>>;
  readonly records?: string;
  readonly letters?: string;
  readonly system?: string;
  readonly base?: 0 | 1;
  readonly ends?: 'inclusive' | 'exclusive';
  readonly types?: Readonly<Record<string, RecordsFieldType>>;
}

/** What is read: a file's text (TSV, CSV, JSON), its bytes (MessagePack), or an inline value (base 0, exclusive). */
export type RecordsInput =
  | { readonly format: 'tsv' | 'csv' | 'json' | 'gff3'; readonly text: string }
  | { readonly format: 'msgpack'; readonly bytes: Uint8Array }
  | { readonly value: Value };

export interface RecordsOptions {
  /** The source's own system's name: `transform.system`, else this (the source node's id). */
  readonly system: string;
  readonly limits: OpenLimits;
}

/** The roles: the protocol's fields, typed by it. */
const Roles = ['target', 'start', 'end', 'position', 'strand', 'id', 'value'] as const;
type Role = (typeof Roles)[number];
const RoleSet = new Set<string>(Roles);

/** A path has at most this many segments (bounded work). */
const MaxPathSegments = 32;
const MaxInt32 = 2 ** 31 - 1;

/** A map's own entry: a transform's keys are any names (`constructor`, `__proto__`), never the prototype's. */
function own<T>(map: Readonly<Record<string, T>> | undefined, key: string): T | undefined {
  return map !== undefined && Object.hasOwn(map, key) ? map[key] : undefined;
}
/** Issues keep the first this many places. */
const MaxPlaces = 10;

class Refused extends Error {
  constructor(readonly issue: Issue) {
    super(issue.message);
  }
}

function refuse(code: Code, message: string, hint?: string): never {
  throw new Refused(issueOf(RecordsIssues, code, message, hint ? { hint } : {}));
}

function tooLarge(message: string, hint: string): never {
  throw new Refused(issueOf(SourceIssues, 'seq.source.too_large', message, { hint }));
}

// ---- documents: records of three kinds behind one reader ----

/** A field's cells: text (a table's) or values (a document's), with the typed column it is when it is one. */
interface Cells {
  get(i: number): unknown;
  readonly text?: boolean;
  readonly column?: Column;
}

interface Records {
  readonly count: number;
  /** Where record `i` is, for an issue: a line, or a record's index. */
  place(i: number): string;
  /** The cells of a column name (a table) or a dotted path (a document); undefined when found in no record. */
  cells(spelling: string): Cells | undefined;
  /** The fields read by default: the scalar members or columns, in order. */
  readonly scalars: readonly string[];
  /** What there is, for an error that names a field found nowhere. */
  readonly names: readonly string[];
  /** A text table: positions and numbers parse past spaces; other cells are categories. */
  readonly text: boolean;
}

const isPlainObject = (v: unknown): v is Record<string, unknown> =>
  typeof v === 'object' && v !== null && !Array.isArray(v) && !ArrayBuffer.isView(v) && !(v instanceof CategoryColumn);

const isScalar = (v: unknown) => v === null || (typeof v !== 'object' && typeof v !== 'function');

function segments(path: string): string[] {
  const s = path.split('.');
  if (s.length > MaxPathSegments)
    refuse('seq.source.records.malformed', `A path of more than ${MaxPathSegments} segments: ${path}`);
  return s;
}

/** A member (own only: `__proto__` and `constructor` are ordinary names), or an index when all digits. */
function step(at: unknown, segment: string): unknown {
  if (Array.isArray(at)) return /^\d+$/.test(segment) ? at[Number(segment)] : undefined;
  if (isPlainObject(at)) return Object.hasOwn(at, segment) ? at[segment] : undefined;
  return undefined;
}

function walk(at: unknown, path: readonly string[]): unknown {
  for (const s of path) {
    at = step(at, s);
    if (at === undefined) return undefined;
  }
  return at;
}

function tableRecords(text: string, delimiter: string, causes: Causes): Records {
  const read = readDelimited(text, { delimiter });
  if (!read.ok) refuse('seq.source.records.malformed', `The table does not read: ${read.message} (line ${read.line})`);
  const [header, ...rows] = read.records;
  if (!header) refuse('seq.source.records.malformed', 'The table has no header');
  const index = new Map<string, number>();
  header.forEach((name, k) => {
    if (index.has(name)) refuse('seq.source.records.malformed', `The header names the column “${name}” twice`);
    index.set(name, k);
  });
  const lines = read.lines.slice(1);
  for (let i = 0; i < rows.length; i++)
    if (rows[i].length !== header.length)
      causes.add(
        'seq.source.records.cells',
        `A row has another count of cells than the header's ${header.length}: missing cells empty, extra ones dropped`,
        `line ${lines[i]}`,
      );
  return {
    count: rows.length,
    place: (i) => `line ${lines[i]}`,
    cells(name) {
      const k = index.get(name);
      return k === undefined ? undefined : { text: true, get: (i) => rows[i][k] ?? '' };
    },
    scalars: header,
    names: header,
    text: true,
  };
}

/** The GFF3 preset's columns, in file order (`seqid` the `target`, `ID` the `id`). */
const Gff3Columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase'] as const;

/** A GFF3 value: percent-decoded (a malformed escape is left as written). */
function gffDecode(v: string): string {
  if (!v.includes('%')) return v;
  try {
    return decodeURIComponent(v);
  } catch {
    return v;
  }
}

/**
 * GFF3 as records: nine tab-separated columns, `#` lines and `##` directives read past, everything after `##FASTA`
 * ignored, `.` missing; the attributes (`key=value;…`) are columns by their key, `seqid` is the `target` and `ID` the
 * `id` (_spec/README.md "Data", Records).
 */
function gff3Records(text: string, causes: Causes): Records {
  const rows: (string | null)[][] = [];
  const attrs: Record<string, string>[] = [];
  const lines: number[] = [];
  const keys = new Set<string>();
  let line = 0;
  let at = 0;
  const n = text.length;
  while (at < n) {
    let end = text.indexOf('\n', at);
    if (end < 0) end = n;
    line++;
    let l = text.slice(at, end);
    at = end + 1;
    if (l.endsWith('\r')) l = l.slice(0, -1);
    if (line === 1 && l.charCodeAt(0) === 0xfeff) l = l.slice(1);
    if (l.trim() === '') continue;
    if (l.startsWith('##FASTA')) break;
    if (l.charCodeAt(0) === 35) continue;
    const cells = l.split('\t');
    if (cells.length !== 9)
      causes.add(
        'seq.source.records.cells',
        'A row has another count of cells than the nine columns: missing cells empty, extra ones dropped',
        `line ${line}`,
      );
    const row: (string | null)[] = [];
    for (let k = 0; k < 8; k++) {
      const c = cells[k] ?? '.';
      row.push(c === '.' || c === '' ? null : gffDecode(c));
    }
    const a: Record<string, string> = Object.create(null);
    const raw = cells[8] ?? '.';
    if (raw !== '.' && raw !== '')
      for (const part of raw.split(';')) {
        const eq = part.indexOf('=');
        if (eq <= 0) continue;
        const key = part.slice(0, eq).trim();
        if (!key) continue;
        // an attribute named as a column leaves the column its own; a key written twice keeps the last
        if ((Gff3Columns as readonly string[]).includes(key)) {
          causes.add(
            'seq.source.records.cell',
            'An attribute is named as a column of the format: the column is the field, the attribute left out',
            `line ${line}`,
          );
          continue;
        }
        if (key in a)
          causes.add(
            'seq.source.records.cell',
            'An attribute key is written twice in a row: the last kept',
            `line ${line}`,
          );
        a[key] = gffDecode(part.slice(eq + 1));
        keys.add(key);
      }
    rows.push(row);
    attrs.push(a);
    lines.push(line);
  }
  const attributeNames = [...keys];
  const names = [...Gff3Columns, ...attributeNames];
  const aliases: Readonly<Record<string, string>> = { target: 'seqid', id: 'ID' };
  return {
    count: rows.length,
    place: (i) => `line ${lines[i]}`,
    cells(name) {
      const k = (Gff3Columns as readonly string[]).indexOf(own(aliases, name) ?? name);
      if (k >= 0) return { text: true, get: (i) => rows[i][k] };
      const key = own(aliases, name) ?? name;
      if (!keys.has(key)) return undefined;
      return { text: true, get: (i) => attrs[i][key] ?? null };
    },
    // `seqid` and `ID` are the roles `target` and `id`, not fields under their own names
    scalars: names.filter((k) => k !== 'seqid' && k !== 'ID'),
    names,
    text: true,
  };
}

function aosRecords(rows: readonly unknown[]): Records {
  rows.forEach((r, i) => {
    if (!isPlainObject(r)) refuse('seq.source.records.malformed', `Record ${i} is not an object`);
  });
  const objects = rows as readonly Record<string, unknown>[];
  // the scalar members: those no record gives an object or a list
  const seen = new Map<string, boolean>();
  for (const r of objects)
    for (const k of Object.keys(r)) {
      const scalar = isScalar(r[k]);
      const was = seen.get(k);
      seen.set(k, was === undefined ? scalar : was && scalar);
    }
  return {
    count: objects.length,
    place: (i) => `record ${i}`,
    cells(path) {
      const p = segments(path);
      let found = false;
      for (const r of objects)
        if (walk(r, p) !== undefined) {
          found = true;
          break;
        }
      return found ? { get: (i) => walk(objects[i], p) } : undefined;
    },
    scalars: [...seen].filter(([, s]) => s).map(([k]) => k),
    names: objects.length ? Object.keys(objects[0]) : [],
    text: false,
  };
}

type ColumnLike = readonly unknown[] | Column;

const isColumnLike = (v: unknown): v is ColumnLike => Array.isArray(v) || columnType(v) !== undefined;

function soaRecords(columns: Record<string, unknown>): Records {
  const names = Object.keys(columns);
  let count: number | undefined;
  for (const k of names) {
    const c = columns[k];
    if (!isColumnLike(c)) refuse('seq.source.records.malformed', `The column “${k}” is not a list`);
    if (count === undefined) count = c.length;
    else if (c.length !== count)
      refuse(
        'seq.source.records.malformed',
        `The columns are of unequal lengths (“${k}” has ${c.length}, not ${count})`,
      );
  }
  const cellsOf = (c: ColumnLike): Cells =>
    c instanceof CategoryColumn
      ? { get: (i) => c.get(i), column: c }
      : Array.isArray(c)
        ? { get: (i) => c[i] }
        : { get: (i) => (c as ArrayLike<number>)[i], column: c as Column };
  return {
    count: count ?? 0,
    place: (i) => `record ${i}`,
    cells(path) {
      const c = walk(columns, segments(path));
      return c !== undefined && isColumnLike(c) && c.length === count ? cellsOf(c) : undefined;
    },
    scalars: names.filter((k) => {
      const c = columns[k] as ColumnLike;
      if (!Array.isArray(c)) return true;
      for (const v of c) if (!isScalar(v)) return false;
      return true;
    }),
    names,
    text: false,
  };
}

/** The records a document's `records` path reaches: an array of records (AoS) or an object of columns (SoA). */
function documentRecords(root: unknown, path: string | undefined): Records {
  const at = path === undefined || path === '' ? root : walk(root, segments(path));
  if (Array.isArray(at)) return aosRecords(at);
  if (isPlainObject(at)) return soaRecords(at);
  return refuse(
    'seq.source.records.malformed',
    path ? `The path “${path}” reaches neither records nor columns` : 'The document holds neither records nor columns',
  );
}

function parseDocument(input: RecordsInput): unknown {
  try {
    if ('value' in input) return input.value;
    if (input.format === 'msgpack') return decodeMsgPack(input.bytes, msgPackColumnDecoder);
    return JSON.parse(input.text);
  } catch (e) {
    // both parsers recurse: a document deeper than the stack is malformed, never a crash
    const deep = e instanceof RangeError && /call stack/i.test(e.message);
    return refuse(
      'seq.source.records.malformed',
      deep ? 'The document is nested too deeply' : `The document does not parse: ${(e as Error)?.message ?? e}`,
    );
  }
}

// ---- cells to values ----

const Decimal = /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/;

/** A position: an integer, or text of one (spaces around it read past); `undefined` missing, `NaN` not a number. */
function positionOf(cell: unknown): number | undefined {
  if (cell === null || cell === undefined) return undefined;
  let n: number;
  if (typeof cell === 'number') n = cell;
  else if (typeof cell === 'string') return integerText(cell);
  else return Number.NaN;
  return Number.isInteger(n) && Math.abs(n) <= MaxInt32 ? n : Number.NaN;
}

const isSpace = (c: number) => c === 32 || c === 9 || c === 10 || c === 13;

/**
 * An integer's text, spaces around it read past (by hand, not a pattern: the open's hot loop reads two a row); `undefined`
 * blank, `NaN` anything else or past int32.
 */
function integerText(s: string): number | undefined {
  let i = 0;
  let j = s.length;
  while (i < j && isSpace(s.charCodeAt(i))) i++;
  while (j > i && isSpace(s.charCodeAt(j - 1))) j--;
  if (i === j) return undefined;
  let sign = 1;
  const c0 = s.charCodeAt(i);
  if (c0 === 45 || c0 === 43) {
    if (c0 === 45) sign = -1;
    i++;
    if (i === j) return Number.NaN;
  }
  let n = 0;
  for (; i < j; i++) {
    const d = s.charCodeAt(i) - 48;
    if (d < 0 || d > 9) return Number.NaN;
    n = n * 10 + d;
    if (n > MaxInt32) return Number.NaN;
  }
  return sign * n;
}

/** A number: one, or text of one; `undefined` missing, `NaN` not a number. */
function numberOf(cell: unknown): number | undefined {
  if (cell === null || cell === undefined) return undefined;
  if (typeof cell === 'number') return cell;
  if (typeof cell === 'string') {
    const s = cell.trim();
    if (s === '') return undefined;
    return Decimal.test(s) ? Number(s) : Number.NaN;
  }
  return Number.NaN;
}

const Strands: Readonly<Record<string, number>> = { '+': 1, '-': -1, '.': 0, '?': 0, '': 0, '1': 1, '-1': -1, '0': 0 };

/** A strand −1, 0, 1; undefined for anything else (read as 0, said once). */
function strandOf(cell: unknown): number | undefined {
  if (cell === null || cell === undefined) return 0;
  if (typeof cell === 'number') return cell === 1 || cell === -1 || cell === 0 ? cell : undefined;
  if (typeof cell === 'string') return Object.hasOwn(Strands, cell.trim()) ? Strands[cell.trim()] : undefined;
  return undefined;
}

/** A cell as text: a scalar's spelling, `null` for none; an object or a list is `false` (said once). */
/** A rows system's name: an alignment's rows (`/data/aln/(rows)`). */
const RowsSystemName = /\/\(rows\)$/;

/** Whether every record's system (its `target`, else the source's own) is a rows system. */
function onRows(target: Cells | undefined, ownSystem: string, count: number): boolean {
  if (!target) return RowsSystemName.test(ownSystem);
  for (let i = 0; i < count; i++) {
    const t = textOf(target.get(i));
    if (!RowsSystemName.test(t ? t : ownSystem)) return false;
  }
  return true;
}

function textOf(cell: unknown): string | null | false {
  if (cell === null || cell === undefined) return null;
  if (typeof cell === 'string') return cell;
  if (typeof cell === 'number' || typeof cell === 'boolean' || typeof cell === 'bigint') return String(cell);
  return false;
}

// ---- issues ----

interface Cause {
  code: Code;
  message: string;
  count: number;
  places: string[];
  hint?: string;
}

/** The causes said, one issue each; a code's severity is its catalog's (EC-87). */
class Causes {
  readonly byKey = new Map<string, Cause>();
  /** A record's place by its index, spelled only for the places kept. */
  placeOf: (i: number) => string = (i) => `record ${i}`;

  /** A cause at a place, or at record `at`. */
  add(code: Code, message: string, at: string | number) {
    const key = `${code}|${message}`;
    let c = this.byKey.get(key);
    if (!c) {
      c = { code, message, count: 0, places: [] };
      this.byKey.set(key, c);
    }
    c.count++;
    if (c.places.length < MaxPlaces) c.places.push(typeof at === 'number' ? this.placeOf(at) : at);
  }

  hint(code: Code, message: string, hint: string) {
    const c = this.byKey.get(`${code}|${message}`);
    if (c) c.hint = hint;
  }

  issues(): Issue[] {
    return [...this.byKey.values()].map((c) => {
      const more = c.count > c.places.length ? ', …' : '';
      return issueOf(
        RecordsIssues,
        c.code,
        `${c.message}: ${c.count} ${c.count === 1 ? 'row' : 'rows'} (${c.places.join(', ')}${more})`,
        c.hint ? { hint: c.hint } : {},
      );
    });
  }
}

// ---- the transform ----

interface Field {
  readonly name: string;
  readonly cells: Cells;
  readonly type: ColumnType;
}

const DeclaredTypes: Readonly<Record<RecordsFieldType, ColumnType>> = {
  int32: 'int32',
  float32: 'float32',
  float64: 'float64',
  category: 'category',
  string: 'category',
};

/** A field's type: declared, else its column's, else a text cell's `category`, a number's `float64`. */
function fieldType(records: Records, cells: Cells, declared: RecordsFieldType | undefined): ColumnType {
  const type = declared === undefined ? undefined : own(DeclaredTypes, declared);
  if (type) return type;
  if (cells.column) return columnType(cells.column)!;
  if (cells.text) return 'category';
  for (let i = 0; i < records.count; i++) {
    const v = cells.get(i);
    if (v === null || v === undefined) continue;
    return typeof v === 'number' ? 'float64' : 'category';
  }
  return 'category';
}

const Constructors = {
  int8: Int8Array,
  uint8: Uint8Array,
  int16: Int16Array,
  uint16: Uint16Array,
  int32: Int32Array,
  uint32: Uint32Array,
  float32: Float32Array,
  float64: Float64Array,
} as const;

/** Whether `rows` is every row in order: a typed column is then kept as it is. */
function isIdentity(rows: Int32Array, count: number): boolean {
  if (rows.length !== count) return false;
  for (let i = 0; i < rows.length; i++) if (rows[i] !== i) return false;
  return true;
}

/** A field's column for `rows`: its typed column kept when it is the type and every row in order, else one copy. */
function fieldColumn(f: Field, rows: Int32Array, identity: boolean, causes: Causes, records: Records): Column {
  const c = f.cells.column;
  if (c && columnType(c) === f.type) {
    if (identity) return c;
    if (c instanceof CategoryColumn) return categoryColumn(Array.from(rows, (r) => c.get(r)));
    const out = new Constructors[f.type as keyof typeof Constructors](rows.length);
    for (let i = 0; i < rows.length; i++) out[i] = (c as ArrayLike<number>)[rows[i]];
    return out;
  }
  if (f.type === 'category') {
    const values: (string | null)[] = new Array(rows.length);
    for (let i = 0; i < rows.length; i++) {
      const t = textOf(f.cells.get(rows[i]));
      if (t === false)
        causes.add('seq.source.records.cell', `The field “${f.name}” reaches an object or a list: null`, rows[i]);
      values[i] = t === false ? null : t;
    }
    return categoryColumn(values);
  }
  const out = new Constructors[f.type as keyof typeof Constructors](rows.length);
  const integer = f.type !== 'float32' && f.type !== 'float64';
  for (let i = 0; i < rows.length; i++) {
    const cell = f.cells.get(rows[i]);
    const n = integer ? positionOf(cell) : numberOf(cell);
    if (n !== undefined && Number.isNaN(n) && !isNaNCell(cell))
      causes.add(
        'seq.source.records.cell',
        `The field “${f.name}” holds a cell that is not ${integer ? 'an integer' : 'a number'}: missing`,
        rows[i],
      );
    out[i] = n === undefined || Number.isNaN(n) ? (integer ? 0 : Number.NaN) : n;
  }
  return out;
}

/** A NaN written as such (MessagePack's NaN, the text `NaN`): a missing value, not a misfit. */
function isNaNCell(cell: unknown): boolean {
  return (typeof cell === 'number' && Number.isNaN(cell)) || (typeof cell === 'string' && /^\s*nan\s*$/i.test(cell));
}

interface SystemRows {
  readonly name: string;
  /** Rows kept for ranges, file order, with their converted start and end. */
  readonly ranges: number[];
  readonly start: number[];
  readonly end: number[];
  /** A row's end was clipped to the letters: its typed column is not the protocol's as given. */
  clipped: boolean;
  /** Values by position (sparse) or in row order. */
  readonly values: Map<number, number>;
  valueRows: number;
}

const OutsideLetters =
  'A row is outside the letters: left out when it starts past their end, clipped when it runs past it';
const BeforeFirst = 'A row starts before the first position';

/**
 * Records read through the transform into `seq` data, or the issues that refuse them (`seq.source.records.malformed`,
 * `.field`, `seq.source.too_large`). Never throws on its input.
 */
export function readRecords(
  input: RecordsInput,
  transform: RecordsTransform,
  options: RecordsOptions,
): { data?: SeqData; issues: Issue[] } {
  const causes = new Causes();
  const issues: Issue[] = [];
  try {
    const file = !('value' in input);
    // the GFF3 preset fixes its own conventions: base 1, ends inclusive
    const gff = 'format' in input && input.format === 'gff3';
    const base = gff ? 1 : (transform.base ?? (file ? 1 : 0));
    const inclusive = gff || (transform.ends ?? (file ? 'inclusive' : 'exclusive')) === 'inclusive';
    const ownSystem = transform.system ?? options.system;

    if ('text' in input && input.text.includes('�'))
      issues.push(
        issueOf(RecordsIssues, 'seq.source.records.encoding', 'The text holds replacement characters (U+FFFD)', {
          hint: RecordsIssues.codes['seq.source.records.encoding'].hint,
        }),
      );
    let records: Records;
    let root: unknown;
    if ('text' in input && input.format === 'gff3') records = gff3Records(input.text, causes);
    else if ('text' in input && input.format !== 'json')
      records = tableRecords(input.text, input.format === 'tsv' ? '\t' : ',', causes);
    else {
      root = parseDocument(input);
      records = documentRecords(root, transform.records);
    }
    causes.placeOf = records.place;
    // bounds before any column is built
    if (records.count > options.limits.maxRows)
      tooLarge(
        `The source holds ${records.count} rows, more than ${options.limits.maxRows}`,
        'Use a format read by region',
      );

    // the letters of the source's own system
    let letters: string | undefined;
    if (transform.letters !== undefined) {
      if (records.text)
        refuse('seq.source.records.field', 'A table holds no letters: `letters` is a path in a document');
      const l = walk(root, segments(transform.letters));
      if (l === undefined)
        refuse('seq.source.records.field', `The letters’ path “${transform.letters}” reaches nothing`);
      if (typeof l !== 'string')
        refuse('seq.source.records.malformed', `The letters’ path “${transform.letters}” reaches no string`);
      letters = l;
    }

    // the roles and the fields
    const map = transform.fields ?? {};
    const used = new Set<string>(Object.values(map));
    const there = () => records.names.slice(0, 20).join(', ') || 'none';
    const missingField = (name: string, spelling: string): never =>
      refuse(
        'seq.source.records.field',
        `The field “${name}” reads “${spelling}”, found in no record; there are: ${there()}`,
        RecordsIssues.codes['seq.source.records.field'].hint,
      );
    const roles: Partial<Record<Role, Cells>> = {};
    for (const role of Roles) {
      const spelling = own(map, role);
      if (spelling !== undefined) roles[role] = records.cells(spelling) ?? missingField(role, spelling);
      else if (!used.has(role)) {
        const c = records.cells(role);
        if (c) roles[role] = c;
      }
    }
    const fields: Field[] = [];
    const named = new Set<string>();
    for (const [name, spelling] of Object.entries(map)) {
      if (RoleSet.has(name)) continue;
      const cells = records.cells(spelling) ?? missingField(name, spelling);
      fields.push({ name, cells, type: fieldType(records, cells, own(transform.types, name)) });
      named.add(name);
    }
    for (const name of records.scalars) {
      if (RoleSet.has(name) || named.has(name) || used.has(name)) continue;
      const cells = records.cells(name);
      if (cells) fields.push({ name, cells, type: fieldType(records, cells, own(transform.types, name)) });
    }

    const hasRanges = roles.start !== undefined || roles.position !== undefined;
    if (roles.start && !roles.end && !roles.position) missingField('end', own(map, 'end') ?? 'end');
    if (roles.end && !roles.start && !roles.position) missingField('start', own(map, 'start') ?? 'start');
    // keyed records: an `id` and no position, each record one at its place in the file, only on an alignment's rows
    // system (`…/(rows)`: its fields, joined by the row's key; _spec/README.md "Data", Records); elsewhere a record
    // needs a position, as before
    const keyed =
      !hasRanges &&
      !roles.value &&
      letters === undefined &&
      roles.id !== undefined &&
      onRows(roles.target, ownSystem, records.count);
    if (!hasRanges && !roles.value && letters === undefined && !keyed)
      refuse(
        'seq.source.records.field',
        `The records have no \`start\`, \`position\` or \`value\`; there are: ${there()}`,
        RecordsIssues.codes['seq.source.records.field'].hint,
      );

    // rows by system: a row's target, else the source's own
    const systems = new Map<string, SystemRows>();
    const systemOf = (name: string) => {
      let s = systems.get(name);
      if (!s) {
        s = { name, ranges: [], start: [], end: [], clipped: false, values: new Map(), valueRows: 0 };
        systems.set(name, s);
      }
      return s;
    };
    // the source's own system when it holds the letters, a row names no target, or no row names any: so a file whose
    // rows all have targets lists no empty system of its own
    if (letters !== undefined) systemOf(ownSystem);
    // values a row from the first position without `position`; with `start` and `end`, `value` is a ranges field only
    const hasValues = roles.value !== undefined && (!hasRanges || roles.position !== undefined);
    const L = letters?.length;
    let pastByOne = 0;
    let pastMore = 0;
    let zeroStart = false;
    for (let i = 0; i < records.count; i++) {
      let target = ownSystem;
      if (roles.target) {
        const t = textOf(roles.target.get(i));
        if (t) target = t;
      }
      const sys = systemOf(target);
      const onOwn = target === ownSystem;

      let keep = true;
      let s = 0;
      let e = 0;
      if (keyed) {
        s = sys.ranges.length;
        e = s + 1;
      } else if (hasRanges) {
        const a = positionOf((roles.position ?? roles.start)!.get(i));
        const b = roles.position ? a : positionOf(roles.end!.get(i));
        if (a === undefined || b === undefined) {
          causes.add('seq.source.records.missing', 'A position is missing', i);
          keep = false;
        } else if (Number.isNaN(a) || Number.isNaN(b)) {
          causes.add('seq.source.records.value', 'A position is not an integer', i);
          keep = false;
        } else if (b < a) {
          causes.add('seq.source.records.order', 'A row ends before it starts', i);
          keep = false;
        } else {
          if (a === 0) zeroStart = true;
          s = a - base;
          e = roles.position || inclusive ? b - base + 1 : b - base;
          if (e > MaxInt32) {
            causes.add('seq.source.records.value', 'An end is past int32', i);
            keep = false;
          } else if (s < 0) {
            causes.add('seq.source.records.range', BeforeFirst, i);
            keep = false;
          } else if (onOwn && L !== undefined && (s >= L || e > L)) {
            causes.add('seq.source.records.range', OutsideLetters, i);
            if (e === L + 1) pastByOne++;
            else pastMore++;
            if (s >= L) keep = false;
            else {
              e = L;
              sys.clipped = true;
            }
          }
        }
      }
      let value = Number.NaN;
      if (roles.value) {
        const cell = roles.value.get(i);
        const v = numberOf(cell);
        if (v !== undefined && Number.isNaN(v) && !isNaNCell(cell)) {
          causes.add('seq.source.records.value', 'A value is not a number', i);
          keep = false;
        } else value = v ?? Number.NaN;
      }
      // a row's value without `position` is at its own index: a row left out is NaN there, never a shift
      if (hasValues && !roles.position) sys.values.set(sys.valueRows++, keep ? value : Number.NaN);
      if (!keep) continue;
      if (hasRanges || keyed) {
        sys.ranges.push(i);
        sys.start.push(s);
        sys.end.push(e);
      }
      if (hasValues && roles.position) {
        if (sys.values.has(s)) causes.add('seq.source.records.cell', 'A position is named twice: the last kept', i);
        sys.values.set(s, value);
      }
    }
    if (systems.size === 0) systemOf(ownSystem);
    // the hints: the other base, when the rows say so
    if (zeroStart && base === 1)
      causes.hint('seq.source.records.range', BeforeFirst, 'The file counts from 0: try `base: 0`');
    if (pastByOne > 0 && pastMore === 0)
      causes.hint(
        'seq.source.records.range',
        OutsideLetters,
        `Every row would fit one position lower: try \`base: ${1 - base}\`${inclusive ? ' or `ends: exclusive`' : ''}`,
      );

    // the tables, a system each
    const systemInfo: SystemInfo[] = [];
    const shapes: ShapeInfo[] = [];
    const pairs: SeqPair[] = [];
    const single = systems.size === 1;
    // a typed column is kept when it is in the protocol's conventions already
    const conventions = base === 0 && !inclusive && !roles.position;
    for (const sys of systems.values()) {
      const tables: { -readonly [K in keyof Tables]: Tables[K] } = {};
      const isOwn = sys.name === ownSystem;
      if (isOwn && letters !== undefined)
        tables.letters = { form: 'letters', start: 0, letters } satisfies LettersTable;
      if (hasRanges || keyed)
        tables.ranges = rangesTable(sys, roles, fields, records, causes, single && conventions && !sys.clipped);
      if (hasValues)
        tables.values = valuesTable(sys, records.count, single && !hasRanges ? roles.value : undefined, options.limits);
      // a keyed system is as long as its records
      const length = keyed ? sys.ranges.length : isOwn ? L : undefined;
      systemInfo.push(length === undefined ? { name: sys.name } : { name: sys.name, length });
      const ranges = tables.ranges as RangesList | undefined;
      const fieldInfo: FieldInfo[] | undefined = ranges?.fields?.map((f) => ({
        name: f.name,
        type: columnType(f.column)!,
      }));
      for (const shape of ['letters', 'ranges', 'values'] as const)
        if (tables[shape])
          shapes.push(
            shape === 'ranges' && fieldInfo
              ? { system: sys.name, shape, fields: fieldInfo }
              : { system: sys.name, shape },
          );
      pairs.push({
        // a system of unknown length: the pair holds every position a column can, the target's length stands
        request: {
          system: sys.name,
          start: 0,
          end: length ?? DefaultMaxColumnLength,
          level: 1,
          shapes: shapesOf(tables),
        },
        answer: { ok: true, tables },
      });
    }
    issues.push(...causes.issues());
    return { data: { coverage: { systems: systemInfo, shapes, levels: [1] }, pairs }, issues };
  } catch (e) {
    if (e instanceof Refused) return { issues: [...issues, e.issue] };
    throw e;
  }
}

function shapesOf(tables: Tables): Shape[] {
  return (['letters', 'ranges', 'values'] as const).filter((s) => tables[s] !== undefined);
}

function rangesTable(
  sys: SystemRows,
  roles: Partial<Record<Role, Cells>>,
  fields: readonly Field[],
  records: Records,
  causes: Causes,
  typed: boolean,
): RangesList {
  const n = sys.ranges.length;
  // sorted by start once, each row keeping its place in the file (`row`)
  const order = Int32Array.from({ length: n }, (_, k) => k);
  let sorted = true;
  for (let k = 1; k < n && sorted; k++) if (sys.start[k] < sys.start[k - 1]) sorted = false;
  if (!sorted) order.sort((a, b) => sys.start[a] - sys.start[b] || a - b);
  const rows = Int32Array.from(order, (k) => sys.ranges[k]);
  const identity = sorted && isIdentity(rows, records.count);

  // a typed column in the protocol's conventions with every row in order: kept as it is
  const keep = (c: Cells | undefined) => (typed && identity && c?.column instanceof Int32Array ? c.column : undefined);
  const start = keep(roles.start) ?? Int32Array.from(order, (k) => sys.start[k]);
  const end = keep(roles.end) ?? Int32Array.from(order, (k) => sys.end[k]);

  let strand: Int8Array | undefined;
  const sc = roles.strand;
  if (sc) {
    if (identity && sc.column instanceof Int8Array && validStrands(sc.column)) strand = sc.column;
    else {
      strand = new Int8Array(n);
      for (let k = 0; k < n; k++) {
        const v = strandOf(sc.get(rows[k]));
        if (v === undefined)
          causes.add('seq.source.records.cell', 'A strand is not +, -, ., ?, 1, -1 or 0: read as 0', rows[k]);
        strand[k] = v ?? 0;
      }
    }
  }
  let id: CategoryColumn | undefined;
  const ic = roles.id;
  if (ic) {
    if (identity && ic.column instanceof CategoryColumn) id = ic.column;
    else
      id = categoryColumn(
        Array.from(rows, (r) => {
          const t = textOf(ic.get(r));
          return t === false ? null : t;
        }),
      );
  }
  const out: RangesField[] = [];
  const vc = roles.value;
  if (vc) {
    if (identity && vc.column instanceof Float32Array) out.push({ name: 'value', column: vc.column });
    else out.push({ name: 'value', column: Float32Array.from(rows, (r) => numberOf(vc.get(r)) ?? Number.NaN) });
  }
  for (const f of fields) out.push({ name: f.name, column: fieldColumn(f, rows, identity, causes, records) });
  return {
    form: 'list',
    start,
    end,
    ...(strand && { strand }),
    ...(id && { id }),
    row: rows,
    ...(out.length && { fields: out }),
  };
}

function validStrands(c: Int8Array): boolean {
  for (let i = 0; i < c.length; i++) if (c[i] < -1 || c[i] > 1) return false;
  return true;
}

/**
 * Values a position from the first one named: sparse by position (NaN between), or a row's each from the first
 * position. A typed float32 column of every row is kept as it is. More positions than `maxRows` are refused.
 */
function valuesTable(sys: SystemRows, count: number, typed: Cells | undefined, limits: OpenLimits): ValuesPositions {
  const c = typed?.column;
  if (c instanceof Float32Array && sys.valueRows === count) return { form: 'positions', start: 0, value: c };
  if (sys.values.size === 0) return { form: 'positions', start: 0, value: new Float32Array(0) };
  let first = Number.POSITIVE_INFINITY;
  let last = 0;
  for (const p of sys.values.keys()) {
    if (p < first) first = p;
    if (p > last) last = p;
  }
  if (last - first + 1 > limits.maxRows)
    tooLarge(
      `The values span ${last - first + 1} positions, more than ${limits.maxRows}`,
      'Use a format read by region',
    );
  const value = new Float32Array(last - first + 1).fill(Number.NaN);
  for (const [p, v] of sys.values) value[p - first] = v;
  return { form: 'positions', start: first, value };
}

// ---- the matrix read ----

/** A matrix: `rows × cols` values, row-major (the value at row i, column j is at `i × cols + j`), NaN missing. */
export interface Matrix {
  readonly rows: number;
  readonly cols: number;
  readonly values: Float32Array;
}

/** The size a matrix is expected to have: its axes' total (a grid's `|x| × |y|`). */
export interface MatrixSize {
  readonly rows: number;
  readonly cols: number;
}

/**
 * An array of arrays (JSON or MessagePack, found by the transform's `records` path) as a grid's values, row-major;
 * `null` and cells that are no number are NaN. A typed flat array (a MessagePack `Float32Array`) is kept as it is when it
 * is a `Float32Array`, its dimensions the expected size's, else a square. Its size is checked against `expected`
 * (`seq.source.records.size`, the source unavailable). Never throws on its input.
 */
export function readMatrix(
  input: RecordsInput,
  transform: RecordsTransform,
  options: { readonly expected?: MatrixSize; readonly limits: OpenLimits },
): { matrix?: Matrix; issues: Issue[] } {
  const found = findMatrix(input, transform);
  return found.issues.length ? { issues: found.issues } : matrixOf(found.at, transform.records, options.expected);
}

/** The value a matrix source's `records` path reaches, parsed once (a source keeps it, not the text it came from). */
export function findMatrix(input: RecordsInput, transform: RecordsTransform): { at?: unknown; issues: Issue[] } {
  try {
    const root = parseDocument(input);
    const path = transform.records;
    const at = path === undefined || path === '' ? root : walk(root, segments(path));
    if (at === undefined) refuse('seq.source.records.field', `The matrix’s path “${path}” reaches nothing`);
    return { at, issues: [] };
  } catch (e) {
    if (e instanceof Refused) return { issues: [e.issue] };
    throw e;
  }
}

/** Whether a value is a flat matrix: a typed array, or an array of numbers and nulls. */
function isFlat(at: unknown): boolean {
  if (ArrayBuffer.isView(at)) return true;
  if (!Array.isArray(at) || at.length === 0) return false;
  for (const v of at) if (v !== null && typeof v !== 'number') return false;
  return true;
}

/** A parsed value as a matrix of the expected size (the size checked on every call: a grid asks with its axes'). */
export function matrixOf(
  at: unknown,
  path: string | undefined,
  expected?: MatrixSize,
): { matrix?: Matrix; issues: Issue[] } {
  const causes = new Causes();
  try {
    const check = (rows: number, cols: number): void => {
      if (rows * cols > DefaultMaxBlockValues)
        tooLarge(
          `The matrix holds ${rows * cols} values, more than ${DefaultMaxBlockValues}`,
          'Use a format read by region',
        );
      if (expected && (rows !== expected.rows || cols !== expected.cols))
        refuse(
          'seq.source.records.size',
          `The matrix is ${rows} × ${cols}, not ${expected.rows} × ${expected.cols}`,
          RecordsIssues.codes['seq.source.records.size'].hint,
        );
    };
    // a flat array: its dimensions the expected size's, else a square
    if (isFlat(at)) {
      const flat = at as ArrayLike<unknown>;
      const side = Math.round(Math.sqrt(flat.length));
      const rows = expected?.rows ?? side;
      const cols = expected?.cols ?? side;
      if (rows * cols !== flat.length)
        refuse(
          'seq.source.records.size',
          `The flat matrix has ${flat.length} values, not ${rows} × ${cols}`,
          RecordsIssues.codes['seq.source.records.size'].hint,
        );
      check(rows, cols);
      if (flat instanceof Float32Array) return { matrix: { rows, cols, values: flat }, issues: [] };
      const values = new Float32Array(flat.length);
      for (let i = 0; i < flat.length; i++) values[i] = cellNumber(flat[i], i, causes);
      return { matrix: { rows, cols, values }, issues: causes.issues() };
    }
    if (!Array.isArray(at))
      refuse(
        'seq.source.records.malformed',
        path ? `The path “${path}” reaches no array of arrays` : 'The document holds no array of arrays',
      );
    const rows = at.length;
    const first = rows ? (at[0] as ArrayLike<unknown>) : undefined;
    const cols = first && typeof first === 'object' ? first.length : 0;
    check(rows, cols);
    const values = new Float32Array(rows * cols);
    for (let i = 0; i < rows; i++) {
      const row = at[i] as ArrayLike<unknown> | null | undefined;
      if (typeof row !== 'object' || row === null || row.length !== cols)
        refuse(
          'seq.source.records.size',
          `Row ${i} of the matrix has ${row && typeof row === 'object' ? row.length : 'no'} values, not ${cols}`,
          RecordsIssues.codes['seq.source.records.size'].hint,
        );
      for (let j = 0; j < cols; j++) values[i * cols + j] = cellNumber(row[j], i, causes);
    }
    return { matrix: { rows, cols, values }, issues: causes.issues() };
  } catch (e) {
    if (e instanceof Refused) return { issues: [e.issue] };
    throw e;
  }
}

/** A matrix cell: a number, NaN for `null`, missing and anything else (said once). */
function cellNumber(cell: unknown, row: number, causes: Causes): number {
  if (typeof cell === 'number') return cell;
  if (cell !== null && cell !== undefined)
    causes.add('seq.source.records.cell', 'A matrix cell is not a number: missing', `row ${row}`);
  return Number.NaN;
}

// ---- the runtime ----

/** The transform an entry carries (the source node's params, as the resolver copies them). */
function transformOf(s: SourceEntry): RecordsTransform {
  return {
    fields: s.fields,
    records: s.records,
    letters: s.letters,
    system: s.system,
    base: s.base,
    ends: s.ends,
    types: s.types,
  };
}

/** An opened records source: the `seq` join over its tables, with what its reading said (row-level issues). */
export interface RecordsSource extends OpenedSource {
  readonly issues: readonly Issue[];
  /** A source of `shape: matrix`: its array of arrays as a grid's values, checked against the axes' total. */
  readonly matrix?: (expected?: MatrixSize) => { matrix?: Matrix; issues: Issue[] };
}

const rank = (i: Issue) => (i.severity === 'error' ? 0 : i.severity === 'warning' ? 1 : 2);

function isFormat(f: string): f is RecordsFormat {
  return (RecordsFormats as readonly string[]).includes(f);
}

async function open(source: SourceEntry, ctx: OpenContext): Promise<RecordsSource> {
  const fail = (issue: Issue): never => {
    throw new SourceOpenError([{ ...issue, address: issue.address ?? source.address }]);
  };
  const format = source.format;
  if (!isFormat(format))
    return fail(issueOf(SourceIssues, 'seq.source.no_data', `The records runtime does not read the format ${format}`));
  let input: RecordsInput;
  if (source.value !== undefined) {
    input =
      typeof source.value === 'string' && format !== 'msgpack'
        ? { format, text: source.value }
        : { value: source.value };
  } else if (source.uri !== undefined) {
    let bytes: Uint8Array;
    try {
      bytes = await ctx.assets.readBytes(source.uri, ctx.signal);
    } catch (e) {
      if (isAbortError(e) || ctx.signal.aborted) throw e;
      return fail(issueOf(SourceIssues, 'seq.source.unreadable', `The file ${source.uri} could not be read`));
    }
    // bounded before it is parsed: a parse makes several times a file's bytes
    if (bytes.byteLength > ctx.limits.maxBytes)
      return fail(
        issueOf(
          SourceIssues,
          'seq.source.too_large',
          `The file has ${bytes.byteLength} bytes, more than ${ctx.limits.maxBytes}`,
          {
            hint: 'Use a format read by region',
          },
        ),
      );
    input = format === 'msgpack' ? { format, bytes } : { format, text: new TextDecoder('utf-8').decode(bytes) };
  } else return fail(issueOf(SourceIssues, 'seq.source.no_data', 'A records source needs a value or a uri'));
  if (source.shape === 'matrix') {
    // a matrix serves no region: a grid reads it whole, once its axes' total is known
    if (format === 'gff3' || format === 'tsv' || format === 'csv')
      return fail(
        issueOf(SourceIssues, 'seq.source.no_data', `A matrix is read from JSON or MessagePack, not ${format}`),
      );
    // parsed once, here: the source keeps the array it reaches, not the file's text or bytes
    const transform = transformOf(source);
    const found = findMatrix(input, transform);
    const at = found.at;
    const empty = openSeq({ coverage: { systems: [], shapes: [], levels: [1] }, pairs: [] });
    const placed = (issues: Issue[]) => issues.map((i) => ({ ...i, address: i.address ?? source.address }));
    return {
      ...empty,
      issues: [],
      matrix: (expected) => {
        const r = found.issues.length ? { issues: found.issues } : matrixOf(at, transform.records, expected);
        return { ...r, issues: placed(r.issues) };
      },
    };
  }
  const id = source.address.slice(source.address.lastIndexOf('/') + 1);
  const read = readRecords(input, transformOf(source), { system: id, limits: ctx.limits });
  const placed = read.issues.map((i) => ({ ...i, address: i.address ?? source.address }));
  // the refusal first: the store's report of an unopened source is its first issue
  if (!read.data) throw new SourceOpenError([...placed].sort((a, b) => rank(a) - rank(b)));
  return { ...openSeq(read.data), issues: placed };
}

/** The runtime of files of records: `tsv`, `csv`, `json`, `msgpack`. */
export const RecordsRuntime: FormatRuntime = { name: 'records', formats: RecordsFormats, open };
