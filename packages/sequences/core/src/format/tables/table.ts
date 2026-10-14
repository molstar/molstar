/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `Table` (seqspec:format#positions) and the annotation columns per geometry (seqspec:format#annotation-tables): rows
// or columns, typed, dictionary-coded and JSON columns; normalized to columns with codes decoded
// (seqspec:format#what-a-change-affects). Which columns a geometry reads is the model's `ColumnRoles`
// (seqspec:model#annotations), one table for inline validation, file typing and the model.
//
// The check costs O(cells) with nothing allocated per cell, O(distinct values) for refs, and paths built on failure
// only; a repeated failure is one issue with its count. A table's columns are read once per table object (`tableOf`),
// on the immutability `validate`'s memo assumes.

import { ColumnRole, ColumnRoles, columnNumbers, type Geometry } from "@molstar/seq-core/model";
import {
  type CheckContext,
  type FieldType,
  isPlainObject,
  isRefSyntax,
  ptr,
  report,
  show,
  typeError,
} from "../schema/field-schema.js";
import { CellKind, CellScan, type Col, ColumnForm, cellAt, cellKind, cellOk, scanCells } from "./cells.js";

export type TypedColumn =
  | Int8Array
  | Uint8Array
  | Int16Array
  | Uint16Array
  | Int32Array
  | Uint32Array
  | Float32Array
  | Float64Array;
/** A column (seqspec:format#values): numbers may be a typed array, and so may a dictionary column's codes. */
export type TableColumn =
  | readonly unknown[]
  | TypedColumn
  | { dict: readonly string[]; codes: readonly (number | null)[] | TypedColumn }
  | { json: readonly unknown[] };
export type Table = Record<string, unknown>[] | Record<string, TableColumn>;

/** Per geometry, the columns it reads (`ColumnRoles`' names); every other column is a field. */
export const GeometryColumns: Readonly<Record<Geometry, readonly string[]>> = Object.fromEntries(
  Object.entries(ColumnRoles).map(([g, roles]) => [g, Object.keys(roles)]),
) as unknown as Record<Geometry, readonly string[]>;
/** Per geometry, sets of which a table with rows has at least one column. */
const REQUIRED: Readonly<Record<Geometry, readonly (readonly string[])[]>> = {
  ranges: [["start", "beg", "boundary"]],
  values: [["value"]],
  pairs: [["start", "beg"]],
  matrix: [["index", "pos"], ["value"]],
  grid: [["value"]],
  rows: [["row", "state", "node"]],
};
/** Target and positional columns: read only where a geometry lists them, `svs.columns.unknown` elsewhere. */
const PLACING = new Set(
  Object.values(ColumnRoles).flatMap((r) => Object.keys(r).filter((k) => r[k] !== ColumnRole.descriptive)),
);
const TARGET_KINDS: readonly string[] = ["sequence", "ligand", "unit", "alignment", "row"];
const EXPECTED = [
  "an index (int >= 0)",
  'Label (315, "27A", "111.1", or a name such as "N13")',
  "a ref",
  "float | null",
  "a list of float | null, one per state",
  '"ab" | "ba" | null',
  "str",
  "a list of float | null (or a Float32Array, a Float64Array), one per pair of positions, a-major",
];
const rows = (n: number) => (n > 1 ? ` (${n} rows)` : "");

/** An array, or a typed array (codes may be one too). */
export const isList = (v: unknown): v is ArrayLike<unknown> => Array.isArray(v) || ArrayBuffer.isView(v);

/** A table's columns, read once per table object; what is malformed is kept for the check to report. */
interface TableInfo {
  readonly cols: Map<string, Col>;
  /** The rows form's length, or the first column's. */
  readonly length: number;
  /** Rows form: rows that are not objects, and the first; columns form: columns of no known shape. */
  readonly badRows: number;
  readonly firstBadRow: number;
  readonly badShapes: readonly string[];
  /** Columns with a value in some row, by name, as asked. */
  readonly present: Map<string, boolean>;
  sides?: LabelSides;
}
const INFO = new WeakMap<object, TableInfo>();

function columnOf(name: string, col: unknown): Col | null {
  if (Array.isArray(col)) return { name, form: ColumnForm.array, data: col, dict: null, length: col.length };
  if (ArrayBuffer.isView(col) && !(col instanceof DataView)) {
    const data = col as unknown as ArrayLike<number>;
    return { name, form: ColumnForm.typed, data, dict: null, length: data.length };
  }
  if (!isPlainObject(col)) return null;
  if (Array.isArray(col.json))
    return { name, form: ColumnForm.json, data: col.json, dict: null, length: col.json.length };
  if (Array.isArray(col.dict) && isList(col.codes)) {
    return { name, form: ColumnForm.dict, data: col.codes, dict: col.dict, length: col.codes.length };
  }
  return null;
}

/** A table's columns; null when it is neither form. */
function tableOf(v: unknown): TableInfo | null {
  if (!Array.isArray(v) && !isPlainObject(v)) return null;
  let t = INFO.get(v);
  if (t) return t;
  const cols = new Map<string, Col>();
  let badRows = 0;
  let firstBadRow = -1;
  const badShapes: string[] = [];
  let length = -1;
  if (Array.isArray(v)) {
    for (let i = 0; i < v.length; i++) {
      const row = v[i];
      if (!isPlainObject(row)) {
        if (badRows++ === 0) firstBadRow = i;
        continue;
      }
      for (const k in row) {
        if (!cols.has(k)) cols.set(k, { name: k, form: ColumnForm.rows, data: v, dict: null, length: v.length });
      }
    }
    length = v.length;
  } else {
    for (const [k, raw] of Object.entries(v)) {
      const c = columnOf(k, raw);
      if (!c) badShapes.push(k);
      else {
        cols.set(k, c);
        if (length < 0) length = c.length;
      }
    }
  }
  t = { cols, length: Math.max(0, length), badRows, firstBadRow, badShapes, present: new Map() };
  INFO.set(v, t);
  return t;
}

/** Whether column `name` holds a value in some row. */
function present(t: TableInfo, name: string): boolean {
  let p = t.present.get(name);
  if (p !== undefined) return p;
  const c = t.cols.get(name);
  p = false;
  for (let i = 0; c && !p && i < c.length; i++) p = cellAt(c, i) !== null;
  t.present.set(name, p);
  return p;
}
const anyPresent = (t: TableInfo, names: readonly string[]) => names.some((k) => present(t, k));

/** Whether an inline table has a column `name` (a split's `by`, seqspec:format#layer-types). */
export const tableHas = (v: unknown, name: string): boolean => tableOf(v)?.cols.has(name) ?? false;

export interface LabelSides {
  readonly a: boolean;
  readonly b: boolean;
  readonly numbering: boolean;
  readonly bNumbering: boolean;
}
const NO_SIDES: LabelSides = { a: false, b: false, numbering: false, bNumbering: false };

/** Which sides of a table are written in labels, and whether it names numberings per row. */
export function labelSides(v: unknown, _ctx?: CheckContext): LabelSides {
  const t = tableOf(v);
  if (!t) return NO_SIDES;
  t.sides ??= {
    a: anyPresent(t, ["beg", "end", "pos"]),
    b: anyPresent(t, ["b_beg", "b_end"]),
    numbering: present(t, "numbering"),
    bNumbering: present(t, "b_numbering"),
  };
  return t.sides;
}

/** A cell's path, built on failure only. */
function cellPath(path: string, c: Col, i: number): string {
  if (c.form === ColumnForm.rows) return ptr(ptr(path, i), c.name);
  const at = ptr(path, c.name);
  if (c.form === ColumnForm.dict) return ptr(ptr(at, "codes"), i);
  return ptr(c.form === ColumnForm.json ? ptr(at, "json") : at, i);
}

const scan = new CellScan();

/** One column of a kind the geometry reads: each failure kind reported once with its count, one `RefUse` per value. */
function checkColumn(c: Col, kind: number, path: string, ctx: CheckContext) {
  scan.reset(c);
  scanCells(c, kind, 0, c.length, scan);
  const at = ptr(path, c.name);
  const code = kind === CellKind.label ? "svs.position.syntax" : "svs.columns.type";
  let firstPath = scan.firstBad < 0 ? "" : cellPath(path, c, scan.firstBad);
  let first = scan.firstBad < 0 ? null : cellAt(c, scan.firstBad);
  const refs: { ref: string; count: number; path: string }[] = [];
  if (c.form === ColumnForm.dict) {
    const dict = c.dict ?? [];
    if (scan.badCodes > 0) {
      const bad = c.data[scan.firstBadCode];
      const message = `expected a code in 0..${dict.length - 1} or null, got ${show(bad)}${rows(scan.badCodes)}`;
      report(ctx, "svs.columns.type", cellPath(path, c, scan.firstBadCode), message);
    }
    // Each entry the codes use, checked once for all its rows.
    for (let e = 0; e < dict.length; e++) {
      const n = scan.perEntry[e];
      const x = dict[e];
      if (n === 0 || x === null) continue;
      if (kind === CellKind.ref && typeof x === "string")
        refs.push({ ref: x, count: n, path: ptr(ptr(at, "dict"), e) });
      else if (!cellOk(kind, x)) {
        if (scan.bad === 0) {
          firstPath = ptr(ptr(at, "dict"), e);
          first = x;
        }
        scan.bad += n;
      }
    }
  } else {
    for (const [ref, slot] of scan.distinct)
      refs.push({ ref, count: scan.counts[slot], path: cellPath(path, c, scan.firsts[slot]) });
  }
  if (scan.bad > 0) report(ctx, code, firstPath, `expected ${EXPECTED[kind]}, got ${show(first)}${rows(scan.bad)}`);
  for (const r of refs) {
    if (!isRefSyntax(r.ref, true)) {
      const message = `"${r.ref}" is not a ref: [A-Za-z0-9][A-Za-z0-9_.-]*, then /suffixes${rows(r.count)}`;
      report(ctx, "svs.ref.syntax", r.path, message);
    } else ctx.refs.push({ ref: r.ref, kinds: TARGET_KINDS, path: r.path, from: ctx.from, count: r.count });
  }
}

export function table(geometry: Geometry): FieldType<Table> {
  const roles = ColumnRoles[geometry];
  return {
    name: `Table<${geometry}>`,
    check(v, path, ctx) {
      const t = tableOf(v);
      if (!t) return typeError(ctx, path, "a Table: row objects or column arrays", v);
      if (t.badRows > 0) {
        const message = `expected a row object, got ${show((v as unknown[])[t.firstBadRow])}${rows(t.badRows)}`;
        report(ctx, "svs.params.type", ptr(path, t.firstBadRow), message);
      }
      for (const k of t.badShapes) {
        const col = (v as Record<string, unknown>)[k];
        typeError(ctx, ptr(path, k), "a column: an array, {dict, codes} or {json}", col, "svs.columns.type");
      }
      const first = t.cols.values().next().value;
      for (const c of t.cols.values()) {
        const at = c.form === ColumnForm.rows ? path : ptr(path, c.name);
        if (first && c.form !== ColumnForm.rows && c.length !== first.length) {
          const message = `column "${c.name}" has ${c.length} values; the first column has ${first.length}`;
          report(ctx, "svs.columns.length", at, message);
        }
        if (c.dict && !c.dict.every((s) => typeof s === "string"))
          typeError(ctx, ptr(at, "dict"), "list<str>", c.dict, "svs.columns.type");
        if (roles[c.name] !== undefined) checkColumn(c, cellKind(c.name, geometry), path, ctx);
        else if (PLACING.has(c.name)) {
          const names = GeometryColumns[geometry].join(", ");
          report(ctx, "svs.columns.unknown", at, `"${c.name}" is not a ${geometry} column; ${geometry} reads ${names}`);
        }
      }
      // A table with no rows needs no columns.
      for (const set of t.length > 0 ? REQUIRED[geometry] : []) {
        if (!set.some((k) => t.cols.has(k))) {
          report(ctx, "svs.params.missing", path, `a ${geometry} table needs a column ${set.join(" or ")}`);
        }
      }
      const form = (idx: string[], lab: string[], side: string) => {
        if (anyPresent(t, idx) && anyPresent(t, lab)) {
          report(ctx, "svs.position.form", path, `${side}: use indices (${idx}) or labels (${lab}), not both`);
        }
      };
      form(["start", "stop", "boundary", "index"], ["beg", "end", "pos"], "positions");
      form(["b_start", "b_stop"], ["b_beg", "b_end"], "side b");
      // An insertion code completes an author number: it needs its side's labels.
      for (const [code, lab] of [
        ["ins_code", ["beg", "pos"]],
        ["b_ins_code", ["b_beg"]],
      ] as const) {
        if (!present(t, code) || anyPresent(t, lab)) continue;
        const idx = code === "ins_code" ? "start, index" : "b_start";
        const m = `${code} completes an author number in ${lab.join(" or ")} ("52" and "A" name 52A); without ${lab.join(" or ")} (indices: ${idx}) there is none: give positions as labels in an author numbering, or leave ${code} out`;
        report(ctx, "svs.position.form", path, m);
      }
    },
    normalize(v) {
      const t = tableOf(v);
      if (!t) return v;
      // A document value, as JSON.parse makes it: a `__proto__` column is an own key (CF-02), never the prototype.
      const out: Record<string, TableColumn> = {};
      const put = (k: string, v: TableColumn) =>
        Object.defineProperty(out, k, { value: v, enumerable: true, writable: true, configurable: true });
      for (const c of t.cols.values()) {
        if (c.form === ColumnForm.array) put(c.name, c.data as unknown[]);
        else if (c.form === ColumnForm.typed) put(c.name, columnNumbers(c.data as Float64Array, null));
        else {
          const values = Array.from({ length: c.length }, (_, i) => cellAt(c, i));
          // JSON only from a `{json}` column.
          put(c.name, c.form === ColumnForm.json ? { json: values } : values);
        }
      }
      return out;
    },
  };
}

/**
 * A file table's columns renamed by `fields`: every entry reads the columns as read, so names can be swapped, chained
 * or fanned out; null ignores a column; a column no entry names or reads keeps its name.
 */
export function renameFields<C>(
  columns: Readonly<Record<string, C>>,
  fields: Readonly<Record<string, string | null>>,
): Record<string, C> {
  const named = new Set<string>();
  for (const [k, from] of Object.entries(fields)) {
    named.add(k);
    if (from != null) named.add(from);
  }
  // Keyed by data: no prototype, so a `__proto__` field is a column (CF-02).
  const out: Record<string, C> = Object.create(null);
  for (const [k, c] of Object.entries(columns)) if (!named.has(k)) out[k] = c;
  for (const [k, from] of Object.entries(fields))
    if (from != null && Object.hasOwn(columns, from)) out[k] = columns[from];
  return out;
}

const NUMBER = /^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$/;
/** Columns some geometry reads as text: refs, labels and the descriptive strings (`value` aside). */
const TEXT = new Set(
  Object.values(ColumnRoles).flatMap((r) =>
    Object.keys(r).filter((k) => {
      const kind = cellKind(k, "ranges");
      return kind !== CellKind.index && kind !== CellKind.value;
    }),
  ),
);

/**
 * A file table's text cells typed after `fields`: a column some geometry reads as text (a ref, a label, an `id`, a
 * `category`) stays text, so `01` stays `01`; an index column of integers is an `Int32Array` (a `Float64Array`, NaN
 * missing, when one is missing); `value` and fields whose every cell is a number are a `Float64Array`. An index or
 * `value` column with a cell that is no number is a list, its numbers as numbers and the rest as written (the loader's
 * check names and leaves out the bad cells: CF-03); any other column stays text.
 */
export function typeColumns(
  columns: Readonly<Record<string, readonly (string | null)[]>>,
): Record<string, readonly (string | number | null)[] | Float64Array | Int32Array> {
  const out: Record<string, readonly (string | number | null)[] | Float64Array | Int32Array> = Object.create(null);
  for (const [k, cells] of Object.entries(columns)) {
    if (TEXT.has(k) || !cells.every((s) => s === null || NUMBER.test(s))) {
      // An index or value column with a cell that is no number: its numbers read as numbers, the rest left as
      // written, for the static check to name and the loader to leave out (CF-03).
      const numeric = cellKind(k, "ranges") === CellKind.index || k === "value";
      out[k] = !TEXT.has(k) && numeric ? cells.map((s) => (s !== null && NUMBER.test(s) ? Number(s) : s)) : cells;
      continue;
    }
    const nums = Float64Array.from(cells, (s) => (s === null ? Number.NaN : Number(s)));
    const int32 = cellKind(k, "ranges") === CellKind.index && nums.every((x) => (x | 0) === x);
    out[k] = int32 ? Int32Array.from(nums) : nums;
  }
  return out;
}
