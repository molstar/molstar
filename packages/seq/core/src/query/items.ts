/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The schema selection's items (_spec/README.md "Queries", Language; addressing.md, Schema selections): `elements`
// takes a table of items, each a record of properties: one item as the form's params, columns as lists (row i across
// them one item, a scalar the same in every row; several lists zip, never a product), or records as
// `:items [(fields …) …]` (`fields` a record, never evaluated); in JSON and the builder an object (columns where its
// values are arrays) or an array of objects, as a records file. Normalized to columns (typed by each property, long
// names read as the short one, a field absent in every row dropped, a row given twice once) and printed in the shorter
// layout (one item as params, else columns or records, ties to columns), so `parse(print(a))` equals `a` as a value.

import { type Atom, isNode, type Node, printAddress, type Value } from '@molstar/core/repr/tree';
import { coerce } from './expressions.js';
import { fail } from './form-kit.js';
import { findProperty, levelRank, type PropertyDef, type PropertyTable } from './properties.js';

/** The range fields: not properties, the sugar of `in-range` over `seq` and over `auth-seq` in the numbering's order. */
export const RangeFields: Readonly<Record<string, string>> = {
  beg: 'beg',
  end: 'end',
  'auth-beg': 'auth-beg',
  'auth-end': 'auth-end',
  beg_label_seq_id: 'beg',
  end_label_seq_id: 'end',
  beg_auth_seq_id: 'auth-beg',
  end_auth_seq_id: 'auth-end',
};

/** A table of items: its rows and, by each field's short name (sorted), a column (`null`: no test on that row). */
export interface Items {
  readonly rows: number;
  readonly fields: readonly string[];
  readonly columns: Readonly<Record<string, readonly Atom[]>>;
}

/** What JSON and the builder give: an object (columns where its values are arrays), or an array of objects. */
export type ItemsInput = Readonly<Record<string, Atom | readonly Atom[]>> | readonly Readonly<Record<string, Atom>>[];

/** A field's short name and its property (none for a range field). */
function fieldOf(t: PropertyTable, name: string): { readonly name: string; readonly def?: PropertyDef } {
  const range = RangeFields[name];
  if (range) return { name: range };
  const f = findProperty(t, name);
  if ('def' in f) return { name: f.def.name, def: f.def };
  if ('missing' in f)
    return fail('seq.query.extension.missing', `No extension provides the property ${name}`, {
      hint: `namespace: ${f.missing}`,
    });
  return fail('seq.query.form.param', `A form of kind elements names no property ${name}`, { hint: f.unknown });
}

/** A value read as its field's type: a property's, an integer for `beg` and `end`, a label for the auth ends. */
function typed(field: { name: string; def?: PropertyDef }, v: Value, row: number): Atom {
  if (v === null) return null;
  if (isNode(v) || Array.isArray(v) || typeof v === 'object')
    return fail('seq.query.form.param', `Row ${row + 1}: :${field.name} takes an atom`);
  const type = field.def?.type ?? (field.name === 'beg' || field.name === 'end' ? 'integer' : 'word');
  const c = coerce(v as Atom, type);
  if (c === undefined) return fail('seq.query.form.type', `Row ${row + 1}: :${field.name} takes a ${type}`);
  // an auth end is a label: an integer read back as one (`100`), a word as itself (`100A`)
  if (field.name === 'auth-beg' || field.name === 'auth-end') return typeof v === 'number' ? v : (c as string);
  return c;
}

/** Normalizes rows of named values: short names, typed columns, rows checked, a duplicate row and an empty field dropped. */
function normalize(t: PropertyTable, rows: readonly Readonly<Record<string, Value>>[]): Items {
  const names = new Map<string, { name: string; def?: PropertyDef }>();
  const cells: Record<string, Atom>[] = [];
  rows.forEach((row, i) => {
    const out: Record<string, Atom> = {};
    for (const k of Object.keys(row)) {
      const f = fieldOf(t, k);
      if (Object.hasOwn(out, f.name))
        fail('seq.query.form.param', `Row ${i + 1} gives ${f.name} twice (a short and a long name)`);
      names.set(f.name, f);
      out[f.name] = typed(f, row[k], i);
    }
    const set = (k: string) => out[k] !== undefined && out[k] !== null;
    if (set('seq') && (set('beg') || set('end')))
      fail('seq.query.form.param', `Row ${i + 1} gives seq with beg or end`, { hint: 'a position, or a range' });
    if (typeof out.beg === 'number' && typeof out.end === 'number' && out.beg > out.end)
      fail('seq.query.form.param', `Row ${i + 1}: beg is past end`);
    cells.push(out);
  });
  const fields = [...names.keys()].filter((f) => cells.some((c) => c[f] !== undefined && c[f] !== null)).sort();
  const seen = new Set<string>();
  const columns: Record<string, Atom[]> = Object.fromEntries(fields.map((f) => [f, []]));
  let count = 0;
  for (const c of cells) {
    const key = JSON.stringify(fields.map((f) => [typeof c[f], c[f] ?? null]));
    if (seen.has(key)) continue;
    seen.add(key);
    count++;
    for (const f of fields) columns[f].push(c[f] ?? null);
  }
  for (const f of fields) Object.freeze(columns[f]);
  return Object.freeze({ rows: count, fields: Object.freeze(fields), columns: Object.freeze(columns) });
}

/** Columns as rows: lists of one length zipped, scalars repeated; unequal lengths are a form error. */
function zip(columns: Readonly<Record<string, Value>>): Record<string, Value>[] {
  let n: number | undefined;
  for (const [k, v] of Object.entries(columns))
    if (Array.isArray(v)) {
      if (n !== undefined && v.length !== n)
        fail('seq.query.form.param', `The columns have unequal lengths (${k}: ${v.length}, before it ${n})`);
      n = v.length;
    }
  const rows: Record<string, Value>[] = [];
  for (let i = 0; i < (n ?? 1); i++) {
    const row: Record<string, Value> = {};
    for (const [k, v] of Object.entries(columns)) row[k] = Array.isArray(v) ? (v[i] as Value) : v;
    rows.push(row);
  }
  return rows;
}

/** An `elements` form's items: its params, or its `:items` records (alone). */
export function readItems(t: PropertyTable, form: Node): Items {
  const params = form.params ?? {};
  if (!Object.hasOwn(params, 'items')) return normalize(t, zip(params));
  if (Object.keys(params).length > 1)
    fail('seq.query.form.param', 'A form of kind elements takes :items alone', {
      hint: '(elements :items [(fields …) …])',
    });
  const list = params.items;
  if (!Array.isArray(list)) return fail('seq.query.form.param', ':items takes a list of (fields …) records');
  return normalize(
    t,
    list.map((r: Value, i) => {
      if (!isNode(r) || r.kind !== 'fields' || r.args?.length || r.children?.length)
        fail('seq.query.form.param', `Row ${i + 1} is not a (fields …) record`);
      return (r as Node).params ?? {};
    }),
  );
}

/** Items from JSON or the builder's object or records. */
export function itemsOf(t: PropertyTable, input: ItemsInput): Items {
  if (Array.isArray(input)) return normalize(t, input as readonly Readonly<Record<string, Value>>[]);
  return normalize(t, zip(input as Readonly<Record<string, Value>>));
}

export function itemsEqual(a: Items, b: Items): boolean {
  if (a.rows !== b.rows || a.fields.length !== b.fields.length) return false;
  for (let i = 0; i < a.fields.length; i++) {
    const f = a.fields[i];
    if (f !== b.fields[i]) return false;
    const x = a.columns[f];
    const y = b.columns[f];
    for (let r = 0; r < a.rows; r++) if (!Object.is(x[r], y[r])) return false;
  }
  return true;
}

/**
 * The `elements` form of items in the shorter layout: one item as params; else columns (a field the same in every
 * row a scalar) or `:items` records, whichever prints shorter, ties to columns.
 */
export function elementsForm(items: Items): Node {
  const { rows, fields, columns } = items;
  if (rows === 1) {
    const params: Record<string, Value> = {};
    for (const f of fields) params[f] = columns[f][0];
    return fields.length ? { kind: 'elements', params } : { kind: 'elements' };
  }
  const cols: Record<string, Value> = {};
  for (const f of fields) {
    const c = columns[f];
    cols[f] = c.every((v) => Object.is(v, c[0])) ? c[0] : [...c];
  }
  const records: Node[] = [];
  for (let r = 0; r < rows; r++) {
    const params: Record<string, Value> = {};
    for (const f of fields) if (columns[f][r] !== null) params[f] = columns[f][r];
    records.push({ kind: 'fields', params });
  }
  const byColumns: Node = { kind: 'elements', params: cols };
  const byRecords: Node = { kind: 'elements', params: { items: records } };
  return printAddress([byRecords]).length < printAddress([byColumns]).length ? byRecords : byColumns;
}

/** The finest level of items' fields (`levelRank`), -1 when none: the range fields are a residue's. */
export function itemsLevel(t: PropertyTable, items: Items): number {
  let level = -1;
  for (const f of items.fields) {
    const d = RangeFields[f] ? undefined : findProperty(t, f);
    level = Math.max(level, d && 'def' in d ? levelRank(d.def.level) : levelRank('residue'));
  }
  return level;
}

const prop = (name: string): Node => ({ kind: 'prop', args: [name] });

/**
 * The predicate items are sugar for, exactly (locations-and-queries.md, The language): a field `:f v` is
 * `(= (prop f) v)`, `:beg b :end e` is `(in-range (prop seq) b e)` (an end alone `>=` or `<=`), and `:auth-beg`,
 * `:auth-end` the same over `auth-seq` for integer labels; a row is `and` of its fields, items `or` of their rows.
 */
export function itemsPredicate(items: Items): Node {
  const rows: Node[] = [];
  for (let r = 0; r < items.rows; r++) {
    const tests: Node[] = [];
    const at = (f: string) => items.columns[f]?.[r] ?? null;
    for (const f of items.fields) {
      const v = at(f);
      if (v === null || RangeFields[f]) continue;
      tests.push({ kind: '=', args: [prop(f), v] });
    }
    for (const [lo, hi, p] of [
      ['beg', 'end', 'seq'],
      ['auth-beg', 'auth-end', 'auth-seq'],
    ] as const) {
      const b = at(lo);
      const e = at(hi);
      if (b !== null && e !== null) tests.push({ kind: 'in-range', args: [prop(p), b, e] });
      else if (b !== null) tests.push({ kind: '>=', args: [prop(p), b] });
      else if (e !== null) tests.push({ kind: '<=', args: [prop(p), e] });
    }
    rows.push({ kind: 'and', args: tests.length ? tests : [true] });
  }
  return { kind: 'or', args: rows.length ? rows : [false] };
}
