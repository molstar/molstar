/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node, Value } from '@molstar/core/repr/tree';

// A layer's data filter read row by row (_spec/README.md "Themes", Filters): the query's form after its annotation,
// as the query sublanguage evaluates it on items (`item` by `:id` or a field, `except`, `intersect`, `merge`), compiled once a
// prepare to a predicate over a row's key and fields. A form it cannot read row by row filters nothing and is said.

/** A row as a filter reads it: its element key and a field's value. */
export interface FilterRow {
  readonly key: string | number;
  field(name: string): string | number | undefined;
}

export type RowFilter = (row: FilterRow) => boolean;

const isNode = (v: unknown): v is Node => typeof v === 'object' && v !== null && !Array.isArray(v) && 'kind' in v;

/** Whether a row's value fits a param: the same atom, or an atom a list param holds (the query's `item` rule). */
function fits(v: Value, x: string | number | undefined): boolean {
  if (x === undefined) return false;
  if (Array.isArray(v)) return v.some((y) => Object.is(y, x));
  return Object.is(v, x);
}

/** The predicate of a filter form; undefined for a form not read row by row (a location, a range). */
export function compileRowFilter(form: Node): RowFilter | undefined {
  const params = form.params ?? {};
  switch (form.kind) {
    case 'item': {
      const keys = Object.keys(params);
      if (!keys.length) return undefined;
      return (row) => {
        for (const k of keys) if (!fits(params[k], k === 'id' ? row.key : row.field(k))) return false;
        return true;
      };
    }
    case 'except': {
      // one operand: the rows it rejects; several: the first's rows the others reject
      const of = params.of;
      const ops = Array.isArray(of) ? of : [of];
      const parts: RowFilter[] = [];
      for (const o of ops) {
        const f = isNode(o) ? compileRowFilter(o) : undefined;
        if (!f) return undefined;
        parts.push(f);
      }
      if (!parts.length) return undefined;
      if (parts.length === 1) return (row) => !parts[0](row);
      const [first, ...rest] = parts;
      return (row) => first(row) && !rest.some((f) => f(row));
    }
    case 'intersect':
    case 'merge': {
      const of = params.of;
      if (!Array.isArray(of)) return undefined;
      const parts: RowFilter[] = [];
      for (const o of of) {
        const f = isNode(o) ? compileRowFilter(o) : undefined;
        if (!f) return undefined;
        parts.push(f);
      }
      return form.kind === 'intersect' ? (row) => parts.every((f) => f(row)) : (row) => parts.some((f) => f(row));
    }
    default:
      return undefined;
  }
}
