/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { isColumn, type Node, type Value } from '@molstar/core/repr/tree';

// The builder's theme forms (_spec/README.md, Language builder): `theme.field('label')`, `theme.scheme('zappo')`,
// rules and palettes as nodes of `seq.theme`, written into a layer's channel slots. A `field` form over inline data
// lists the column's values (`list`, true by default) when the layer is added: the order is written, not inferred.

type P = Record<string, Value>;

const form = (kind: string, params: P = {}, args?: Value[]): Node =>
  ({ kind, ...(Object.keys(params).length ? { params } : {}), ...(args ? { args } : {}) }) as Node;

const clean = (p: Record<string, unknown>): P => {
  const out: P = {};
  for (const [k, v] of Object.entries(p)) if (v !== undefined) out[k] = v as Value;
  return out;
};

/** Marks a field form whose values the builder lists (removed when it writes them). */
export const ListMarker = '__list';

export interface FieldOptions {
  readonly palette?: Node;
  /** A scale node's address (or handle), shared with the layers naming it. */
  readonly scale?: string | { readonly address: string };
  readonly domain?: readonly [number | null, number | null];
  readonly type?: 'linear' | 'log' | 'sqrt';
  readonly mid?: number;
  readonly range?: readonly [number, number];
  readonly data?: string | { readonly address: string };
  /** Over inline data, write the column's values in first-occurrence order (true). */
  readonly list?: boolean;
}

const addr = (v: string | { readonly address: string } | undefined) =>
  v === undefined ? undefined : typeof v === 'string' ? v : v.address;

/** Measures for clauses. */
export const measure = {
  field: (name: string): Node => form('field', { name }),
  value: (): Node => form('value'),
  letter: (): Node => form('letter'),
};

type Measure = Node | string;
const m = (x: Measure) => (typeof x === 'string' ? measure.field(x) : x);

export const theme = {
  uniform: (value: string | number): Node => form('uniform', { value }),
  field: (name: string, options: FieldOptions = {}): Node => {
    const { list, scale, data, ...rest } = options;
    const params = clean({ name, ...rest, scale: addr(scale), data: addr(data) });
    if (list !== false && !params.palette && !params.scale) (params as Record<string, Value>)[ListMarker] = true;
    return form('field', params);
  },
  scheme: (name: string, params: Record<string, unknown> = {}): Node => form('scheme', clean({ name, ...params })),
  letters: (params: { colors?: Record<string, string>; groups?: Node[]; default?: string }): Node =>
    form('letters', clean(params)),
  group: (letters: string, value: string, label?: string): Node => form('group', clean({ letters, value, label })),
  rules: (rules: readonly Node[], otherwise?: Value): Node => form('rules', clean({ else: otherwise }), [...rules]),
  rule: (when: readonly Node[], value: Value, label?: string): Node =>
    form('rule', clean({ when: [...when], value, label })),
  eq: (x: Measure, v: Value): Node => form('eq', {}, [m(x), v]),
  ne: (x: Measure, v: Value): Node => form('ne', {}, [m(x), v]),
  lt: (x: Measure, v: number): Node => form('lt', {}, [m(x), v]),
  le: (x: Measure, v: number): Node => form('le', {}, [m(x), v]),
  gt: (x: Measure, v: number): Node => form('gt', {}, [m(x), v]),
  ge: (x: Measure, v: number): Node => form('ge', {}, [m(x), v]),
  oneOf: (x: Measure, ...v: Value[]): Node => form('one-of', {}, [m(x), ...v]),
  in: (location: string | { readonly address: string }): Node => form('in', {}, [addr(location)!]),
  within: (location: string | { readonly address: string }): Node => form('within', {}, [addr(location)!]),
  categorical: (params: Record<string, unknown> = {}): Node => form('categorical', clean(params)),
  direct: (): Node => form('direct'),
  binned: (params: Record<string, unknown>): Node => form('binned', clean(params)),
  continuous: (params: Record<string, unknown> = {}): Node => form('continuous', clean(params)),
};

/** A column's values as text in first-occurrence order (a list, a category column, a typed array). */
function distinct(column: unknown): string[] | undefined {
  const values: unknown[] = isColumn(column as Value)
    ? Array.from({ length: (column as { length: number }).length }, (_, i) =>
        typeof (column as { get?: (i: number) => unknown }).get === 'function'
          ? (column as { get: (i: number) => unknown }).get(i)
          : (column as ArrayLike<unknown>)[i],
      )
    : Array.isArray(column)
      ? column
      : [];
  if (!values.length) return undefined;
  const seen = new Set<string>();
  for (const v of values) if (v !== null && v !== undefined) seen.add(String(v));
  return [...seen];
}

/**
 * A layer's params as the builder writes them, given the inline rows of the data it names (`rows`, a ranges node's):
 * a `ranges` layer with no `color` over rows with a `category` column is colored by it; a `field` form marked to list
 * gets its text column's values written (a categorical palette over the default `categories`).
 */
export function writeLayerTheme(
  params: Record<string, unknown>,
  rows: Readonly<Record<string, unknown>> | undefined,
): Record<string, unknown> {
  let color = params.color as Node | undefined;
  if (color === undefined && params.type === 'ranges' && rows && rows.category !== undefined)
    color = theme.field('category');
  // a variants layer: its consequence (else its category), as the layer read it before themes
  if (color === undefined && params.type === 'variants' && rows) {
    const name = rows.consequence !== undefined ? 'consequence' : rows.category !== undefined ? 'category' : undefined;
    if (name) color = theme.field(name);
  }
  const out = { ...params };
  for (const key of ['color', 'ink', 'size', 'opacity'] as const) {
    const v = key === 'color' ? color : (params[key] as Node | undefined);
    if (v === undefined) continue;
    out[key] = listed(v, rows);
  }
  return out;
}

function listed(v: Node, rows: Readonly<Record<string, unknown>> | undefined): Node {
  if (typeof v !== 'object' || v === null || v.kind !== 'field' || !v.params?.[ListMarker]) return unmark(v);
  const { [ListMarker]: _, ...params } = v.params as P;
  const values = rows ? distinct(rows[String(params.name)]) : undefined;
  const text = values && !values.every((x) => x !== '' && !Number.isNaN(Number(x)));
  if (!values || !text) return { ...v, params } as Node;
  return {
    ...v,
    params: { ...params, palette: theme.categorical({ colors: 'categories', values }) },
  } as Node;
}

/** A form with the listing marks of nested `field` forms taken out (only a slot's own form is listed). */
function unmark<T>(v: T): T {
  if (Array.isArray(v)) return v.map(unmark) as T;
  if (!v || typeof v !== 'object' || isColumn(v as Value)) return v;
  const out: Record<string, unknown> = {};
  for (const [k, x] of Object.entries(v)) if (k !== ListMarker) out[k] = unmark(x);
  return out as T;
}
