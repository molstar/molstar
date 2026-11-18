/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Normalization on read (_spec/README.md "Normalization"; language-description.md "Codecs stay schema-free"): a codec
// decodes a column as a plain array, and normalizing turns every value of a `column` param (a table's declared columns
// among them) into its declared typed array. Equality, hashing and printing are defined on normalized trees. The read
// pipeline is decode, migrate, then validate and normalize (`validateAndNormalize`). Unchanged subtrees are shared:
// a tree that is already normalized comes back as the same object.

import {
  type Column,
  columnType,
  hasErrors,
  type Issue,
  isColumn,
  isNode,
  isRecord,
  type Node,
  own,
  toColumn,
  type Value,
} from '../tree/index.js';
import type { Fields, KindSpec, Language, ParamType } from './language.js';
import { type ValidateOptions, validate, validateValue } from './validate.js';

export interface NormalizeOptions {
  /** The languages slots may name, resolved late, as validation's. */
  readonly languages?: (id: string) => Language | undefined;
}

/** Whether values of a type may hold something to normalize: a column, a table, a slot's node, a union's member. */
const Needs = new WeakMap<object, boolean>();
function needs(t: ParamType): boolean {
  let n = Needs.get(t);
  if (n !== undefined) return n;
  switch (t.type) {
    case 'column':
    case 'table':
    case 'node':
      n = true;
      break;
    case 'list':
    case 'map':
      n = needs(t.of);
      break;
    case 'record':
      n = Object.values(t.fields).some(needs);
      break;
    case 'union':
      n = t.of.some(needs);
      break;
    default:
      n = false;
  }
  Needs.set(t, n);
  return n;
}

/**
 * The tree with every column param as its declared typed array; best effort: a value that does not fit its type is
 * left as it is (validation says why). Kinds the language lacks are left as they are.
 */
export function normalize(language: Language, root: Node, options: NormalizeOptions = {}): Node {
  const resolve = (id: string) => (id === language.id ? language : options.languages?.(id));
  return normalizer(resolve).node(language, root);
}

/**
 * A value normalized against its param type, as `normalize` reads a param of that type (what a patch's `set` writes);
 * slots resolve through `options.languages`. Best effort: a value that does not fit is left as it is.
 */
export function normalizeValue(type: ParamType, value: Value, options: NormalizeOptions = {}): Value {
  return normalizer((id) => options.languages?.(id)).value(type, value);
}

function normalizer(resolve: (id: string) => Language | undefined) {
  function node(lang: Language, n: Node): Node {
    // an unknown (extension) kind is kept as it is; its children are still the language's
    const spec: KindSpec | undefined = own(lang.kinds, n.kind);
    const params = spec?.params && n.params ? record(spec.params, n.params) : n.params;
    let args = n.args;
    if (spec?.data && n.args) {
      const specs = spec.data.args;
      const out = n.args.map((v, i) => (specs.length ? value(specs[Math.min(i, specs.length - 1)], v) : v));
      if (out.some((v, i) => v !== n.args![i])) args = out;
    }
    let children = n.children;
    if (n.children) {
      const out = n.children.map((c) => node(lang, c));
      if (out.some((c, i) => c !== n.children![i])) children = out;
    }
    if (params === n.params && args === n.args && children === n.children) return n;
    const next: { -readonly [K in keyof Node]: Node[K] } = { ...n };
    if (params) next.params = params;
    if (args) next.args = args;
    if (children) next.children = children;
    return next;
  }

  function record(specs: Fields, values: Readonly<Record<string, Value>>): Readonly<Record<string, Value>> {
    let out: Record<string, Value> | undefined;
    for (const key of Object.keys(values)) {
      const spec = own(specs, key);
      const v = values[key];
      if (!spec || v === undefined || !needs(spec)) continue;
      const w = value(spec, v);
      if (w !== v) (out ??= { ...values })[key] = w;
    }
    return out ?? values;
  }

  function value(t: ParamType, v: Value): Value {
    if (!needs(t)) return v;
    switch (t.type) {
      case 'column':
        return toColumn(t.of, v, t.missing) ?? v;
      case 'table':
        return isRecord(v) ? record(t.columns, v) : v;
      case 'record':
        return isRecord(v) ? record(t.fields, v) : v;
      case 'list': {
        if (!Array.isArray(v)) return v;
        const out = v.map((x) => value(t.of, x));
        return out.some((x, i) => x !== v[i]) ? out : v;
      }
      case 'map': {
        if (!isRecord(v)) return v;
        let out: Record<string, Value> | undefined;
        for (const key of Object.keys(v)) {
          const w = value(t.of, v[key]);
          if (w !== v[key]) (out ??= { ...v })[key] = w;
        }
        return out ?? v;
      }
      case 'union': {
        // the first member the value fits, as validation reads it
        const m = t.of.find((x) => validateValue(x, v, { languages: resolve, maxIssues: 1 }).length === 0);
        return m ? value(m, v) : v;
      }
      case 'node': {
        const lang = isNode(v) ? resolve(t.language) : undefined;
        return lang && isNode(v) ? node(lang, v) : v;
      }
      default:
        return v;
    }
  }

  return { node, value };
}

/**
 * Validation, then normalization when there is no error (warnings and infos kept): the read pipeline's last step.
 * `node` is absent when the tree is invalid.
 */
export function validateAndNormalize(
  language: Language,
  root: Node,
  options: ValidateOptions = {},
): { readonly node?: Node; readonly issues: readonly Issue[] } {
  const issues = validate(language, root, options);
  if (hasErrors(issues)) return { issues };
  return { node: normalize(language, root, options), issues };
}

/**
 * Where a tree's columns sit in declared column params (a table's declared columns included): each container (a
 * params record, a table, a list, a data form's args) with the keys or indices that hold one. A codec writes a column
 * there as a plain array (normalization types it on read) and tags every other column (`$column`), so a column in
 * `custom`, an `any` param or an undeclared table column comes back typed (_spec/README.md "Codecs").
 */
export type ColumnPlaces = WeakMap<object, Set<string | number>>;

export function declaredColumnPlaces(language: Language, root: Node, options: NormalizeOptions = {}): ColumnPlaces {
  const places: ColumnPlaces = new WeakMap();
  const resolve = (id: string) => (id === language.id ? language : options.languages?.(id));
  const mark = (container: object, key: string | number) => {
    let keys = places.get(container);
    if (!keys) places.set(container, (keys = new Set()));
    keys.add(key);
  };

  function node(lang: Language, n: Node) {
    const spec: KindSpec | undefined = own(lang.kinds, n.kind);
    if (spec?.params && n.params) record(spec.params, n.params);
    if (spec?.data && n.args) {
      const specs = spec.data.args;
      if (specs.length) n.args.forEach((v, i) => value(specs[Math.min(i, specs.length - 1)], v, n.args!, i));
    }
    for (const c of n.children ?? []) node(lang, c);
  }

  function record(specs: Fields, values: Readonly<Record<string, Value>>) {
    for (const key of Object.keys(values)) {
      const spec = own(specs, key);
      if (spec && values[key] !== undefined && needs(spec)) value(spec, values[key], values, key);
    }
  }

  function value(t: ParamType, v: Value, container: object, key: string | number) {
    switch (t.type) {
      case 'column':
        if (columnType(v) === t.of) mark(container, key);
        return;
      case 'table':
      case 'record':
        if (isRecord(v)) record(t.type === 'table' ? t.columns : t.fields, v);
        return;
      case 'list':
        if (Array.isArray(v) && needs(t.of)) v.forEach((x, i) => value(t.of, x, v, i));
        return;
      case 'map':
        if (isRecord(v) && needs(t.of)) for (const k of Object.keys(v)) value(t.of, v[k], v, k);
        return;
      case 'union': {
        const m = isColumn(v)
          ? t.of.find((x) => x.type === 'column' && x.of === columnType(v as Column))
          : t.of.find((x) => needs(x) && validateValue(x, v, { languages: resolve, maxIssues: 1 }).length === 0);
        if (m) value(m, v, container, key);
        return;
      }
      case 'node': {
        const lang = isNode(v) ? resolve(t.language) : undefined;
        if (lang) node(lang, v as Node);
        return;
      }
    }
  }

  node(language, root);
  return places;
}
