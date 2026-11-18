/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Validation against a language's description (_spec/README.md "Validation"): every issue, each with a stable code,
// its severity, a path, the address where the place has one, and an actionable hint. Messages name only the language's
// own vocabulary (its kinds, its params, their types); what the spec wrote is said by the place, never in the text.
// Ids are local: unique among siblings. A slot's node is validated against the language the slot names, looked up
// when it is met, so languages that embed each other validate each other with no import between them. Written from
// meso*'s schema core; columns, tables, extension kinds and the units hook are the delta.

import {
  type CategoryColumn,
  type Column,
  columnMisfit,
  columnType,
  columnValues,
  defineIssues,
  IdPattern,
  type Issue,
  isColumn,
  isNode,
  isRecord,
  kindNamespace,
  type Node,
  own,
  parseAddress,
  printAddress,
  type TreePath,
  type Value,
} from '../tree/index.js';
import { childIndex } from '../tree/address.js';
import type { Changes } from './changes.js';
import type {
  ColumnSpec,
  DataSpec,
  Fields,
  KindSpec,
  Language,
  ParamType,
  RefType,
  TableSpec,
  UnitsHook,
} from './language.js';

export const SchemaIssues = defineIssues('repr.schema', {
  'repr.schema.root.kind': { severity: 'error', description: "The root is not the language's root kind" },
  'repr.schema.kind.unknown': { severity: 'error', description: 'A kind the language does not have' },
  'repr.schema.kind.extension': {
    severity: 'warning',
    description: "A kind in a namespace the language does not have (an extension's): kept and skipped, not checked",
    hint: 'register the extension that defines it',
  },
  'repr.schema.parent.invalid': { severity: 'error', description: 'A node sits under a kind it may not sit under' },
  'repr.schema.id.missing': { severity: 'error', description: 'A node of a kind that needs an id has none' },
  'repr.schema.id.forbidden': { severity: 'error', description: 'A node of a kind that takes no id has one' },
  'repr.schema.id.duplicate': { severity: 'error', description: 'Two sibling nodes have the same id' },
  'repr.schema.param.unknown': { severity: 'error', description: 'A param (or record field) the kind does not have' },
  'repr.schema.param.missing': { severity: 'error', description: 'A required param (or record field) is absent' },
  'repr.schema.param.type': { severity: 'error', description: 'A value is not of its param type' },
  'repr.schema.param.enum': { severity: 'error', description: 'A string is not one of the values its param allows' },
  'repr.schema.param.pattern': { severity: 'error', description: 'A string does not match its param pattern' },
  'repr.schema.param.range': {
    severity: 'error',
    description: 'A number, or the length of a list or column, is outside its range',
  },
  'repr.schema.param.unit': {
    severity: 'error',
    description: "A quantity the units hook does not read in its slot's unit, or a quantity where a pure number is",
  },
  'repr.schema.column.type': {
    severity: 'error',
    description:
      'A column value is not of its declared element type: a typed array of another, or an item that does not fit',
  },
  'repr.schema.table.length': { severity: 'error', description: "A table's columns are not of equal length" },
  'repr.schema.table.column': {
    severity: 'error',
    description: "A table's extra column is neither a column nor a list",
  },
  'repr.schema.table.key': { severity: 'error', description: "A table's key column holds a value twice" },
  'repr.schema.slot.language': { severity: 'error', description: 'A slot names a language that is not registered' },
  'repr.schema.slot.kind': { severity: 'error', description: 'A slot holds a node of a kind it does not allow' },
  'repr.schema.ref.address': {
    severity: 'error',
    description: 'A reference is not an address (/a/b), or not a bare id where its slot resolves ids (:in)',
  },
  'repr.schema.args.count': {
    severity: 'error',
    description: "A data form holds fewer or more positional values than its kind's",
  },
  'repr.schema.args.shape': {
    severity: 'error',
    description: 'A data form has an id or children, or a node of another kind has args',
  },
  'repr.schema.limit_exceeded': {
    severity: 'warning',
    description: 'More issues than asked for: the rest are not reported',
  },
});
type Code = keyof typeof SchemaIssues.codes;

export interface ValidateOptions {
  /** The languages slots may name, resolved late; the validated language itself is always known. */
  readonly languages?: (id: string) => Language | undefined;
  /** Stop after this many issues (default 100). */
  readonly maxIssues?: number;
  /** The tie to a units table; without it a number slot takes no quantity string. */
  readonly units?: UnitsHook;
}

class Limit extends Error {}

/** Every issue of a spec against its language. */
export function validate(language: Language, root: Node, options: ValidateOptions = {}): Issue[] {
  const c = checker(language, options);
  c.run(() => {
    if (root.kind !== language.root)
      c.report('repr.schema.root.kind', `The root must be a ${language.root}`, undefined, 'kind');
    c.node(language, root, undefined, []);
  });
  return c.issues;
}

/** Every issue of a value against a param type (a default, a value from elsewhere); slots resolve through `options`. */
export function validateValue(type: ParamType, value: Value, options: ValidateOptions = {}): Issue[] {
  const c = checker(undefined, options);
  c.run(() => c.value(type, value));
  return c.issues;
}

/** A node found by address: its step path from the root, and its parent (none for the root). */
interface Located {
  readonly node: Node;
  readonly parent: Node | undefined;
  readonly path: TreePath;
}

function locate(root: Node, address: readonly string[]): Located | undefined {
  let node = root;
  let parent: Node | undefined;
  const path: (string | number)[] = [];
  for (const id of address) {
    const i = childIndex(node.children, id);
    if (i < 0) return undefined;
    parent = node;
    node = node.children![i];
    path.push('children', i);
  }
  return { node, parent, path };
}

const addressKey = (a: readonly string[]) => a.join('/');
const isUnder = (a: readonly string[], prefix: readonly string[]) =>
  a.length >= prefix.length && prefix.every((s, i) => a[i] === s);

/**
 * The issues a change introduced, and only those (_spec/README.md "Patches", Validation). `root` is the new root,
 * `changes` the net difference that made it. A made node is checked whole (its subtree, its parent rule, a root's
 * kind); a node whose params changed, for the changed keys and that no unset key is required; a moved node, for its
 * parent rule; a parent whose child list changed, for the uniqueness of the ids the change brought into it. Nothing
 * else is read, so an error present before the change is not reported again.
 */
export function validateChanges(
  language: Language,
  root: Node,
  changes: Changes,
  options: ValidateOptions = {},
): Issue[] {
  const c = checker(language, options);
  // under a kind the language lacks there is no parent rule, as `validate` reads it
  const kindOf = (n: Node | undefined) => (n && own(language.kinds, n.kind) ? n.kind : undefined);
  const made = changes.made;
  const inMade = (a: readonly string[]) => made.some((m) => isUnder(a, m));
  c.run(() => {
    for (const a of made) {
      const at = locate(root, a);
      if (!at) continue;
      c.at(at.path, a, () => {
        if (!at.parent && at.node.kind !== language.root)
          c.report('repr.schema.root.kind', `The root must be a ${language.root}`, undefined, 'kind');
        c.node(language, at.node, kindOf(at.parent), [...a]);
      });
    }
    for (const change of changes.params) {
      if (change.where !== 'params' || inMade(change.at)) continue;
      const at = locate(root, change.at);
      const spec = at && own(language.kinds, at.node.kind);
      if (!at || !spec) continue; // an extension's kind: its params are not checked
      c.at([...at.path, 'params'], change.at, () => {
        for (const k of change.keys) {
          const v = at.node.params ? own(at.node.params, k) : undefined;
          const p = spec.params ? own(spec.params, k) : undefined;
          if (v === undefined) {
            if (p?.required) c.report('repr.schema.param.missing', `The param "${k}" is required`, p.doc, k);
          } else if (!p) {
            const names = Object.keys(spec.params ?? {});
            c.report('repr.schema.param.unknown', 'Not a param it has', names.slice(0, 20).join(', ') || 'none', k);
          } else {
            c.path.push(k);
            c.value(p, v);
            c.path.pop();
          }
        }
      });
    }
    for (const { to } of changes.moved) {
      if (inMade(to)) continue;
      const at = locate(root, to);
      const spec = at && own(language.kinds, at.node.kind);
      if (!at || !spec) continue;
      c.at(at.path, to, () => c.parentRule(at.node.kind, spec, kindOf(at.parent)));
    }
    // the ids the change brought into a child list are unique among their siblings
    const fresh = new Set<string>([...made.map(addressKey), ...changes.moved.map((m) => addressKey(m.to))]);
    for (const p of changes.children) {
      const at = locate(root, p);
      if (!at?.node.children) continue;
      const seen = new Set<string>();
      at.node.children.forEach((child, i) => {
        if (child.id === undefined) return;
        if (!seen.has(child.id)) seen.add(child.id);
        else if (fresh.has(addressKey([...p, child.id])))
          c.at([...at.path, 'children', i], undefined, () =>
            c.report('repr.schema.id.duplicate', 'An earlier sibling has the same id', undefined, 'id'),
          );
      });
    }
  });
  return c.issues;
}

const infinite = (x: number) => x === Number.POSITIVE_INFINITY || x === Number.NEGATIVE_INFINITY;

/** The length of a table's column: a column's or a list's; undefined for anything else. */
function lengthOf(v: Value): number | undefined {
  if (isColumn(v)) return v.length;
  return Array.isArray(v) ? v.length : undefined;
}

function checker(language: Language | undefined, options: ValidateOptions) {
  const issues: Issue[] = [];
  const max = options.maxIssues ?? 100;
  const path: (string | number)[] = [];
  /** The address of the node being checked, or undefined where it has none (a node without an id, a slot). */
  let address: string[] | undefined;

  const report = (code: Code, message: string, hint?: string, ...at: (string | number)[]) => {
    if (issues.length === max) {
      issues.push({ code: 'repr.schema.limit_exceeded', severity: 'warning', message: `Stopped after ${max} issues` });
      throw new Limit();
    }
    const issue: { -readonly [K in keyof Issue]: Issue[K] } = {
      code,
      severity: SchemaIssues.codes[code].severity,
      message,
      path: [...path, ...at],
    };
    if (hint) issue.hint = hint;
    if (address) issue.address = printAddress(address);
    issues.push(issue);
  };
  const resolve = (id: string) => (id === language?.id ? language : options.languages?.(id));
  const list = (names: readonly string[]) => names.slice(0, 20).join(', ') + (names.length > 20 ? ', …' : '');
  /**
   * The value being read is a list's item (a union's members read it as one): its numbers are finite unless their type
   * says `infinite`. Set by a list per item, taken and cleared by `value` on entry.
   */
  let listItem = false;

  function node(lang: Language, n: Node, parent: string | undefined, at: string[] | undefined) {
    const outer = address;
    address = at;
    try {
      nodeIn(lang, n, parent);
    } finally {
      address = outer;
    }
  }

  function nodeIn(lang: Language, n: Node, parent: string | undefined) {
    const spec = own(lang.kinds, n.kind);
    if (!spec) {
      // an extension's kind (format.md "Unknown kinds"): kept, reported, its params skipped; its children are the
      // language's, checked with no parent rule (what may sit under it is the extension's)
      if (kindNamespace(n.kind) !== undefined && kindNamespace(n.kind) !== lang.id)
        report('repr.schema.kind.extension', 'A kind of an extension this reader lacks', undefined, 'kind');
      else report('repr.schema.kind.unknown', `Not a kind of ${lang.id}`, list(Object.keys(lang.kinds)), 'kind');
      children(lang, n, undefined);
      return;
    }
    parentRule(n.kind, spec, parent);
    if (spec.data) {
      if (n.id !== undefined || n.children)
        report('repr.schema.args.shape', `${n.kind} is a data form: values, no id or children`, undefined);
    } else if (n.id !== undefined) {
      if (spec.ids === 'forbidden') report('repr.schema.id.forbidden', `${n.kind} takes no id`, undefined, 'id');
    } else if (spec.ids === 'required') report('repr.schema.id.missing', `${n.kind} needs an id`);
    if (n.args && !spec.data)
      report('repr.schema.args.shape', `${n.kind} takes no positional values`, undefined, 'args');
    path.push('params');
    fields(spec.params ?? {}, n.params ?? {}, 'param');
    path.pop();
    if (spec.data) args(spec.data, n.args ?? []);
    children(lang, n, n.kind);
  }

  function parentRule(kind: string, spec: KindSpec, parent: string | undefined) {
    if (parent !== undefined && !spec.parents.includes(parent)) {
      const where = spec.parents.length ? `under ${list(spec.parents)}` : 'only at the root';
      report('repr.schema.parent.invalid', `${kind} may not sit under ${parent}`, `${kind} sits ${where}`, 'kind');
    }
  }

  /** Runs `f` at a place in the tree, its step path and address (where `validateChanges` starts each check). */
  function at(prefix: readonly (string | number)[], where: readonly string[] | undefined, f: () => void) {
    const outer = address;
    const saved = path.length;
    path.push(...prefix);
    address = where ? [...where] : undefined;
    try {
      f();
    } finally {
      address = outer;
      path.length = saved;
    }
  }

  function children(lang: Language, n: Node, parent: string | undefined) {
    if (!n.children) return;
    path.push('children');
    // ids are local: unique among siblings, so a node's address (its ids from the root) is unique
    const seen = new Set<string>();
    for (let i = 0; i < n.children.length; i++) {
      path.push(i);
      const child = n.children[i];
      const id = child.id;
      if (id !== undefined) {
        if (seen.has(id)) report('repr.schema.id.duplicate', 'An earlier sibling has the same id', undefined, 'id');
        seen.add(id);
      }
      const childAddress = address && id !== undefined && IdPattern.test(id) ? [...address, id] : undefined;
      node(lang, child, parent, childAddress);
      path.pop();
    }
    path.pop();
  }

  function args(d: DataSpec, values: readonly Value[]) {
    const min = d.min ?? d.args.filter((a) => a.required).length;
    const max = d.max ?? d.args.length;
    const names = d.args.map((a) => a.name).join(' ') + (d.rest ? ' …' : '');
    if (values.length < min || values.length > max)
      report('repr.schema.args.count', `Expected ${min} to ${max} positional values`, names, 'args');
    path.push('args');
    for (let i = 0; i < Math.min(values.length, max); i++) {
      if (i >= d.args.length && !d.rest) break;
      const spec = d.args[Math.min(i, d.args.length - 1)];
      path.push(i);
      value(spec, values[i]);
      path.pop();
    }
    path.pop();
  }

  function fields(specs: Fields, values: Readonly<Record<string, Value>>, what: 'param' | 'field' | 'column') {
    for (const key of Object.keys(values)) {
      if (values[key] === undefined) continue;
      if (!own(specs, key)) {
        // a table's extra columns are its custom fields, checked by `table`
        if (what !== 'column')
          report('repr.schema.param.unknown', `Not a ${what} it has`, list(Object.keys(specs)) || 'none', key);
      }
    }
    for (const [key, spec] of Object.entries(specs)) {
      const v = own(values, key);
      if (v === undefined) {
        if (spec.required) report('repr.schema.param.missing', `The ${what} "${key}" is required`, spec.doc, key);
        continue;
      }
      path.push(key);
      value(spec, v);
      path.pop();
    }
  }

  /** Reports and returns whether `v` is of type `t`. */
  function value(t: ParamType, v: Value): boolean {
    const item = listItem;
    listItem = false;
    const wrong = () => (report('repr.schema.param.type', `Expected ${typeName(t)}`, undefined), false);
    switch (t.type) {
      case 'any':
        return true;
      case 'string':
        if (typeof v !== 'string') return wrong();
        if (t.enum && !t.enum.includes(v))
          return report('repr.schema.param.enum', 'Not an allowed value', list(t.enum)), false;
        if (t.pattern !== undefined && !pattern(t.pattern).test(v))
          return report('repr.schema.param.pattern', `Does not match ${t.pattern}`), false;
        return true;
      case 'number':
      case 'integer': {
        const finite = item && !t.infinite;
        if (typeof v === 'string') return quantity(t, v, finite);
        if (typeof v !== 'number' || (t.type === 'integer' && !Number.isInteger(v))) return wrong();
        // a list's item is finite unless its type takes the infinities (NaN is the range's)
        if (finite && infinite(v)) return notFinite();
        return inRange(v, t.min, t.max, 'value');
      }
      case 'boolean':
        return typeof v === 'boolean' || wrong();
      case 'list': {
        if (!Array.isArray(v)) return wrong();
        let ok = inRange(v.length, t['min-length'], t['max-length'], 'length');
        for (let i = 0; i < v.length; i++) {
          path.push(i);
          listItem = true;
          ok = value(t.of, v[i]) && ok;
          path.pop();
        }
        return ok;
      }
      case 'record':
        if (!isRecord(v)) return wrong();
        fields(t.fields, v, 'field');
        return true;
      case 'map': {
        if (!isRecord(v)) return wrong();
        let keys = 0;
        if (t['max-keys'] !== undefined) for (const key in v) if (Object.hasOwn(v, key)) keys++;
        let ok = inRange(keys, 0, t['max-keys'], 'key count');
        for (const key in v) {
          if (!Object.hasOwn(v, key)) continue;
          path.push(key);
          ok = value(t.of, v[key]) && ok;
          path.pop();
        }
        return ok;
      }
      case 'union': {
        // its members read the value as it came: a list's item stays one
        const member = (m: ParamType) => ((listItem = item), value(m, v));
        // a node is read by the member that takes its kind, its issues its own
        if (isNode(v) && t.of.some((m) => m.type === 'node')) {
          const m = unionNodeMember(t.of, v.kind, resolve);
          if (m) return member(m);
          const kinds = t.of.flatMap((x) => (x.type === 'node' ? (x.kinds ?? []) : []));
          if (kinds.length && !t.of.some((x) => x.type === 'node' && !x.kinds))
            return report('repr.schema.slot.kind', 'Not a kind this slot takes', list(kinds)), false;
        }
        // the first member it fits, checked quietly; else one issue for the union
        if (t.of.some((m) => ((listItem = item), fits(m, v)))) return true;
        // a string where a reference is taken: what a reference's spelling says helps more than "fits no member"
        const ref = typeof v === 'string' ? t.of.find((m) => m.type === 'ref') : undefined;
        if (ref) return member(ref);
        return wrong();
      }
      case 'ref':
        if (typeof v !== 'string') return wrong();
        // the spelling only: what it names is a resolve pass's
        if (t.in !== undefined ? IdPattern.test(v) : v.startsWith('/') && parseAddress(v).pattern !== undefined)
          return true;
        return (
          report('repr.schema.ref.address', t.in !== undefined ? 'Expected an id' : 'Expected an address', refHint(t)),
          false
        );
      case 'node': {
        if (!isNode(v)) return wrong();
        const lang = resolve(t.language);
        if (!lang) return report('repr.schema.slot.language', `No language "${t.language}" is registered`), false;
        if (t.kinds && !t.kinds.includes(v.kind))
          return report('repr.schema.slot.kind', 'Not a kind this slot takes', list(t.kinds)), false;
        node(lang, v, undefined, undefined);
        return true;
      }
      case 'column':
        return column(t, v);
      case 'table':
        return table(t, v);
    }
  }

  /** A typed array of the declared element type, or a plain array whose every item fits it. */
  function column(t: ColumnSpec, v: Value): boolean {
    const have = columnType(v);
    if (have === undefined && !Array.isArray(v))
      return report('repr.schema.param.type', `Expected ${typeName(t)}`), false;
    if (have !== undefined && have !== t.of) {
      report(
        'repr.schema.column.type',
        `Expected a ${t.of} column, not a ${have} one`,
        'a plain array, or its typed array',
      );
      return false;
    }
    const length = (v as Column | readonly Value[]).length;
    let ok = inRange(length, t['min-length'], t['max-length'], 'length');
    if (have === 'category' && !t.missing && (v as CategoryColumn).codes.includes(-1)) {
      report('repr.schema.column.type', 'A missing value where none is allowed', 'every value a string');
      ok = false;
    }
    if (have === undefined) {
      const i = columnMisfit(t.of, v as readonly unknown[], t.missing);
      if (i >= 0) {
        report(
          'repr.schema.column.type',
          `Not a ${t.of} value${t.of === 'category' && !t.missing ? ' (no missing)' : ''}`,
          columnHint(t),
          i,
        );
        ok = false;
      }
    }
    return ok;
  }

  function table(t: TableSpec, v: Value): boolean {
    if (!isRecord(v)) return report('repr.schema.param.type', `Expected ${typeName(t)}`), false;
    const before = issues.length;
    fields(t.columns, v, 'column');
    let rows: number | undefined;
    for (const key of Object.keys(v)) {
      const x = v[key];
      if (x === undefined) continue;
      const n = lengthOf(x);
      if (n === undefined) {
        // a declared column's type issue is said by `fields`; an extra one is a column or a list
        if (!own(t.columns, key))
          report('repr.schema.table.column', 'An extra column is a column or a list', undefined, key);
        continue;
      }
      if (rows === undefined) rows = n;
      else if (n !== rows)
        report('repr.schema.table.length', `Expected ${rows} rows, as the first column`, undefined, key);
    }
    if (t.key !== undefined && own(v, t.key) !== undefined) {
      const k = own(v, t.key)!;
      const values = isColumn(k) ? columnValues(k) : Array.isArray(k) ? k : [];
      const seen = new Set<unknown>();
      for (let i = 0; i < values.length; i++) {
        if (seen.has(values[i])) {
          report('repr.schema.table.key', 'A key value given twice', undefined, t.key, i);
          break;
        }
        seen.add(values[i]);
      }
    }
    return issues.length === before;
  }

  /** A string in a number slot: a quantity the units hook reads, in the slot's range once in its unit. */
  function quantity(t: Extract<ParamType, { type: 'number' | 'integer' }>, v: string, finite: boolean): boolean {
    const q = options.units?.quantity(v, t.unit);
    if (q === undefined) return report('repr.schema.param.type', `Expected ${typeName(t)}`), false;
    if ('message' in q) return report('repr.schema.param.unit', q.message, q.hint), false;
    if (finite && infinite(q.value)) return notFinite();
    return Number.isNaN(q.value) || inRange(q.value, t.min, t.max, 'value');
  }

  function notFinite(): false {
    report('repr.schema.param.range', 'The value is not finite', 'a finite number');
    return false;
  }

  function inRange(v: number, min: number | undefined, max: number | undefined, what: string) {
    if ((min === undefined || v >= min) && (max === undefined || v <= max) && !Number.isNaN(v)) return true;
    report('repr.schema.param.range', `The ${what} is outside [${min ?? '-∞'}, ${max ?? '∞'}]`);
    return false;
  }

  /** Whether `v` fits `t`, reporting nothing. */
  function fits(t: ParamType, v: Value): boolean {
    const before = issues.length;
    const saved = path.length;
    try {
      value(t, v);
    } catch (e) {
      if (!(e instanceof Limit)) throw e;
    }
    const ok = issues.length === before;
    issues.length = before;
    path.length = saved;
    return ok;
  }

  const run = (f: () => void) => {
    try {
      f();
    } catch (e) {
      if (!(e instanceof Limit)) throw e;
    }
  };
  return { issues, report, node, value, run, at, parentRule, path };
}

/**
 * The member of a union validation reads a value by: a node by the member that takes its kind; else the first member
 * it fits quietly; else, for a string, the union's ref (its spelling is then what is wrong). Normalization and the
 * passes over refs read a union the same way, so they never disagree with validation on which member a value is.
 */
export function unionMember(of: readonly ParamType[], v: Value, options: ValidateOptions = {}): ParamType | undefined {
  if (isNode(v) && of.some((m) => m.type === 'node')) {
    const m = unionNodeMember(of, v.kind, (id) => options.languages?.(id));
    if (m) return m;
  }
  const quiet = { ...options, maxIssues: 1 };
  const fit = of.find((m) => validateValue(m, v, quiet).length === 0);
  if (fit) return fit;
  return typeof v === 'string' ? of.find((m) => m.type === 'ref') : undefined;
}

/** The member of a union that holds a node of `kind`: the node member whose `kinds` name it, else whose language has it. */
export function unionNodeMember(
  of: readonly ParamType[],
  kind: string,
  resolve: (id: string) => Language | undefined,
): Extract<ParamType, { type: 'node' }> | undefined {
  const slots = of.filter((m): m is Extract<ParamType, { type: 'node' }> => m.type === 'node');
  return (
    slots.find((m) => m.kinds?.includes(kind)) ??
    slots.find((m) => !m.kinds && Object.hasOwn(resolve(m.language)?.kinds ?? {}, kind))
  );
}

function columnHint(t: ColumnSpec): string {
  if (t.of === 'category') return t.missing ? 'a string or nil' : 'a string';
  if (t.of === 'float32' || t.of === 'float64') return 'a number, or nil (NaN) for a missing value';
  return `an integer of ${t.of}`;
}

/** What a reference's slot takes, for a hint: `an address of a sequence`, `an id of a state in states`. */
export function refHint(t: RefType): string {
  const of = t.kinds?.length ? ` of a ${t.kinds.slice(0, 8).join(' or ')}` : '';
  return t.in !== undefined ? `an id${of} in ${t.in}` : `an address${of}: /a/b`;
}

const Patterns = new Map<string, RegExp>();
function pattern(source: string) {
  let re = Patterns.get(source);
  if (!re) Patterns.set(source, (re = new RegExp(source)));
  return re;
}

/** A param type as people read it (`integer ≥ 1`, `"a" | "b"`, `list of string`, `int32 column`). */
export function typeName(t: ParamType): string {
  switch (t.type) {
    case 'string':
      if (t.enum) return t.enum.map((e) => JSON.stringify(e)).join(' | ');
      return t.pattern !== undefined ? `string matching ${t.pattern}` : 'string';
    case 'number':
    case 'integer': {
      const r = [t.min !== undefined ? `≥ ${t.min}` : '', t.max !== undefined ? `≤ ${t.max}` : ''].filter(Boolean);
      return [t.type, ...r].join(' ') + (t.unit !== undefined ? ` in ${t.unit}` : '');
    }
    case 'list':
      return `list of ${typeName(t.of)}`;
    case 'record':
      return `record (${Object.keys(t.fields).join(', ')})`;
    case 'map':
      return `map of ${typeName(t.of)}`;
    case 'union':
      return t.of.map(typeName).join(' or ');
    case 'node':
      return `${t.language} node${t.kinds ? ` (${t.kinds.join(', ')})` : ''}`;
    case 'ref': {
      const kinds = t.kinds ? ` (${t.kinds.join(', ')})` : '';
      return t.in !== undefined ? `id in ${t.in}${kinds}` : `address${kinds}`;
    }
    case 'column':
      return `${t.of} column`;
    case 'table':
      return `table (${Object.keys(t.columns).join(', ')})`;
    default:
      return t.type;
  }
}
