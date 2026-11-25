/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Predicates (_spec/README.md "Queries", Language): the expression forms a
// test is written in, with positional operands (the author, 2026-11-20, "Accept, but allow positional operands"; held
// in the node's `args`): the operators `=`, `!=`, `<`, `<=`, `>`, `>=`, `in-range`, `in`, `and`, `or`, `not`, and
// `prop` and `value` reading properties (_spec/README.md "Queries", Properties). Compiled once a form, typed when read:
// an operator over a word and a number, a test that is not boolean, is `seq.query.form.type`; a property unknown is a
// param error, a namespaced one nobody registered `seq.query.extension.missing`. An element lacking a property gives no
// value, and a test without a value is false: the element yields nothing, counted for one warning a property.

import { type Atom, isNode, type Node, type Value } from '@molstar/core/repr/tree';
import { fail } from './form-kit.js';
import { findProperty, levelRank, type PropertyDef, type PropertyTable, type PropertyType } from './properties.js';
import type { PropValue, Reader } from './structure.js';

/** A value's type: a property's, a list of atoms (`in`'s second operand), or `nil`. */
export type ExprType = PropertyType | 'list' | 'nil';

/** A value an expression gives an element: what a property gives, a list, or undefined (lacking). */
export type ExprValue = PropValue | readonly Atom[];

/** An extension's operator (a `query-form` entry of role `operator`): over values, its result a value of `type`. */
export interface QueryOperator {
  readonly name: string;
  readonly type: PropertyType;
  evaluate(args: readonly ExprValue[]): Atom | undefined;
}

/** The readers of one frame, by property; undefined where the frame lacks it. */
export type Readers = (d: PropertyDef) => Reader | undefined;

/** What a predicate counts as it runs: elements lacking each property, by its name. */
export type Absent = Map<string, number>;

interface Compiled {
  readonly type: ExprType;
  /** The finest level it reads (`levelRank`), -1 when none. */
  readonly level: number;
  /** A literal's value (coerced to the other operand's type by comparisons). */
  readonly literal?: Atom | readonly Atom[];
  ev(rs: Readers, r: number, a: number): ExprValue;
}

/** A compiled test: the level it reads, the properties, and its value per element. */
export interface Predicate {
  readonly level: number;
  readonly uses: readonly PropertyDef[];
  test(rs: Readers, r: number, a: number): boolean;
}

export interface CompileContext {
  readonly properties: PropertyTable;
  readonly operators?: ReadonlyMap<string, QueryOperator>;
}

/** The expression forms' kinds, built in. */
export const ExpressionKinds = ['=', '!=', '<', '<=', '>', '>=', 'in-range', 'in', 'and', 'or', 'not', 'prop', 'value'];

const isNumeric = (t: ExprType) => t === 'integer' || t === 'number';

function literalType(v: Value): ExprType {
  if (v === null) return 'nil';
  if (typeof v === 'boolean') return 'boolean';
  if (typeof v === 'number') return Number.isInteger(v) ? 'integer' : 'number';
  if (typeof v === 'string') return 'word';
  if (Array.isArray(v)) return 'list';
  return fail('seq.query.form.type', 'An operand is not an atom, a list or an expression');
}

/** A literal read as a property's type, or undefined when its spelling is not one. */
export function coerce(v: Atom, type: PropertyType): Atom | undefined {
  if (v === null) return null;
  if (type === 'word' || type === 'words')
    return typeof v === 'string' ? v : typeof v === 'number' ? `${v}` : undefined;
  if (type === 'boolean') return typeof v === 'boolean' ? v : undefined;
  const n = typeof v === 'number' ? v : typeof v === 'string' && v.trim() !== '' ? Number(v) : Number.NaN;
  if (Number.isNaN(n) || (type === 'integer' && !Number.isInteger(n))) return undefined;
  return n;
}

const args = (form: Node) => form.args ?? [];

function arity(form: Node, min: number, max = min) {
  const n = args(form).length;
  if (n < min || n > max)
    fail(
      'seq.query.form.param',
      `A form of kind ${form.kind} takes ${min === max ? min : `${min} to ${max}`} operands`,
    );
}

/** Compiles an expression; the properties it reads are added to `uses`, what lacks counted in `absent`. */
function compile(
  cx: CompileContext,
  v: Value,
  uses: Set<PropertyDef>,
  absent: Absent,
  values: Map<string, PropertyDef>,
): Compiled {
  if (!isNode(v)) {
    const type = literalType(v);
    const lit = v as Atom | readonly Atom[];
    return { type, level: -1, literal: lit, ev: () => lit };
  }
  const k = v.kind;
  const sub = (x: Value) => compile(cx, x, uses, absent, values);
  if (k === 'prop' || k === 'value') {
    const d = k === 'prop' ? propOf(cx, v) : valueOf(v, values);
    uses.add(d);
    const name = d.name;
    return {
      type: d.type,
      level: levelRank(d.level),
      ev: (rs, r, a) => {
        const rd = rs(d);
        const x = rd?.(r, a);
        if (x === undefined) absent.set(name, (absent.get(name) ?? 0) + 1);
        return x;
      },
    };
  }
  const ops = args(v).map(sub);
  const level = Math.max(-1, ...ops.map((o) => o.level));
  switch (k) {
    case 'and':
    case 'or':
    case 'not': {
      if (k === 'not') arity(v, 1);
      else arity(v, 1, Number.MAX_SAFE_INTEGER);
      for (const o of ops) if (o.type !== 'boolean') typeError(k, o.type);
      const and = k === 'and';
      if (k === 'not') return { type: 'boolean', level, ev: (rs, r, a) => bool(ops[0].ev(rs, r, a), (x) => !x) };
      return {
        type: 'boolean',
        level,
        ev: (rs, r, a) => {
          for (const o of ops) {
            const x = o.ev(rs, r, a);
            if (x === undefined) return undefined;
            if (x === !and) return !and;
          }
          return and;
        },
      };
    }
    case '=':
    case '!=':
      arity(v, 2);
      return equality(k, ops[0], ops[1], level);
    case '<':
    case '<=':
    case '>':
    case '>=': {
      arity(v, 2);
      const [x, y] = numeric(k, ops);
      const cmp = Comparisons[k];
      return { type: 'boolean', level, ev: (rs, r, a) => both(x.ev(rs, r, a), y.ev(rs, r, a), cmp) };
    }
    case 'in-range': {
      arity(v, 3);
      const [x, lo, hi] = numeric(k, ops);
      return {
        type: 'boolean',
        level,
        ev: (rs, r, a) => {
          const p = x.ev(rs, r, a);
          const l = lo.ev(rs, r, a);
          const h = hi.ev(rs, r, a);
          return typeof p === 'number' && typeof l === 'number' && typeof h === 'number' ? p >= l && p <= h : undefined;
        },
      };
    }
    case 'in': {
      arity(v, 2);
      const [x, list] = ops;
      if (list.type !== 'list' || !Array.isArray(list.literal)) return typeError(k, list.type);
      const t = x.type;
      if (t === 'list' || t === 'nil') return typeError(k, t);
      const set = new Set<Atom>();
      for (const y of list.literal as readonly Atom[]) {
        const c = isProp(t) ? coerce(y, t) : y;
        if (c === undefined) typeError(k, literalType(y));
        set.add(c as Atom);
      }
      return { type: 'boolean', level, ev: (rs, r, a) => member(x.ev(rs, r, a), set) };
    }
    default: {
      const op = cx.operators?.get(k);
      if (op)
        return { type: op.type, level, ev: (rs, r, a) => op.evaluate(ops.map((o) => o.ev(rs, r, a))) ?? undefined };
      if (k.includes(':'))
        return fail('seq.query.extension.missing', `No extension provides the operator ${k}`, {
          hint: `namespace: ${k.slice(0, k.indexOf(':'))}`,
        });
      return fail('seq.query.form.param', `A predicate takes expression forms; ${k} is none`, {
        hint: `forms: ${ExpressionKinds.join(' ')}`,
      });
    }
  }
}

const isProp = (t: ExprType): t is PropertyType => t !== 'list' && t !== 'nil';

function typeError(kind: string, t: ExprType): never {
  return fail('seq.query.form.type', `A form of kind ${kind} is given a ${t}`);
}

const bool = (x: ExprValue, f: (b: boolean) => boolean) => (typeof x === 'boolean' ? f(x) : undefined);

const Comparisons: Readonly<Record<string, (a: number, b: number) => boolean>> = {
  '<': (a, b) => a < b,
  '<=': (a, b) => a <= b,
  '>': (a, b) => a > b,
  '>=': (a, b) => a >= b,
};

function both(x: ExprValue, y: ExprValue, cmp: (a: number, b: number) => boolean): boolean | undefined {
  return typeof x === 'number' && typeof y === 'number' ? cmp(x, y) : undefined;
}

/** Numeric operands, a literal read as a number. */
function numeric(kind: string, ops: readonly Compiled[]): Compiled[] {
  return ops.map((o) => {
    if (isNumeric(o.type)) return o;
    if (o.literal !== undefined && !Array.isArray(o.literal)) {
      const n = coerce(o.literal as Atom, 'number');
      if (typeof n === 'number') return { type: 'number', level: -1, literal: n, ev: () => n };
    }
    return typeError(kind, o.type);
  });
}

function member(x: ExprValue, set: ReadonlySet<Atom>): boolean | undefined {
  if (x === undefined) return undefined;
  if (Array.isArray(x)) return x.some((w) => set.has(w));
  return set.has(x as Atom);
}

/** `=` and `!=`: one type each side, a literal read as the other side's; on a set of words, membership. */
function equality(kind: string, a: Compiled, b: Compiled, level: number): Compiled {
  let [x, y] = [a, b];
  if (x.literal !== undefined && y.literal === undefined) [x, y] = [y, x];
  let lit: Atom | undefined;
  if (y.literal !== undefined) {
    if (Array.isArray(y.literal)) return typeError(kind, 'list');
    lit = isProp(x.type) ? coerce(y.literal as Atom, x.type) : (y.literal as Atom);
    if (lit === undefined) return typeError(kind, y.type);
  } else if (
    !(x.type === y.type || (isNumeric(x.type) && isNumeric(y.type)) || x.type === 'words' || y.type === 'words')
  )
    return fail('seq.query.form.type', `A form of kind ${kind} compares a ${x.type} and a ${y.type}`);
  const neg = kind === '!=';
  const eq = (p: ExprValue, q: ExprValue): boolean | undefined => {
    if (p === undefined || q === undefined) return undefined;
    const ps = Array.isArray(p) ? p : [p];
    const qs = Array.isArray(q) ? q : [q];
    const hit = ps.some((u) => qs.some((w) => u === w));
    return neg ? !hit : hit;
  };
  if (lit !== undefined) {
    const l = lit;
    return { type: 'boolean', level, ev: (rs, r, a2) => eq(x.ev(rs, r, a2), l) };
  }
  return { type: 'boolean', level, ev: (rs, r, a2) => eq(x.ev(rs, r, a2), y.ev(rs, r, a2)) };
}

function propOf(cx: CompileContext, form: Node): PropertyDef {
  arity(form, 1);
  const name = args(form)[0];
  if (typeof name !== 'string') return fail('seq.query.form.param', 'A form of kind prop takes a property’s name');
  const f = findProperty(cx.properties, name);
  if ('def' in f) return f.def;
  if ('missing' in f)
    return fail('seq.query.extension.missing', `No extension provides the property ${name}`, {
      hint: `namespace: ${f.missing}`,
    });
  return fail('seq.query.form.param', `No property is named ${name}`, { hint: f.unknown });
}

/** `(value <values annotation> :field f)`: a values channel on the element's residues, as a residue property. */
function valueOf(form: Node, values: Map<string, PropertyDef>): PropertyDef {
  arity(form, 1);
  const at = args(form)[0];
  const field = form.params?.field;
  if (typeof at !== 'string' || !at.startsWith('/'))
    return fail('seq.query.form.param', 'A form of kind value takes a values annotation’s address');
  if (field !== undefined && typeof field !== 'string')
    return fail('seq.query.form.param', 'A form of kind value takes a word for :field');
  for (const k of Object.keys(form.params ?? {}))
    if (k !== 'field') fail('seq.query.form.param', 'A form of kind value takes :field only');
  const key = `${at} ${field ?? ''}`;
  let d = values.get(key);
  if (!d) {
    d = {
      name: `(value ${at})`,
      aliases: [],
      level: 'residue',
      type: 'number',
      read: { values: at, field: field ?? '' },
    };
    values.set(key, d);
  }
  return d;
}

/**
 * A test compiled: boolean, else `seq.query.form.type`; reading no level finer than `maxLevel` when given, else
 * `seq.query.property.level`.
 */
export function compilePredicate(
  cx: CompileContext,
  test: Value,
  absent: Absent,
  maxLevel?: number,
  what = 'A test',
): Predicate {
  const uses = new Set<PropertyDef>();
  const c = compile(cx, test, uses, absent, new Map());
  if (c.type !== 'boolean') typeError('test', c.type);
  if (maxLevel !== undefined && c.level > maxLevel)
    fail('seq.query.property.level', `${what} reads a property finer than its elements`, {
      hint: `it reads ${[...uses].map((d) => `${d.name} (${d.level})`).join(', ')}`,
    });
  return { level: c.level, uses: [...uses], test: (rs, r, a) => c.ev(rs, r, a) === true };
}
