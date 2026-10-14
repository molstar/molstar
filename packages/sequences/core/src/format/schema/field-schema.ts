/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The field layer: our own field types with MVS's names (Mol*'s src/extensions/mvs/tree/
// generic/field-schema.ts), no io-ts. Each type checks a value, reporting issues with a code and a JSON pointer, may
// normalize it (defaults filled), and names itself for the generated reference.

import { Issue } from '@molstar/seq-core/model';
import { issueOf } from '../load/issues.js';
import type { SvsRegistries } from '../registry.js';

/** A typed ref met while checking; the validator resolves it once the tree is indexed. */
export interface RefUse {
  readonly ref: string;
  readonly kinds: readonly string[];
  readonly path: string;
  /** The node whose params hold it (cycles). */
  readonly from: string;
  /** A table column's value: the rows that hold it (one use per distinct value, seqspec:format#validation). */
  readonly count?: number;
}

export interface CheckContext {
  readonly issues: Issue[];
  readonly registries: SvsRegistries;
  readonly refs: RefUse[];
  /** A newer minor version (seqspec:format#versioning): unknown params, kinds, types and themes are warnings. */
  readonly lenient: boolean;
  /** The ref of the node being checked. */
  from: string;
  /** ThemeRef nesting (seqspec:format#themes: at most 4). */
  depth: number;
  /**
   * The document is read from an archive, its main file in this directory of it (seqspec:sources#archives): `uri`s may
   * be relative.
   */
  readonly base?: string;
}

export interface FieldType<V = unknown> {
  /** The type as the reference prints it: `list<int>`, `ref<sequence | ligand>`. */
  readonly name: string;
  check(value: unknown, path: string, ctx: CheckContext): void;
  /** Defaults filled; absent: the value is its own normal form. */
  normalize?(value: V, reg: SvsRegistries): V;
  /** Phantom: the value's TypeScript type. */
  readonly __value?: V;
}

export type ValueOf<T> = T extends FieldType<infer V> ? V : never;

export function isPlainObject(v: unknown): v is Record<string, unknown> {
  return v !== null && typeof v === 'object' && !Array.isArray(v);
}

/** A JSON pointer one step down (RFC 6901 escapes). */
export function ptr(path: string, key: string | number): string {
  return `${path}/${String(key).replaceAll('~', '~0').replaceAll('/', '~1')}`;
}

/** A value for a message, cut short. */
export function show(v: unknown): string {
  const s = v === undefined ? 'nothing' : JSON.stringify(v);
  return s.length > 40 ? `${s.slice(0, 37)}...` : s;
}

/** Reports an issue at `path`, its severity the registered code's. */
export function report(ctx: CheckContext, code: string, path: string, message: string): void {
  ctx.issues.push(issueOf(ctx.registries, code, path, message));
}

/** Something this version does not know: an error, or a warning when the document is a newer minor. */
export function reportUnknown(ctx: CheckContext, code: string, path: string, message: string): void {
  ctx.issues.push(
    ctx.lenient
      ? Issue.warning('svs.version.newer_minor', path, `${message} (a newer minor version: ignored)`)
      : Issue.error(code, path, message),
  );
}

export function typeError(ctx: CheckContext, path: string, expected: string, v: unknown, code = 'svs.params.type') {
  report(ctx, code, path, `expected ${expected}, got ${show(v)}`);
}

/** A context whose issues and refs are kept apart (unions try each alternative). */
export function scratch(ctx: CheckContext): CheckContext {
  return { issues: [], registries: ctx.registries, refs: [], lenient: ctx.lenient, from: ctx.from, depth: ctx.depth };
}

/** A context that reports nowhere (normalizing reads shapes it does not report on). */
export function quiet(registries: SvsRegistries): CheckContext {
  return { issues: [], registries, refs: [], lenient: true, from: '', depth: 0 };
}

export function scalar<V>(name: string, test: (v: unknown) => boolean, code = 'svs.params.type'): FieldType<V> {
  return {
    name,
    check(v, path, ctx) {
      if (!test(v)) typeError(ctx, path, name, v, code);
    },
  };
}

/** A string that `valid` accepts: a type error when not a string, else `code`. */
export function refined(name: string, valid: (s: string) => boolean, code: string, hint: string): FieldType<string> {
  return {
    name,
    check(v, path, ctx) {
      if (typeof v !== 'string') return typeError(ctx, path, name, v);
      if (!valid(v)) report(ctx, code, path, `${show(v)} is not ${hint}`);
    },
  };
}

export const str = scalar<string>('str', (v) => typeof v === 'string');
export const int = scalar<number>('int', (v) => Number.isInteger(v));
export const float = scalar<number>('float', (v) => typeof v === 'number' && Number.isFinite(v));
export const bool = scalar<boolean>('bool', (v) => typeof v === 'boolean');
/** Any JSON value. */
export const json = scalar<unknown>('json', () => true);

export function literal<V extends string | number | boolean>(...values: V[]): FieldType<V> {
  const set = new Set<unknown>(values);
  return scalar(values.map((v) => JSON.stringify(v)).join(' | '), (v) => set.has(v));
}

export function list<V>(item: FieldType<V>): FieldType<V[]> {
  const norm = item.normalize;
  return {
    name: `list<${item.name}>`,
    check(v, path, ctx) {
      if (!Array.isArray(v)) return typeError(ctx, path, this.name, v);
      for (let i = 0; i < v.length; i++) item.check(v[i], ptr(path, i), ctx);
    },
    normalize: norm && ((v, reg) => (Array.isArray(v) ? v.map((x) => norm(x, reg)) : v)),
  };
}

export function tuple<V extends unknown[]>(...items: { [K in keyof V]: FieldType<V[K]> }): FieldType<V> {
  const name = `[${items.map((t) => t.name).join(', ')}]`;
  return {
    name,
    check(v, path, ctx) {
      if (!Array.isArray(v) || v.length !== items.length) return typeError(ctx, path, name, v);
      items.forEach((t, i) => {
        t.check(v[i], ptr(path, i), ctx);
      });
    },
  };
}

export function nullable<V>(type: FieldType<V>): FieldType<V | null> {
  const norm = type.normalize;
  return {
    name: `${type.name} | null`,
    check: (v, path, ctx) => (v === null ? undefined : type.check(v, path, ctx)),
    normalize: norm && ((v, reg) => (v === null ? null : norm(v, reg))),
  };
}

/**
 * The first alternative that accepts the value. When none does, the issues of the alternative that got past the
 * shape test (its issues are not just a type error at `path`) with the fewest, so a bad color says so; else a type
 * error naming every alternative.
 */
export function union<T extends FieldType<unknown>[]>(...types: T): FieldType<ValueOf<T[number]>> {
  const name = types.map((t) => t.name).join(' | ');
  const pick = (v: unknown, ctx: CheckContext) => {
    const tries = types.map((t) => ({ t, c: scratch(ctx) }));
    for (const x of tries) {
      x.t.check(v, '', x.c);
      if (!x.c.issues.some((i) => i.severity === 'error')) return x;
    }
    return null;
  };
  return {
    name,
    check(v, path, ctx) {
      const ok = pick(v, ctx);
      if (ok) return ok.t.check(v, path, ctx);
      let best: CheckContext | null = null;
      for (const t of types) {
        const c = scratch(ctx);
        t.check(v, path, c);
        const shapeOnly = c.issues.every((i) => i.path === path && i.code === 'svs.params.type');
        if (!shapeOnly && (!best || c.issues.length < best.issues.length)) best = c;
      }
      if (best) for (const i of best.issues) ctx.issues.push(i);
      else typeError(ctx, path, name, v);
    },
    normalize(v, reg) {
      const ok = pick(v, quiet(reg));
      return (ok?.t.normalize ? ok.t.normalize(v, reg) : v) as ValueOf<T[number]>;
    },
  };
}

export function dict<V>(value: FieldType<V>): FieldType<Record<string, V>> {
  const norm = value.normalize;
  return {
    name: `dict<str, ${value.name}>`,
    check(v, path, ctx) {
      if (!isPlainObject(v)) return typeError(ctx, path, this.name, v);
      for (const k of Object.keys(v)) value.check(v[k], ptr(path, k), ctx);
    },
    normalize:
      norm &&
      ((v, reg) => (isPlainObject(v) ? Object.fromEntries(Object.entries(v).map(([k, x]) => [k, norm(x, reg)])) : v)),
  };
}

/** The one namespace syntax (seqspec:format#extensions): a short name or reverse-DNS, and `<namespace>:<name>`. */
const NAMESPACE = /^[a-z][a-z0-9_-]*(\.[a-z0-9_-]+)*$/;
const NAMESPACED = /^[a-z][a-z0-9_-]*(\.[a-z0-9_-]+)*:[A-Za-z0-9][A-Za-z0-9_.-]*$/;
/** No extension registers these; `mvs` is the kernel's own table, and the only one a document's `custom` may name. */
export const ReservedNamespaces: ReadonlySet<string> = new Set(['svs', 'seqv', 'mvs']);
export const isNamespace = (ns: string) => NAMESPACE.test(ns);
/** `<namespace>.<code>`: dotted lower-case segments, the namespace being all but (at least) the last. */
export const isIssueCode = (code: string) => /^[a-z][a-z0-9_-]*(\.[a-z0-9_-]+)+$/.test(code);
export const isNamespaced = (name: unknown): name is string => typeof name === 'string' && NAMESPACED.test(name);

const AUTHOR_REF = /^[A-Za-z0-9][A-Za-z0-9_.-]*$/;
const DERIVED_REF = /^[A-Za-z0-9][A-Za-z0-9_.-]*(\/[A-Za-z0-9_.:@-]+)+$/;

/** What a track template's layers write for the node an instance stands for (seqspec:format#placement). */
export const TEMPLATE_NODE = '$node';

export function isRefSyntax(s: string, derived: boolean): boolean {
  return AUTHOR_REF.test(s) || (derived && DERIVED_REF.test(s));
}

/** A ref's value type: a string, carrying (as a phantom) the kinds it may name, for the builder's handles. */
export type RefValue<K extends string = string> = string & { readonly __refKinds?: K };
/** A `ref<kind>` field: the kinds it takes, for the builder (B1). */
export interface RefFieldType<K extends string = string> extends FieldType<RefValue<K>> {
  readonly kinds: readonly K[];
}

/**
 * A ref to a node of one of `kinds`, or a derived ref (`1brs/A`, `msa/3`) a node of those yields
 * (seqspec:format#nodes-and-refs).
 */
export function ref<K extends string>(...kinds: K[]): RefFieldType<K> {
  const name = `ref<${kinds.join(' | ')}>`;
  return {
    name,
    kinds,
    check(v, path, ctx) {
      if (typeof v !== 'string') return typeError(ctx, path, name, v);
      // A template's node: resolved against its track (`svs.ref.syntax` outside a template).
      if (!isRefSyntax(v, true) && v !== TEMPLATE_NODE) {
        return report(ctx, 'svs.ref.syntax', path, `"${v}" is not a ref: [A-Za-z0-9][A-Za-z0-9_.-]*, then /suffixes`);
      }
      ctx.refs.push({ ref: v, kinds, path, from: ctx.from });
    },
  };
}
