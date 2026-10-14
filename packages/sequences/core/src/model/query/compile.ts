/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The compiler (seqspec:seqql#overview, seqspec:seqql#extension): an expression, checked and with its macros expanded,
// to a query function over a context; a host's env (its language, its runtimes, its letter classes) from registry
// entries, with no module state. A word where a selection goes compiles to a ref.

import {
  type Apply,
  checkQuery,
  createLanguage,
  type Expression,
  isApply,
  type Language,
  type LanguageExtension,
  printQuery,
  type QueryBinding,
  type QueryIssue,
  queryIssue,
  type SymbolDef,
} from "@molstar/seq-core/query";
import { F64, QueryAbort, QueryContext, type QueryOptions } from "./context.ts";
import { MODIFIER_RUNTIMES } from "./modifiers.ts";
import { ONTO_RUNTIMES } from "./onto.ts";
import { PARTNERS_RUNTIMES } from "./partners.ts";
import { PROPERTY_RUNTIMES } from "./properties.ts";
import { BASE_RUNTIMES, refQuery } from "./runtime.ts";
import { SCAN_RUNTIMES } from "./scans.ts";
import { SeqSelection } from "./selection.ts";

export type QueryFn<T> = (ctx: QueryContext) => T;

/** What a runtime's `compile` reads its call through: arguments compiled by the kind it needs, keys as in the JSON. */
export interface CompileContext {
  readonly call: Apply;
  /** The call's JSON pointer, for the issues its runs report. */
  readonly path: string;
  readonly env: QueryEnv;
  /** A selection argument (a word is a ref); `key` may go deeper, as `1/0/args/0`. Undefined when absent. */
  query(key: string): QueryFn<SeqSelection> | undefined;
  /** Each selection of a rest argument. */
  queries(key: string): QueryFn<SeqSelection>[];
  /** A number argument that writes its value to `F64[0]`, never boxed: a numeric constant, or a call whose
   * runtime has `number`; undefined for anything else (then read it with `value`). */
  number(key: string): QueryFn<void> | undefined;
  /** A value or test argument: a constant, or a call. */
  value(key: string): QueryFn<unknown> | undefined;
  values(key: string): QueryFn<unknown>[];
  /** A written argument as it is in the JSON (a ref, a label, an enum, an index). */
  raw<T = Expression>(key: string): T | undefined;
  /** A text argument compared with constant strings: a call whose runtime has `equals`; undefined for anything else. */
  equals(key: string, values: readonly string[]): QueryFn<number> | undefined;
}

export interface SymbolRuntime {
  /** The definition object it runs: a runtime is added with its symbol. */
  readonly symbol: SymbolDef;
  compile(c: CompileContext): QueryFn<unknown>;
  /**
   * A property of numbers: its value written to `F64[0]` (NaN: missing), for the relations (seqspec:seqql#performance:
   * no boxing).
   */
  number?(c: CompileContext): QueryFn<void>;
  /**
   * A property of text compared with constant strings, without making its text per element
   * (seqspec:seqql#performance): 1 where its text is one of `values`, 0 where not, -1 where it is missing.
   */
  equals?(c: CompileContext, values: readonly string[]): QueryFn<number>;
}

/** One registry entry: what an extension package ships. */
export interface QueryExtension extends LanguageExtension {
  readonly namespace: string;
  readonly runtimes?: readonly SymbolRuntime[];
  /** Letter-class tables by name, data. */
  readonly classes?: Readonly<Record<string, Readonly<Record<string, string>>>>;
}

export interface QueryEnv {
  readonly language: Language;
  readonly runtimes: ReadonlyMap<string, SymbolRuntime>;
  readonly classes: ReadonlyMap<string, Readonly<Record<string, string>>>;
  /** Registration issues (`createLanguage`'s). */
  readonly issues: readonly QueryIssue[];
}

/** A host's env from its entries: the base language and runtimes plus the extensions'. */
export function createQueryEnv(ext: readonly QueryExtension[] = []): QueryEnv {
  const { language, issues } = createLanguage(ext);
  const runtimes = new Map<string, SymbolRuntime>();
  const base = [
    ...BASE_RUNTIMES,
    ...MODIFIER_RUNTIMES,
    ...ONTO_RUNTIMES,
    ...PARTNERS_RUNTIMES,
    ...SCAN_RUNTIMES,
    ...PROPERTY_RUNTIMES,
  ];
  for (const r of [...base, ...ext.flatMap((e) => e.runtimes ?? [])]) {
    if (language.symbols.get(r.symbol.name) === r.symbol) runtimes.set(r.symbol.name, r);
  }
  return { language, runtimes, classes: new Map(ext.flatMap((e) => Object.entries(e.classes ?? {}))), issues };
}

export interface QueryResult {
  readonly selection: SeqSelection;
  /** The run's warnings, and the error that ended it. */
  readonly issues: readonly QueryIssue[];
}
export interface CompiledQuery {
  /** The checker's and the compiler's; with an error, every run gives an empty selection and these. */
  readonly issues: readonly QueryIssue[];
  run(o: QueryOptions): QueryResult;
}

function compileNode(e: Expression, path: string, env: QueryEnv, query: boolean, number = false): QueryFn<unknown> {
  if (!isApply(e)) return query ? watched(refQuery(String(e), path), path, false) : () => e;
  const r = env.runtimes.get(e.head);
  if (!r) {
    const s = env.language.symbols.get(e.head);
    const message = `${s?.alias ?? e.head} cannot run here: this host has no runtime for ${e.head}`;
    throw new QueryAbort(queryIssue("seqql.symbol.unavailable", message, { path, signature: s?.signature ?? null }));
  }
  const c = contextOf(e, path, env);
  if (number && r.number) return r.number(c);
  const f = r.compile(c);
  return query ? watched(f as QueryFn<SeqSelection>, path, r.symbol.kind === "generator" && !GIVEN.has(e.head)) : f;
}

/** How a runtime reads the call `e` at `path`: its arguments compiled by the kind it asks for. */
function contextOf(e: Apply, path: string, env: QueryEnv): CompileContext {
  const get = (key: string) => key.split("/").reduce<unknown>((x, k) => (x as Record<string, unknown>)?.[k], e.args);
  const at = (key: string) => `${path}/args/${key}`;
  const one = (key: string, q: boolean) => {
    const x = get(key) as Expression | undefined;
    return x === undefined ? undefined : compileNode(x, at(key), env, q);
  };
  const many = (key: string, q: boolean) =>
    ((get(key) ?? []) as Expression[]).map((x, i) => compileNode(x, at(`${key}/${i}`), env, q));
  const c: CompileContext = {
    call: e,
    path,
    env,
    query: (key) => one(key, true) as QueryFn<SeqSelection> | undefined,
    queries: (key) => many(key, true) as QueryFn<SeqSelection>[],
    number: (key) => {
      const x = get(key);
      if (typeof x === "number")
        return () => {
          F64[0] = x;
        };
      if (!isApply(x) || !env.runtimes.get(x.head)?.number) return undefined;
      return compileNode(x, at(key), env, false, true) as QueryFn<void>;
    },
    value: (key) => one(key, false),
    values: (key) => many(key, false),
    raw: <T>(key: string) => get(key) as T | undefined,
    equals: (key, values) => {
      const x = get(key);
      const rt = isApply(x) ? env.runtimes.get(x.head) : undefined;
      return isApply(x) && rt?.equals ? rt.equals(contextOf(x, at(key), env), values) : undefined;
    },
  };
  return c;
}

/** Generators that give nothing on purpose, or what `each` gave. */
const GIVEN = new Set(["seq.generator.empty", "seq.generator.current"]);
/**
 * A selection node, watched for seqspec:seqql#empty-results: the first, in evaluation order, to give nothing from
 * something (a selection inside it that gave groups, or, for a generator, its data) is named if the whole query selects
 * nothing.
 */
function watched(f: QueryFn<SeqSelection>, path: string, source: boolean): QueryFn<SeqSelection> {
  return (ctx) => {
    const mark = ctx.nonEmpty;
    const out = f(ctx);
    if (out.groupCount > 0) {
      ctx.nonEmpty++;
      ctx.lastNonEmpty = out;
    } else if (ctx.emptyAt === null && (source || ctx.nonEmpty > mark)) {
      ctx.emptyAt = path;
      ctx.emptyFrom = ctx.nonEmpty > mark ? ctx.lastNonEmpty : null;
    }
    return out;
  };
}

/** The path of the call an issue's path is in: the call itself, or the call a written argument belongs to. */
function callOf(expr: Expression, path: string): string {
  let [x, at, call] = [expr as unknown, "", ""];
  for (const k of path.split("/").slice(1)) {
    x = (x as Record<string, unknown> | undefined)?.[k];
    at += `/${k}`;
    if (isApply(x)) call = at;
  }
  return call;
}
/** The empty-result note when the query selects nothing and no issue of the call named already says why. */
function diagnose(ctx: QueryContext, out: SeqSelection, expr: Expression): void {
  const at = ctx.emptyAt;
  if (out.groupCount > 0 || at === null) return;
  if (ctx.issues.some((i) => callOf(expr, i.path) === at)) return;
  const why = ctx.why?.get(at);
  const from =
    ctx.emptyFrom && !why ? ` from ${SeqSelection.summary(ctx.emptyFrom, { maxParts: 3, indices: true })}` : "";
  ctx.warn("seqql.empty", at, `${ctx.say(at)} selected nothing${from}${why ? `: ${why}` : ""}`);
}

/**
 * Checks an expression (a bound one, or stored JSON), expands its macros and compiles it. With the binding's
 * `source`, every issue, the runtime's included, points into the text (at the call or the `sel` it came from).
 */
export function compileQuery(expr: Expression, env: QueryEnv, source?: QueryBinding["source"]): CompiledQuery {
  const place = (i: QueryIssue): QueryIssue => {
    if (i.at || !source) return i;
    let p = i.path;
    while (p && !source.has(p)) p = p.slice(0, p.lastIndexOf("/"));
    const at = source.get(p)?.at;
    return at ? { ...i, at: { line: at.line, column: at.column } } : i;
  };
  const checked = checkQuery(expr, env.language, source);
  /**
   * A runtime issue gains its symbol's signature, and the symbol's hint for its code (seqspec:seqql#errors), as checker
   * issues do.
   */
  const help = (i: QueryIssue): QueryIssue => {
    let x: unknown = checked.expr;
    let s = isApply(x) ? env.language.symbols.get(x.head) : undefined;
    for (const k of i.path.split("/").slice(1)) {
      x = (x as Record<string, unknown> | undefined)?.[k];
      if (isApply(x)) s = env.language.symbols.get(x.head) ?? s;
    }
    const hint = s?.hints.find((h) => h.code === i.code && h.value === undefined)?.text;
    return place({ ...i, signature: s?.signature ?? null, message: hint ? `${i.message}; ${hint}` : i.message });
  };
  const errors = checked.issues.filter((i) => i.severity === "error");
  let fn: QueryFn<unknown> | null = null;
  const issues = [...checked.issues];
  if (checked.expr !== null && errors.length === 0) {
    try {
      fn = compileNode(checked.expr, "", env, true);
    } catch (e) {
      if (!(e instanceof QueryAbort)) throw e;
      issues.push(place(e.issue));
    }
  }
  return {
    issues,
    run(o: QueryOptions): QueryResult {
      if (!fn) return { selection: SeqSelection.empty(o.set), issues };
      const ctx = new QueryContext(o);
      ctx.say = (path) => {
        let x: unknown = checked.expr;
        for (const k of path.split("/").slice(1)) x = (x as Record<string, unknown> | undefined)?.[k];
        const reads = isApply(x) ? env.language.symbols.get(x.head)?.reads : undefined;
        return x === undefined ? path : printQuery(x as Expression, env.language, reads ? [reads] : null);
      };
      try {
        ctx.tick();
        const selection = fn(ctx) as SeqSelection;
        diagnose(ctx, selection, checked.expr as Expression);
        return { selection, issues: [...issues, ...ctx.issues.map(help)] };
      } catch (e) {
        if (!(e instanceof QueryAbort)) throw e;
        return { selection: SeqSelection.empty(o.set), issues: [...issues, ...ctx.issues.map(help), help(e.issue)] };
      }
    },
  };
}
