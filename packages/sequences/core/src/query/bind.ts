/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The binder (seqspec:seqql#binding-and-checking): read nodes to the JSON form, short names resolved by operand kind
// and slot chain, positionals named, defaults omitted; and macro expansion for checking and compiling.
// The printer resolves with the same functions, so what it prints binds back to what it printed.
import { queryIssue, suggest } from "./issue.ts";
import type { Language } from "./language.ts";
import { MAX_DEPTH, type QueryNode } from "./read.ts";
import { slotMessage } from "./table.ts";
import type { ArgDef, Expression, MacroReport, QueryIssue, Slot, SymbolDef, Type } from "./types.ts";
import { isApply } from "./types.ts";

/** What an operand is, for the binder; a word or string is a ref, so a selection. */
export type OperandKind = "selection" | "test" | "value" | "unknown";

/** What a call gives the resolver: its positionals' kinds and its keywords' names. */
export interface CallShape {
  readonly kinds: readonly OperandKind[];
  readonly keywords: readonly string[];
}

type Problem = Required<Pick<QueryIssue, "code" | "message" | "signature">> & Pick<QueryIssue, "known">;
/** The slots a test reads where it sits; null outside a `:where` or inside a selection. */
type Chain = readonly Slot[] | null;

/** The slot chains, by the slot a binder's `:where` binds: the slots its test reads, in resolution order. */
export const SLOT_CHAINS: Readonly<Record<Slot, readonly Slot[]>> = {
  position: ["position", "system"],
  system: ["system"],
  item: ["item"],
  group: ["group", "item", "system"],
};

const kindOfType = (t: Type): OperandKind =>
  t === "query" ? "selection" : t === "bool" || t === "test" ? "test" : "value";
export const kindOfValue = (v: string | number | boolean): OperandKind =>
  typeof v === "string" ? "selection" : typeof v === "boolean" ? "test" : "value";
const takes = (s: SymbolDef, shape: CallShape) =>
  shape.keywords.every((k) => s.args.some((a) => !a.positional && a.name === k)) &&
  (shape.kinds.length <= s.args.filter((a) => a.positional).length || s.args.some((a) => a.rest));

/** The chain an argument's value is bound in: a `:where` its binder's; a selection none (it is computed once). */
export const argChain = (a: ArgDef, chain: Chain) =>
  a.binds ? SLOT_CHAINS[a.binds] : a.type === "query" ? null : chain;

/** What operands are together: all selections, all tests or values, a mix, or nothing known. */
function together(kinds: readonly OperandKind[]): OperandKind | "mixed" {
  const known = kinds.filter((k) => k !== "unknown");
  const selections = known.filter((k) => k === "selection").length;
  if (known.length === 0) return "unknown";
  return selections === known.length ? "selection" : selections === 0 ? "test" : "mixed";
}

/** The kind of a call written with this name: its symbol's result, or for and, or, not its operands'. */
export function calledKind(language: Language, written: string, operands: () => OperandKind[]): OperandKind {
  const full = language.symbols.get(written);
  if (full) return kindOfType(full.result);
  const kinds = new Set((language.aliases.get(written) ?? []).map((s) => kindOfType(s.result)));
  if (kinds.size < 2) return kinds.values().next().value ?? "unknown";
  const kind = together(operands());
  return kind === "mixed" ? "unknown" : kind;
}

/** A full name is its symbol; a short name by its candidates, operand kind, then the innermost slot chain. */
export function resolveName(language: Language, written: string, chain: Chain, shape: CallShape): SymbolDef | Problem {
  const full = language.symbols.get(written);
  if (full) return full;
  const candidates = language.aliases.get(written);
  if (!candidates) {
    const names = [...language.aliases.keys(), ...(written.includes(".") ? language.symbols.keys() : [])];
    const known = suggest(written, names);
    const best = known[0] && (language.aliases.get(known[0])?.[0] ?? language.symbols.get(known[0]));
    const message = `unknown symbol "${written}"${known[0] ? `; did you mean "${known[0]}"?` : ": no symbol has a similar name"}`;
    return { code: "seqql.symbol.unknown", message, signature: best ? best.signature : null, known };
  }
  if (candidates.every((s) => s.kind === "property")) {
    const fits = candidates.filter((s) => takes(s, shape));
    for (const slot of chain ?? []) {
      const s = fits.find((f) => f.reads === slot);
      if (s) return s;
    }
    const shown = fits.length > 0 ? fits : candidates;
    const message = slotMessage(
      written,
      shown.map((s) => s.reads as Slot),
      !chain && shown.length > 1,
    );
    return { code: "seqql.slot", message, signature: [...new Set(shown.map((s) => s.signature))].join(" or ") };
  }
  if (candidates.length === 1) return candidates[0] as SymbolDef;
  const kind = together(shape.kinds);
  if (kind === "mixed") {
    const message = `${written} mixes a selection and a test: filter with (pick Q :where T)`;
    return { code: "seqql.type", message, signature: candidates.map((s) => s.signature).join(" or ") };
  }
  const want = kind === "unknown" ? (chain ? "test" : "selection") : kind;
  return candidates.find((s) => kindOfType(s.result) === want) ?? (candidates[0] as SymbolDef);
}

/** `(field label)` and the like bind as their own symbols: the property it names, or null. */
export function fieldRewrite(language: Language, s: SymbolDef, args: Readonly<Record<string, unknown>> = {}) {
  const name = args["0"];
  if (s.name !== `seq.${s.reads}.field` || typeof name !== "string" || Object.keys(args).length !== 1) return null;
  const own = language.aliases.get(name) ?? [];
  return own.find((t) => t.kind === "property" && t.reads === s.reads && !t.args.some((a) => a.positional)) ?? null;
}

/** A list's arguments: positionals, `:keyword value` pairs (no value: the next item is a keyword or none), and the
 * first positional written after a keyword. */
function split(items: readonly QueryNode[]) {
  const positionals: QueryNode[] = [];
  const keywords: [QueryNode & { kind: "keyword" }, QueryNode | undefined][] = [];
  let late: QueryNode | null = null;
  for (let k = 1; k < items.length; k++) {
    const n = items[k] as QueryNode;
    if (n.kind === "keyword") {
      const next = items[k + 1];
      keywords.push([n, next && next.kind !== "keyword" ? next : undefined]);
      if (next && next.kind !== "keyword") k++;
    } else {
      if (keywords.length > 0) late ??= n;
      positionals.push(n);
    }
  }
  return { positionals, keywords, late };
}

export interface QueryBinding {
  /** The JSON form, full names; null when an error was found. */
  readonly expr: Expression | null;
  readonly issues: readonly QueryIssue[];
  /** The node each JSON path was bound from: positions for later issues, a string's raw text for macros. */
  readonly source: ReadonlyMap<string, QueryNode>;
}

const hintOr = (s: SymbolDef, code: string, text: string) =>
  s.hints.find((h) => h.code === code && h.value === undefined)?.text ?? text;

/**
 * A label written as a number keeps its written form where the number would lose it (CM-01): IMGT's label `111.10` is
 * not `111.1`, nor `10.0` the label `10`. Only a label: another text argument written as a number is `seqql.type`
 * whatever its form (LA-23: `1.0` passed where `10` was refused). Where positionals and keywords meet the signature,
 * the one place that knows the argument's type.
 */
const asWritten = (def: ArgDef, n: QueryNode, v: Expression | null): Expression | null =>
  typeof v === "number" && def.type === "label" && n.kind === "atom" && n.raw !== String(v) ? n.raw : v;

/** Binds a read tree with a language's symbols. Every problem is an issue with its path, position and the
 * symbol's signature; the binder goes on past one, so siblings are checked too. */
export function bindQuery(node: QueryNode, language: Language): QueryBinding {
  const issues: QueryIssue[] = [];
  const source = new Map<string, QueryNode>();
  const report = (code: string, message: string, n: QueryNode, path: string, more: Partial<QueryIssue> = {}) =>
    issues.push(queryIssue(code, message, { path, at: n.at, ...more }));
  const kindOf = (n: QueryNode): OperandKind =>
    n.kind === "atom"
      ? kindOfValue(n.value)
      : n.kind === "keyword"
        ? "unknown"
        : calledKind(language, headOf(n), () => split(n.items).positionals.map(kindOf));

  function bind(n: QueryNode, path: string, chain: Chain): Expression | null {
    source.set(path, n);
    if (n.kind !== "list") return n.kind === "atom" ? n.value : null;
    const written = headOf(n);
    const { positionals, keywords, late } = split(n.items);
    const shape = { kinds: positionals.map(kindOf), keywords: keywords.map(([k]) => k.name) };
    const s = resolveName(language, written, chain, shape);
    if (typeof s !== "function") {
      report(s.code, s.message, n.items[0] as QueryNode, path, s);
      return null;
    }
    const before = issues.length;
    const say = (code: string, message: string, at: QueryNode, known?: string[]) =>
      report(code, message, at, path, { signature: s.signature, ...(known ? { known } : {}) });
    const order = `positional arguments come before keywords: ${s.signature}`;
    if (late) say("seqql.arg.order", hintOr(s, "seqql.arg.order", order), late);

    const args: Record<string, Expression | Expression[]> = {};
    const defs = s.args.filter((a) => a.positional);
    for (let k = 0; k < positionals.length; k++) {
      const p = positionals[k] as QueryNode;
      const at = Math.min(k, defs.length - 1);
      const def = defs[at];
      if (!def || (k > at && !def.rest)) {
        const count = `${defs.length} positional argument${defs.length === 1 ? "" : "s"}`;
        say("seqql.arg.count", hintOr(s, "seqql.arg.count", `${written} takes ${count}: ${s.signature}`), p);
        break;
      }
      const v = asWritten(
        def,
        p,
        bind(p, def.rest ? `${path}/args/${at}/${k - at}` : `${path}/args/${at}`, argChain(def, chain)),
      );
      const list = args[at];
      if (v !== null && def.rest) args[at] = Array.isArray(list) ? [...list, v] : [v];
      else if (v !== null && v !== def.default) args[at] = v;
    }
    const given = new Map<string, Expression | null>();
    for (const [kw, value] of keywords) {
      const def = s.args.find((a) => !a.positional && a.name === kw.name);
      const names = s.args.filter((a) => !a.positional).map((a) => `:${a.name}`);
      const known = suggest(`:${kw.name}`, names);
      const takes = names.length > 0 ? `it takes ${names.join(", ")}` : "it takes no keywords";
      const fix = known[0] ? `did you mean ${known[0]}?` : hintOr(s, "seqql.arg.unknown", takes);
      if (!def) say("seqql.arg.unknown", `${written} has no :${kw.name}; ${fix}`, kw, known);
      else if (given.has(kw.name)) say("seqql.arg.repeated", `:${kw.name} is given twice`, kw);
      else if (!value) say("seqql.arg.missing", `:${kw.name} needs a value: ${s.signature}`, kw);
      else given.set(kw.name, asWritten(def, value, bind(value, `${path}/args/${kw.name}`, argChain(def, chain))));
    }
    if (issues.length > before) return null;
    for (const a of s.args) {
      const v = given.get(a.name);
      if (!a.positional && v != null && v !== a.default) args[a.name] = v;
    }
    const own = written === s.name ? null : fieldRewrite(language, s, args);
    if (own) return { head: own.name };
    const out = Object.keys(args).length > 0 ? { head: s.name, args } : { head: s.name };
    return s.normalize ? s.normalize(out) : out;
  }

  const expr = bind(node, "", null);
  return { expr: issues.length > 0 ? null : expr, issues, source };
}

const headOf = (n: QueryNode & { kind: "list" }) => String((n.items[0] as QueryNode & { kind: "atom" }).value);

/** A decoded index into a string atom as a column offset from its start: escapes count as written. */
function columnOf(raw: string, index: number): number {
  if (!raw.startsWith('"')) return index;
  let r = 1;
  for (let k = 0; k < index && r < raw.length - 1; k++) r += raw[r] !== "\\" ? 1 : raw[r + 1] === "u" ? 6 : 2;
  return r;
}

/**
 * Replaces every macro call by its expansion, recursively, for checking and compiling. A macro's reports
 * become issues at the call's path; with the binding's `source`, positions point into the string argument. `base` is
 * the expression's own path in the bound one, when it is a part of it.
 */
export function expandMacros(expr: Expression, language: Language, source?: QueryBinding["source"], base = "") {
  const issues: QueryIssue[] = [];
  const place = (path: string, r: MacroReport, index: number | undefined) => {
    const n = source?.get(`${path}/args/${r.arg ?? "0"}`);
    if (!n || index === undefined) return null;
    return { line: n.at.line, column: n.at.column + columnOf(n.kind === "atom" ? n.raw : "", index) };
  };
  function walk(e: Expression, path: string, depth: number): Expression {
    if (!isApply(e)) return e;
    const s = language.symbols.get(e.head);
    if (s?.expand) {
      if (depth >= MAX_DEPTH) {
        issues.push(queryIssue("seqql.macro.depth", `${s.alias ?? s.name} expands into itself`, { path }));
        return e;
      }
      const out = s.expand(e, (r) => {
        const at = place(path, r, r.index) ?? source?.get(path)?.at ?? null;
        const open = place(path, r, r.open);
        const message = open ? `${r.message} at ${open.line}:${open.column}` : r.message;
        issues.push(queryIssue(r.code, message, { path, at, open, signature: s.signature }));
      });
      return walk(out, path, depth + 1);
    }
    if (!e.args) return e;
    const args: Record<string, Expression | Expression[]> = {};
    for (const [k, v] of Object.entries(e.args)) {
      args[k] = Array.isArray(v)
        ? v.map((x: Expression, j) => walk(x, `${path}/args/${k}/${j}`, depth))
        : walk(v as Expression, `${path}/args/${k}`, depth);
    }
    return { head: e.head, args };
  }
  return { expr: walk(expr, base, 0), issues };
}
