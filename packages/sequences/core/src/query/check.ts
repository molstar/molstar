/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The checker (seqspec:seqql#types, seqspec:seqql#errors): types, arity, slots and (current) on a bound expression, or
// on stored JSON. Macros are expanded as it goes: a macro's own arguments are checked first, then its expansion, whose
// issues are placed at the call. The result is the expansion the runtime compiles.
import { argChain, bindQuery, expandMacros, type QueryBinding } from "./bind.ts";
import { queryIssue, suggest } from "./issue.ts";
import { createLanguage, type Language } from "./language.ts";
import { readQuery, type TextPosition } from "./read.ts";
import { LABEL_RANGE } from "./sel.ts";
import { slotMessage } from "./table.ts";
import type { ArgDef, Expression, QueryIssue, Slot, SymbolDef, Type } from "./types.ts";
import { expressionTooDeep, isApply, MAX_EXPRESSION_DEPTH } from "./types.ts";

const BASE = createLanguage().language;
const WANT: Record<string, string> = {
  query: "a selection",
  test: "a test, as (> (length) 50)",
  bool: "true or false",
  num: "a number",
  int: "a whole number from 0",
  str: "text",
  json: "a value",
  ref: "a ref",
  label: "labels or label ranges: 10-20, 45A, -5-10",
};
const GOT: Record<string, string> = { query: "a selection", bool: "a test", str: "text", json: "a value" };
/** The results a type takes from a call; `ref`, `label` and an enum take written values only. */
const TAKES: Record<string, readonly string[]> = {
  query: ["query"],
  test: ["bool", "json"],
  bool: ["bool", "json"],
  num: ["num", "int", "json"],
  int: ["int", "json"],
  str: ["str", "json"],
  json: ["str", "num", "int", "bool", "json"],
};
const LITERAL: Record<string, string> = {
  query: "string",
  str: "string",
  ref: "string",
  test: "boolean",
  num: "number",
};
/** A system's own fields, read without `:of` or `:ns`. */
const SYSTEM_FIELDS = ["molecule", "complex", "alignment", "row"];

const fits = (t: Type, v: string | number | boolean) =>
  typeof t === "object"
    ? t.enum.includes(v as string)
    : t === "int"
      ? Number.isInteger(v) && (v as number) >= 0
      : t === "label"
        ? typeof v !== "boolean" && LABEL_RANGE.test(String(v))
        : t === "json" || typeof v === (LITERAL[t] ?? "boolean");
const hintOf = (s: SymbolDef, code: string) => s.hints.find((h) => h.code === code && h.value === undefined)?.text;

export interface QueryCheck {
  /** The expression with its macros expanded, which the runtime compiles; null when an error was found. */
  readonly expr: Expression | null;
  readonly issues: readonly QueryIssue[];
}

/**
 * Checks an expression against a language: argument types, required and surplus arguments, `slice`'s pairs,
 * the slot each property reads, `(current)` outside `each`, a system's own fields. With the binding's `source`,
 * issues point into the text. A hint of the symbol's for the code is the message, except about a selection argument.
 */
export function checkQuery(expr: Expression, language: Language = BASE, source?: QueryBinding["source"]): QueryCheck {
  const issues: QueryIssue[] = [];
  let errors = 0;
  // A stored expression never passed the reader's bound: refused before anything recurses into it (CM-03).
  if (expressionTooDeep(expr)) {
    const message = `the expression nests deeper than ${MAX_EXPRESSION_DEPTH}; write it flatter`;
    return { expr: null, issues: [queryIssue("seqql.read.depth", message)] };
  }

  function walk(e: Expression, path: string, chain: readonly Slot[] | null, inEach: boolean, macro: string | null) {
    if (!isApply(e)) return;
    const s = language.symbols.get(e.head);
    const say = (code: string, message: string, p = path, more: Partial<QueryIssue> = {}) => {
      const where = macro ?? p;
      const place = { path: where, at: source?.get(where)?.at ?? null, signature: s?.signature ?? null };
      const issue = { ...queryIssue(code, message, place), ...more };
      issues.push(issue);
      if (issue.severity === "error") errors++;
    };
    if (!s) {
      const known = suggest(e.head, language.symbols.keys());
      return say(
        "seqql.symbol.unknown",
        `unknown symbol "${e.head}"${known[0] ? `; did you mean "${known[0]}"?` : ""}`,
      );
    }
    const before = errors;
    const name = s.alias ?? s.name;
    const args = e.args ?? {};
    const positional = s.args.filter((a) => a.positional);
    const label = (a: ArgDef) => (a.positional ? `${name}'s ${a.name}` : `:${a.name}`);
    if (s.reads && !chain?.includes(s.reads)) {
      say("seqql.slot", hintOf(s, "seqql.slot") ?? slotMessage(s.name, [s.reads]));
    }
    if (s.name === "seq.generator.current" && !inEach) {
      say("seqql.current", hintOf(s, "seqql.current") ?? "(current) is only meaningful inside each");
    }

    const check = (a: ArgDef, v: Expression, p: string) => {
      const result = isApply(v) ? language.symbols.get(v.head)?.result : undefined;
      const ok = isApply(v)
        ? result === undefined || (typeof a.type === "string" && (TAKES[a.type] ?? []).includes(result as string))
        : fits(a.type, v);
      if (!ok) {
        const known = typeof a.type === "object" && typeof v === "string" ? suggest(v, a.type.enum) : [];
        const want = typeof a.type === "object" ? `one of ${a.type.enum.join(", ")}` : WANT[a.type];
        const head = isApply(v) && `(${language.symbols.get(v.head)?.alias ?? v.head}${v.args ? " …" : ""})`;
        const quote = typeof v === "number" && (a.type === "query" || a.type === "ref");
        const got = head ? `${GOT[result as string] ?? "a number"} ${head}` : JSON.stringify(v);
        const fix = quote ? `; a ref that reads as a number is quoted: "${v}"` : "";
        // Capitals where the values are lower-case words (`:side A`, QE1a) name something else, a chain: the
        // argument's doc says what it takes, where a did-you-mean would steer to the wrong value.
        const named = typeof a.type === "object" && typeof v === "string" && /[A-Z]/.test(v);
        const mean = named ? ` (the values are lower-case): ${a.doc}` : known[0] ? `; did you mean "${known[0]}"?` : "";
        const word = typeof v === "string" ? v.toLowerCase().replace(/s$/, "") : "";
        const valued =
          typeof a.type === "object" ? s.hints.find((h) => h.code === "seqql.type" && h.value === word) : undefined;
        const hint = valued?.text ?? (a.type !== "query" && a.type !== "test" ? hintOf(s, "seqql.type") : undefined);
        say("seqql.type", hint ?? `${label(a)} takes ${want}, not ${got}${mean}${fix}`, p, known[0] ? { known } : {});
      }
      const each = inEach || (s.name === "seq.modifier.each" && a === positional[1]);
      walk(v, p, argChain(a, chain), each, macro);
    };
    for (const [k, v] of Object.entries(args)) {
      const numbered = /^\d+$/.test(k);
      const a = numbered ? positional[Number(k)] : s.args.find((x) => !x.positional && x.name === k);
      const p = `${path}/args/${k}`;
      if (!a) say(numbered ? "seqql.arg.count" : "seqql.arg.unknown", `${name} has no argument ${k}`, p);
      else if (Array.isArray(v) !== Boolean(a.rest)) say("seqql.type", `${label(a)} is one value, or a rest a list`, p);
      else if (Array.isArray(v)) for (const [j, x] of (v as Expression[]).entries()) check(a, x, `${p}/${j}`);
      else check(a, v as Expression, p);
    }
    for (const a of s.args) {
      const key = a.positional ? String(positional.indexOf(a)) : a.name;
      if (a.optional || a.default !== undefined || key in args) continue;
      const hint = a.type !== "query" ? hintOf(s, "seqql.arg.missing") : undefined;
      say("seqql.arg.missing", hint ?? `${name} needs ${a.positional ? a.doc : `:${a.name}`}: ${s.signature}`);
    }

    // A text property (a label, an id) compared with a number checks clean and never matches (CM-06): an error that
    // says to quote it.
    if (s.name === "core.rel.eq" || s.name === "core.rel.neq" || s.name === "core.rel.in") {
      const list = (k: string) =>
        Array.isArray(args[k]) ? (args[k] as Expression[]) : args[k] === undefined ? [] : [args[k] as Expression];
      const sides = [
        ...list("0").map((v) => [v, `${path}/args/0`] as const),
        ...list("1").map((v, j) => [v, Array.isArray(args["1"]) ? `${path}/args/1/${j}` : `${path}/args/1`] as const),
      ];
      const text = sides.find(([v]) => isApply(v) && language.symbols.get(v.head)?.result === "str")?.[0];
      if (text && isApply(text)) {
        const head = `(${language.symbols.get(text.head)?.alias ?? text.head})`;
        for (const [v, p] of sides)
          if (typeof v === "number") say("seqql.type", `${head} is text: quote the number to compare, "${v}"`, p);
      }
    }

    const bounds = args["1"];
    if (s.name === "seq.modifier.slice" && Array.isArray(bounds) && bounds.every((b) => typeof b === "number")) {
      const last = `${path}/args/1/${bounds.length - 1}`;
      if (bounds.length % 2) say("seqql.arg.count", hintOf(s, "seqql.arg.count") ?? "", last);
      for (let j = 0; j + 1 < bounds.length; j += 2) {
        const [start, stop] = [bounds[j] as number, bounds[j + 1] as number];
        // A start after the stop wraps a cyclic system's origin, known at run time (CM-10); equal ends select nothing.
        if (stop !== start) continue;
        const message = `slice ${start} ${stop} selects nothing: the stop comes after the start, as (slice Q 9 20)`;
        say("seqql.slice.range", message, `${path}/args/1/${j}`);
      }
    }
    const field = args["0"];
    if (typeof field === "string" && s.name === "seq.system.field" && !args.of && !args.ns) {
      const known = suggest(field, SYSTEM_FIELDS);
      if (!SYSTEM_FIELDS.includes(field)) {
        const fix = known[0] ? `did you mean "${known[0]}"?` : hintOf(s, "seqql.field.unknown");
        say("seqql.field.unknown", `a system has no field "${field}"; ${fix}`, `${path}/args/0`, { known });
      }
    } else if (typeof field === "string" && s.name === "seq.item.field") {
      const own = [...language.aliases].flatMap(([n, c]) => (c.some((t) => t.reads === "item") ? [n] : []));
      // a typo by swap or substitution: a column may well be a property's name and one more letter ("evalue")
      const known = suggest(field, own).filter((n) => n.length === field.length);
      const message = `an item has no field "${field}"; did you mean "${known[0]}"?`;
      if (known[0]) say("seqql.field.near", message, `${path}/args/0`, { known });
    }

    if (s.expand && errors === before) {
      const x = expandMacros(e, language, source, path);
      issues.push(...x.issues);
      errors += x.issues.length;
      if (x.issues.length === 0) walk(x.expr, path, chain, inEach, macro ?? path);
    }
  }

  walk(expr, "", null, false, null);
  return { expr: errors > 0 ? null : expandMacros(expr, language).expr, issues };
}

/** Reads, binds and checks a text: the stored JSON (macros as written), its expansion (null on an error) and every
 * issue, the reader's, the binder's or the checker's. `start` places the text in a larger one: a reply's block
 * (seqspec:seqql#queries-in-text). */
export function checkText(text: string, language: Language = BASE, start?: TextPosition) {
  const read = readQuery(text, start);
  const bound = read.node ? bindQuery(read.node, language) : null;
  const checked = bound?.expr != null ? checkQuery(bound.expr, language, bound.source) : null;
  const issues = [...read.issues, ...(bound?.issues ?? []), ...(checked?.issues ?? [])];
  return { expr: bound?.expr ?? null, expanded: checked?.expr ?? null, issues };
}
