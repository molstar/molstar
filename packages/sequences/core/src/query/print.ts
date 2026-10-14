/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The canonical printer (seqspec:seqql#the-printer): JSON to the S-expression documents store beside it. A short name
// wherever it binds back to the same symbol at that place (resolved with the binder's own functions, from what is
// printed), positionals then keywords in signature order, defaults omitted, one line when it fits in 80 columns.
// `printQuery(bindQuery(readQuery(t)))` is a fixed point, and binding a print gives the JSON back.
import { argChain, calledKind, fieldRewrite, kindOfValue, type OperandKind, resolveName } from "./bind.ts";
import { createLanguage, type Language } from "./language.ts";
import { NUMBER } from "./read.ts";
import { type Expression, isApply, type Slot } from "./types.ts";

const WIDTH = 80;
/** A word the reader reads back as the same string (seqspec:seqql#the-reader): no delimiter, no leading ":" or "$". */
const WORD = /^[^ \t\r\n()";[\],:$][^ \t\r\n()";[\],]*$/;
const BASE = createLanguage().language;

type Doc = string | { readonly head: string; readonly parts: readonly (readonly [string | null, Doc])[] };

const atomText = (v: string | number | boolean, quote: boolean) =>
  typeof v !== "string"
    ? String(v)
    : !quote && WORD.test(v) && !NUMBER.test(v) && v !== "true" && v !== "false"
      ? v
      : JSON.stringify(v);

const list = (v: Expression | readonly Expression[] | undefined): readonly Expression[] =>
  v === undefined ? [] : Array.isArray(v) ? v : [v as Expression];

/**
 * The document of an expression, and the operand kind the binder will read from it (seqspec:seqql#binding-and-checking
 * step 2).
 */
function build(language: Language, e: Expression, chain: readonly Slot[] | null, quote = false): [Doc, OperandKind] {
  if (!isApply(e)) return [atomText(e, quote), kindOfValue(e)];
  const s = language.symbols.get(e.head);
  const args = e.args ?? {};
  const defs = s?.args ?? [];
  const positional = defs.filter((a) => a.positional);
  const parts: [string | null, Doc][] = [];
  const kinds: OperandKind[] = [];
  const keywords: string[] = [];
  const numbered = Object.keys(args).filter((k) => /^\d+$/.test(k));
  for (const k of numbered.sort((a, b) => Number(a) - Number(b))) {
    const def = positional[Number(k)];
    if (def && !def.rest && args[k] === def.default) continue;
    for (const v of list(args[k])) {
      const [doc, kind] = build(language, v, def ? argChain(def, chain) : chain, s?.name === "seq.macro.sel");
      parts.push([null, doc]);
      kinds.push(kind);
    }
  }
  const named = defs.filter((a) => !a.positional && a.name in args).map((a) => a.name);
  for (const k of [...named, ...Object.keys(args).filter((k) => !numbered.includes(k) && !named.includes(k))]) {
    const def = defs.find((a) => !a.positional && a.name === k);
    if (def && args[k] === def.default) continue;
    keywords.push(k);
    for (const v of list(args[k])) parts.push([k, build(language, v, def ? argChain(def, chain) : chain)[0]]);
  }
  if (!s) return [{ head: e.head, parts }, "unknown"];
  const short =
    s.alias !== null &&
    !fieldRewrite(language, s, e.args) &&
    resolveName(language, s.alias, chain, { kinds, keywords }) === s;
  const head = short ? (s.alias as string) : s.name;
  return [{ head, parts }, calledKind(language, head, () => kinds)];
}

const flat = (d: Doc): string =>
  typeof d === "string"
    ? d
    : `(${d.head}${d.parts.map(([k, v]) => `${k === null ? " " : ` :${k} `}${flat(v)}`).join("")})`;

/** One line when it fits; else the head and the positionals that fit, then one argument a line, two spaces in. */
function layout(d: Doc, column: number): string {
  const line = flat(d);
  if (typeof d === "string" || column + line.length <= WIDTH) return line;
  let out = `(${d.head}`;
  let k = 0;
  for (; k < d.parts.length; k++) {
    const [kw, v] = d.parts[k] as readonly [string | null, Doc];
    const text = flat(v);
    if (kw !== null || column + out.length + 1 + text.length > WIDTH) break;
    out += ` ${text}`;
  }
  const pad = " ".repeat(column + 2);
  for (; k < d.parts.length; k++) {
    const [kw, v] = d.parts[k] as readonly [string | null, Doc];
    const lead = kw === null ? "" : `:${kw} `;
    out += `\n${pad}${lead}${layout(v, column + 2 + lead.length)}`;
  }
  return `${out})`;
}

/**
 * The canonical S-expression of an expression, short names as the language resolves them (default: the base); `chain`,
 * the slot chain of the `:where` it sits in, for a part of a query printed alone (a property's short name).
 */
export function printQuery(expr: Expression, language: Language = BASE, chain: readonly Slot[] | null = null): string {
  return layout(build(language, expr, chain)[0], 0);
}
