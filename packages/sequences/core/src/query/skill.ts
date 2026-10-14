/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The generated skill and reference: the fixed guidance plus a language's table, so a host's skill lists exactly what
// its checker accepts.
import { SLOT_CHAINS } from "./bind.js";
import { SeqqlIssueCodes } from "./codes.js";
import type { Language } from "./language.js";
import type { ArgDef, Slot, SymbolDef, SymbolKind, Type } from "./types.js";

/** Appendix A.1: generic, naming base symbols only; `{refs}`, `{inputs}` and `{table}` are filled in. */
export const GUIDANCE = `SEQUENCE QUERIES (SeqQL)

A query selects positions and annotation items in the sequences. Write one S-expression:
(name positional... :keyword value ...). Words need no quotes; quote "..." text with spaces,
brackets, commas or parentheses, and refs that look like numbers ("123"). ; starts a comment.

Wrap each query in [SeqQL]...[/SeqQL] and send each query in its own block: several blocks are
several queries, and text around them is fine. A query passed as a tool's argument takes no tags.

A word where a selection goes is a ref and selects what it names: a unit's positions, a
molecule's or complex's units, an alignment's columns, a row's residues, an annotation's items.
{refs}
A selection is a list of groups: a sequence, a run of positions, an item, a hit. Filters test
groups; results step through them. Names are short, and a test's name reads what its :where
tests: (label) is an item's label in (items :where ...), a position's in (positions :where ...).

POSITIONS: labels for what a person says, indices for what a tool returned.
Labels are numbered as the ruler shows them (10-20, 45A), inclusive: (range Q 10-20), [10-20].
Indices count from 0 and exclude the end, as tool results give them: (slice Q 9 20), @9:20.

SHORTHAND, anywhere a selection goes: (sel "...")
  R[10-20,45]  labels; R[uniprot: 120-140] in another numbering
  R@9:20  indices 9 to 19
  A[#id1,#id2]  items by id, not label
  $name, $name[10-20]  an input, narrowed
  a, b  union
{inputs}
PATTERNS
Narrow by intersecting: (and Q (letters C)); remove with (and Q (not B)).
Keep groups with (pick Q :where T); test single elements with a generator's :where.
A selection inside a test is computed once: (pick Q :where (overlaps (sel "...")))
Carry positions to another sequence with onto, then narrow or widen there.
(and A B) intersects within one sequence: carry one side over with onto first.

{table}

RESULTS AND ERRORS
A result starts with a summary ("3 groups, 412 positions: 1brs/A 10-20 @9:20 ..."): labels, then
indices. Check it against the request. "Selects nothing" names the part that emptied it and why;
fix that part. An error gives line:column, what is wrong, the symbol's signature and a
suggestion; resend the corrected query.
`;

export interface SkillOptions {
  /** The host's line on what refs name: "Refs: 1brs/A (a chain's sequence), …, ctx.site.hinge (a Context entity)." */
  readonly refs?: string;
  /** The named inputs the host passes: `{ name: "selection", says: "the person's selection" }`. */
  readonly inputs?: readonly { readonly name: string; readonly says: string }[];
}

const SLOTS: readonly Slot[] = ["position", "system", "item", "group"];
const OWNER: Record<Slot, string> = { position: "position's", system: "sequence's", item: "item's", group: "group's" };
const SECTIONS: readonly (readonly [string, readonly SymbolKind[]])[] = [
  ["SELECTIONS", ["generator", "macro"]],
  ["COMBINE", ["combinator", "modifier", "filter"]],
  ["TESTS", ["property", "core"]],
];
const WIDTH = 96;
/** EXAMPLES' first line: models copied the examples' refs (QE1a, q16 and q17). */
const ILLUSTRATIVE = "Refs in the examples are illustrative; use the refs you were given.";
const nameOf = (s: SymbolDef) => s.alias ?? s.name;
const takesPositional = (s: SymbolDef) => s.args.some((a) => a.positional);
const family = (s: SymbolDef) => (/^core\.(rel|logic)\./.exec(s.name)?.[1] ?? "") as "" | "rel" | "logic";

/** Words wrapped greedily at the guidance's width. */
const wrap = (text: string) =>
  text
    .split(" ")
    .reduce<string[]>((lines, w) => {
      const last = lines[lines.length - 1];
      if (last !== undefined && last.length + 1 + w.length <= WIDTH) lines[lines.length - 1] = `${last} ${w}`;
      else lines.push(w);
      return lines;
    }, [])
    .join("\n");

/** TESTS: per slot, its binders, the properties without positionals folded, the slots its chain reads on, then the
 * properties taking positionals one a line; the relations folded, other core symbols one a line, the logic folded. */
function tests(symbols: readonly SymbolDef[], all: readonly SymbolDef[]): string[] {
  const lines: string[] = [];
  for (const slot of SLOTS) {
    const props = symbols.filter((s) => s.kind === "property" && s.reads === slot);
    if (props.length === 0) continue;
    const binders = all.filter((s) => s.args.some((a) => a.binds === slot)).map(nameOf);
    const folded = props.filter((s) => !takesPositional(s)).map((s) => s.signature);
    const chain = SLOT_CHAINS[slot].slice(1).map((t) => `then its ${OWNER[t]}`);
    const tail = chain.length > 0 ? `; ${chain.join(", ")}` : "";
    lines.push(`${slot} (${binders.join(", ")} :where): ${folded.join(" ")}${tail}`);
    for (const s of props.filter(takesPositional)) lines.push(`  ${s.signature}  ${s.doc}`);
  }
  const core = symbols.filter((s) => s.kind === "core");
  const folded = (f: string) => core.filter((s) => family(s) === f).map((s) => s.signature);
  lines.push(`${folded("rel").join(" ")}  a missing value fails every comparison`);
  for (const s of core.filter((s) => family(s) === "")) lines.push(`${s.signature}  ${s.doc}`);
  lines.push(folded("logic").join(" "));
  return lines;
}

/** The agent's skill: the guidance with the host's lines and the language's table, rendered compactly. */
export function skillText(language: Language, opts: SkillOptions = {}): string {
  const symbols = [...language.symbols.values()];
  const examples = symbols
    .flatMap((s) => s.examples.filter((e) => e.skill !== undefined))
    .sort((a, b) => (a.skill as number) - (b.skill as number))
    .map((e) => `${e.text}  ${e.says}`);
  const table = SECTIONS.map(([title, kinds]) => {
    const own = symbols.filter((s) => kinds.includes(s.kind));
    return [title, ...(title === "TESTS" ? tests(own, symbols) : own.map((s) => `${s.signature}  ${s.doc}`))].join(
      "\n",
    );
  });
  const inputs = (opts.inputs ?? []).map((i) => `${i.name} (${i.says})`).join(", ");
  return GUIDANCE.replace("{refs}\n", () => (opts.refs ? `${wrap(opts.refs)}\n` : ""))
    .replace("{inputs}\n", () => (inputs ? `${wrap(`Inputs: ${inputs}.`)}\n\n` : "\n"))
    .replace("{table}", () => [...table, ["EXAMPLES", ILLUSTRATIVE, ...examples].join("\n")].join("\n\n"));
}

const TYPE: Record<string, string> = {
  num: "number",
  int: "whole number",
  str: "text",
  bool: "boolean",
  json: "value",
  ref: "ref",
  label: "label",
  query: "selection",
  test: "test",
};
const typeText = (t: Type, a?: ArgDef) =>
  typeof t === "object"
    ? t.enum.map((v) => `\`${v}\``).join(" \\| ")
    : `${TYPE[t]}${a?.binds ? ` of ${a.binds === "item" ? "an" : "a"} ${a.binds}` : ""}`;

function entry(s: SymbolDef): string {
  const binds = s.args.filter((a) => a.binds).map((a) => ` · \`:${a.name}\` binds: ${a.binds}`);
  const out = [
    `### ${nameOf(s)}`,
    `\`${s.name}\` · ${s.kind} · result: ${typeText(s.result)}${s.reads ? ` · reads: ${s.reads}` : ""}${binds.join("")}`,
    `\`${s.signature}\`: ${s.doc}.`,
  ];
  if (s.args.length > 0) out.push("", "| Argument | Type | Default | |", "|---|---|---|---|");
  for (const a of s.args) {
    const arg = a.positional ? `\`${a.name}${a.rest ? "…" : ""}\` (positional)` : `\`:${a.name}\``;
    const d = a.default !== undefined ? `\`${a.default}\`` : a.optional ? "none" : a.rest ? "one or more" : "required";
    out.push(`| ${arg} | ${typeText(a.type, a)} | ${d} | ${a.doc} |`);
  }
  if (s.examples.length > 0) out.push("", "Examples:");
  for (const e of s.examples) out.push(`- \`${e.text}\`: ${e.says}${e.fixture ? ` (${e.fixture})` : ""}`);
  if (s.hints.length > 0) out.push("", "Hints:");
  for (const h of s.hints) out.push(`- \`${h.code}\`${h.value ? ` (${h.value})` : ""}: ${h.text}`);
  return out.join("\n");
}

/** The reference: every symbol's full entry in markdown, in the skill's sections, then the issue codes. */
export function referenceText(language: Language): string {
  const symbols = [...language.symbols.values()];
  const sections = SECTIONS.map(([title, kinds]) => {
    const own = symbols.filter((s) => kinds.includes(s.kind));
    return [`## ${title[0]}${title.slice(1).toLowerCase()}`, ...own.map(entry)].join("\n\n");
  });
  const codes = ["## Issue codes", "", "| Code | Severity |", "|---|---|"];
  for (const c of SeqqlIssueCodes) codes.push(`| \`${c.code}\` | ${c.severity} |`);
  return `# SeqQL reference\n\n${sections.join("\n\n")}\n\n${codes.join("\n")}\n`;
}
