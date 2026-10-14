/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The v1 base table: 62 symbols, 53 short names, with the metadata the checker, the skill and the reference read. The
// order is the skill's: generators and macros as A.3's SELECTIONS lists them, combinators, modifiers and filters as its
// COMBINE, properties by slot as its TESTS.
import { expandSel, normalizeSel, SEL_STRING } from "./sel.ts";
import { defineSymbol } from "./symbol.ts";
import type { ArgDef, Slot, SymbolExample, SymbolHint, SymbolKind, SymbolSpec, Type } from "./types.ts";

type Meta = Pick<SymbolSpec, "doc"> &
  Partial<Pick<SymbolSpec, "args" | "examples" | "hints" | "result" | "expand" | "normalize">>;

const SYSTEM_KINDS: Type = { enum: ["unit", "molecule", "row", "columns"] };
const ENDS: Type = { enum: ["local_best_last", "local_best_first", "all"] };
/** Its place in the skill's EXAMPLES (A.3: simple first, then composed). */
const skill = (n: number) => ({ skill: n });

type More = Partial<ArgDef>;
const arg = (name: string, type: Type, doc: string, more?: More): ArgDef => ({ name, type, doc, ...more });
const pos = (name: string, type: Type, doc: string, more?: More) => arg(name, type, doc, { positional: true, ...more });
const rest = (name: string, type: Type, doc: string): ArgDef => pos(name, type, doc, { rest: true });
const opt = (name: string, type: Type, doc: string, more?: More) => arg(name, type, doc, { optional: true, ...more });
const where = (slot: Slot, doc: string): ArgDef => opt("where", "test", doc, { binds: slot });
const q = (doc = "the selection"): ArgDef => pos("Q", "query", doc);
const ex = (text: string, says: string, more?: Partial<SymbolExample>): SymbolExample => ({ text, says, ...more });
const hint = (code: string, text: string, value?: string): SymbolHint => ({ code, text, ...(value ? { value } : {}) });

function sym<const N extends string>(name: N, alias: string, kind: SymbolKind, signature: string, m: Meta) {
  return defineSymbol({ name, alias, kind, signature, args: [], result: "query", examples: [], hints: [], ...m });
}

const BINDERS: Record<Slot, string> = {
  position: "(positions :where …)",
  system: "(systems :where …) or (rows :where …)",
  item: "(items :where …) or (partners :where …)",
  group: "(pick Q :where …)",
};
const A_SLOT: Record<Slot, string> = { position: "a position", system: "a system", item: "an item", group: "a group" };
const orList = (xs: readonly string[]) =>
  xs.length < 2 ? (xs[0] ?? "") : `${xs.slice(0, -1).join(", ")} or ${xs.at(-1)}`;

/** `seqql.slot`'s message: what a property reads and where it can be used, or "inside a :where" when several. */
export function slotMessage(name: string, slots: readonly Slot[], anyWhere = false): string {
  const use = anyWhere ? "use it inside a :where" : `use it in ${orList(slots.map((s) => BINDERS[s]))}`;
  return `${name} reads ${orList(slots.map((s) => A_SLOT[s]))}: ${use}`;
}

/**
 * A property `seq.<slot>.<short name>`: reads its slot (seqspec:seqql#elements-and-groups); its first hint says where
 * it can be used.
 */
function prop<const N extends `seq.${Slot}.${string}`>(
  name: N,
  result: Type,
  m: Meta & { signature?: string; alias?: string },
) {
  const [, reads, last] = name.split(".") as [string, Slot, string];
  const alias = m.alias ?? last;
  const where = hint("seqql.slot", slotMessage(alias, [reads]));
  return defineSymbol({
    name,
    alias,
    kind: "property",
    reads,
    signature: `(${alias})`,
    args: [],
    ...m,
    result,
    hints: [where, ...(m.hints ?? [])],
    examples: m.examples ?? [],
  });
}

/** A relation `(op a b)` over JSON values or numbers. */
const relation = <const N extends string>(name: N, alias: string, type: Type, m: Meta) =>
  sym(name, alias, "core", `(${alias} a b)`, {
    result: "bool",
    args: [pos("a", type, "the left side"), pos("b", type, "the right side")],
    hints: type === "num" ? [hint("seqql.type", `${alias} compares numbers: (${alias} (length) 50)`)] : [],
    ...m,
  });

export const BASE_SYMBOLS = [
  // Core: logic over tests, relations, matching.
  sym("core.logic.and", "and", "core", "(and T T...)", {
    doc: "every test holds; stops at the first that fails",
    result: "bool",
    args: [pos("T", "bool", "a test"), rest("T", "bool", "more tests")],
    examples: [
      ex("(items :of regions :where (and (> (length) 30) (matches (label) binding)))", "long binding regions", {
        fixture: "case1",
      }),
    ],
    hints: [hint("seqql.type", "and mixes a selection and a test: filter with (pick Q :where T)")],
  }),
  sym("core.logic.or", "or", "core", "(or T T...)", {
    doc: "some test holds; stops at the first that does",
    result: "bool",
    args: [pos("T", "bool", "a test"), rest("T", "bool", "more tests")],
    examples: [
      ex("(positions :in p53 :where (or (= (letter) C) (= (letter) H)))", "cysteines and histidines of p53", {
        fixture: "case1",
      }),
    ],
    hints: [hint("seqql.type", "or mixes a selection and a test: filter with (pick Q :where T)")],
  }),
  sym("core.logic.not", "not", "core", "(not T)", {
    doc: "the test fails; true where its value is missing",
    result: "bool",
    args: [pos("T", "bool", "a test")],
    examples: [
      ex("(items :of regions :where (not (matches (label) binding)))", "regions other than binding ones", {
        fixture: "case1",
      }),
    ],
  }),
  relation("core.rel.eq", "=", "json", {
    doc: "equal as JSON values; a missing value equals nothing",
    examples: [ex("(items :where (= (type) uniprot:Domain))", "every domain item")],
  }),
  relation("core.rel.neq", "!=", "json", {
    doc: "not equal as JSON values",
    examples: [ex("(items :of variants :where (!= (category) benign))", "variants not benign")],
  }),
  relation("core.rel.lt", "<", "num", {
    doc: "less than; false when a side is missing or not a number",
    examples: [ex("(positions :in p53 :where (< (index) 10))", "the first ten positions of p53", { fixture: "case1" })],
  }),
  relation("core.rel.lte", "<=", "num", {
    doc: "at most; false when a side is missing or not a number",
    examples: [
      ex("(items :of regions :where (<= (length) 40))", "regions of at most 40 positions", { fixture: "case1" }),
    ],
  }),
  relation("core.rel.gr", ">", "num", {
    doc: "greater than; false when a side is missing or not a number",
    examples: [ex("(systems :where (> (length) 100))", "sequences longer than 100", { fixture: "case1" })],
  }),
  relation("core.rel.gre", ">=", "num", {
    doc: "at least; false when a side is missing or not a number",
    examples: [ex("(pick (letters C) :where (>= (length) 2))", "runs of two or more cysteines", { fixture: "case14" })],
  }),
  sym("core.rel.in", "in", "core", "(in x v...)", {
    doc: "x equals one of the v",
    result: "bool",
    args: [pos("x", "json", "the value"), rest("v", "json", "the values it may equal")],
    examples: [
      ex("(positions :in p53 :where (in (letter) C H))", "cysteines and histidines of p53", { fixture: "case1" }),
    ],
    hints: [hint("seqql.arg.missing", "in needs a value and what it may equal: (in (letter) C H)")],
  }),
  sym("core.str.matches", "matches", "core", "(matches s pattern :case b)", {
    doc: "a regular expression, case-insensitive unless :case true",
    result: "bool",
    args: [
      pos("s", "str", "the text"),
      pos("pattern", "str", "an ECMAScript regular expression"),
      opt("case", "bool", "case-sensitive", { default: false }),
    ],
    examples: [ex("(rows :of fam :where (matches (label) HUMAN :case true))", "the human members' rows")],
    hints: [hint("seqql.arg.missing", "matches needs the text and a pattern: (matches (label) kinase)")],
  }),

  // Combining selections: the same short names as the logic, told apart by operand kind
  // (seqspec:seqql#binding-and-checking).
  sym("seq.combinator.intersect", "and", "combinator", "(and Q Q...)", {
    doc: "each group of the first, intersected with the others",
    args: [q("the groups kept"), rest("Q", "query", "the selections they are intersected with")],
    examples: [ex("(and 1brs/A domains)", "positions of 1brs/A inside any domain")],
    hints: [hint("seqql.type", "and mixes a selection and a test: filter with (pick Q :where T)")],
  }),
  sym("seq.combinator.merge", "or", "combinator", "(or Q Q...)", {
    doc: "all groups, duplicates dropped",
    args: [q(), rest("Q", "query", "more selections")],
    examples: [
      ex('(or (sel "1brs/A[10-20]") (sel "1brs/B[30-35]"))', "two sites, one group each", { fixture: "case7" }),
    ],
  }),
  sym("seq.modifier.complement", "not", "modifier", "(not Q :in U)", {
    doc: "the positions of U (default all) not in Q",
    args: [q("the positions removed"), opt("in", "query", "the positions kept from; default (all)")],
    examples: [ex("(not domains :in 1brs/A)", "positions of 1brs/A outside every domain")],
    hints: [hint("seqql.arg.count", "not takes one selection; to remove B from A, intersect: (and A (not B))")],
  }),

  // Generators and macros (SELECTIONS).
  sym("seq.generator.all", "all", "generator", "(all)", {
    doc: "every unit and alignment row, one group each",
    examples: [ex("(and (all) (letters W))", "every tryptophan", { fixture: "case1" })],
  }),
  sym("seq.generator.empty", "empty", "generator", "(empty)", {
    doc: "nothing",
    examples: [ex("(empty)", "no groups")],
  }),
  sym("seq.macro.sel", "sel", "macro", "(sel S)", {
    doc: "the shorthand",
    expand: expandSel,
    normalize: normalizeSel,
    args: [
      pos("S", "str", "refs with label ranges [10-20], index slices @9:20, item ids [#id], inputs $name; a, b a union"),
    ],
    examples: [
      ex('(sel "1brs/A[10-20]")', "residues 10-20 of 1brs/A", { skill: 1, fixture: "case7" }),
      ex('(sel "p53@0:50")', "the first 50 positions, as indices", { skill: 2, fixture: "case1" }),
    ],
    hints: [
      hint("seqql.arg.missing", SEL_STRING),
      hint("seqql.type", SEL_STRING),
      hint("seqql.sel.syntax", 'labels in brackets, indices after @: (sel "1brs/A[10-20]"), (sel "1brs/A@9:20")'),
      hint("seqql.sel.mixed", "a bracket holds labels or #ids, not both"),
    ],
  }),
  sym("seq.generator.systems", "systems", "generator", "(systems :of R :kind K :where T)", {
    doc: "whole sequences of kind unit, molecule, row or columns, one group each",
    args: [
      opt("of", "ref", "a molecule or complex for units, an alignment for rows; default every system"),
      opt("kind", SYSTEM_KINDS, "the kind of system", { default: "unit" }),
      where("system", "the systems kept"),
    ],
    examples: [ex("(systems :of 1brs)", "the units of complex 1brs, one group each", { fixture: "case7" })],
    hints: [hint("seqql.arg.unknown", "it takes :of, :kind and :where; narrow positions with (and Q ...)")],
  }),
  sym("seq.generator.rows", "rows", "generator", "(rows :of A :where T)", {
    doc: "an alignment's rows, one group each",
    args: [opt("of", "ref", "the alignment; default every alignment"), where("system", "the rows kept")],
    examples: [
      ex("(rows :of msa :where (>= (field identity :of traits) 0.9))", "rows at least 90% identical", skill(9)),
    ],
    hints: [hint("seqql.ref.kind", ":of takes an alignment (msa); a row alone is msa/3")],
  }),
  sym("seq.generator.items", "items", "generator", "(items :of A :where T)", {
    doc: "items of annotation A (all annotations if no :of), one group each",
    args: [opt("of", "ref", "the annotation; default every annotation with items"), where("item", "the items kept")],
    examples: [ex("(items :of regions)", "p53's regions, one group each", { fixture: "case1" })],
    hints: [
      hint(
        "seqql.arg.unknown",
        "it takes :of and :where; items touching a selection: (pick (items :of A) :where (overlaps Q))",
      ),
      hint("seqql.ref.kind", "values and matrices have no items: select by value with (values :of A :above x)"),
    ],
  }),
  sym("seq.generator.values", "values", "generator", "(values :of A :min x :max x :above x :below x)", {
    doc: "positions where values annotation A passes",
    args: [
      arg("of", "ref", "a values annotation"),
      opt("min", "num", "at least"),
      opt("max", "num", "at most"),
      opt("above", "num", "greater than"),
      opt("below", "num", "less than"),
    ],
    examples: [
      ex("(values :of cons :above 0.8)", "columns with conservation above 0.8", { skill: 4, fixture: "case3" }),
    ],
    hints: [
      hint(
        "seqql.ref.kind",
        ":of takes a values annotation (cons, plddt); for items by value: (pick (items :of A) :where (> (value) 0.8))",
      ),
    ],
  }),
  sym("seq.generator.letters", "letters", "generator", "(letters L :in Q)", {
    doc: "positions with one of the letters L",
    args: [
      pos("L", "str", "the letters, case-insensitive: C, CHW"),
      opt("in", "query", "where to look; default (all)"),
    ],
    examples: [
      ex("(and 1brs/A (letters C))", "cysteines of 1brs/A", skill(8)),
      ex("(letters CHW :in p53)", "cysteines, histidines and tryptophans of p53", { fixture: "case1" }),
    ],
    hints: [
      hint("seqql.arg.missing", "letters needs the letters: (letters C)"),
      hint("seqql.limit.positions", "scan less: (and 1brs/A (letters C)) reads 1brs/A only"),
    ],
  }),
  sym("seq.generator.positions", "positions", "generator", "(positions :in Q :where T)", {
    doc: "positions passing a test",
    args: [opt("in", "query", "where to look; default (all)"), where("position", "the test each position passes")],
    examples: [ex("(positions :in p53 :where (< (index) 10))", "the first ten positions of p53", { fixture: "case1" })],
    hints: [
      hint("seqql.limit.positions", "a test runs once per position: give it :in Q, or (and Q (positions :where T))"),
    ],
  }),
  sym("seq.generator.search", "search", "generator", "(search P :in Q :edits n :classes C :regex b :ends e)", {
    doc: "hits of a pattern, one group each; plain letters (ATG) match themselves, :classes only reads codes such as X or N",
    args: [
      pos("P", "str", "the pattern"),
      opt("in", "query", "where to look; default every unit"),
      opt("edits", "int", "edits allowed", { default: 0 }),
      opt("classes", "str", "a letter-class table: the pattern's letters are classes (protein, iupac)"),
      opt("regex", "bool", "the pattern is a regular expression", { default: false }),
      opt("ends", ENDS, "which ends an approximate search reports", { default: "local_best_last" }),
    ],
    examples: [
      ex("(search GXGXXG :in 1brs/A :classes protein)", "a motif", skill(11)),
      ex(
        '(search "[KR]{2}.{10,12}[KR]{3,}" :in p53 :regex true)',
        "a regular expression, quoted: the bipartite nuclear localization signal",
        { fixture: "case1" },
      ),
    ],
    hints: [
      hint("seqql.arg.missing", "search needs a pattern: (search GXGXXG)"),
      hint("seqql.limit.positions", "search less: (search GXGXXG :in 1brs/A) reads 1brs/A only"),
    ],
  }),
  sym("seq.generator.input", "input", "generator", "(input name)", {
    doc: "a selection the app passes",
    args: [pos("name", "str", "the input's name: selection, focus, structure")],
    examples: [ex("(input selection)", "the person's selection")],
    hints: [hint("seqql.input.unknown", "inputs are named by the app: (input selection)")],
  }),
  sym("seq.generator.current", "current", "generator", "(current)", {
    doc: "inside each, the group being run on",
    examples: [
      ex("(each (items :of regions) (expand (current) 5))", "every region widened by 5", { fixture: "case1" }),
    ],
    hints: [hint("seqql.current", "(current) is only meaningful inside each")],
  }),

  // Modifiers and filters (COMBINE).
  sym("seq.modifier.range", "range", "modifier", "(range Q part... :numbering N)", {
    doc: "narrow to labels, inclusive",
    args: [
      q(),
      rest("part", "label", "a label or label range: 10-20, 45A, -5-10"),
      opt("numbering", "str", "the numbering; default the system's"),
    ],
    examples: [
      ex("(range 1brs/A 10-20 45)", "labels 10 to 20 and 45 of 1brs/A", { fixture: "case7" }),
      ex("(range p53 100-120)", "residues 100 to 120 of p53", { fixture: "case1" }),
    ],
    hints: [
      hint("seqql.arg.order", "positional arguments come before keywords: (range 1brs/A 10-20 :numbering uniprot)"),
      hint("seqql.label.unknown", "labels are as the ruler shows them; for indices from 0 use (slice Q start stop)"),
    ],
  }),
  sym("seq.modifier.slice", "slice", "modifier", "(slice Q start stop...)", {
    doc: "narrow to indices, from 0, end excluded",
    args: [q(), rest("bounds", "int", "start and stop indices, in pairs: 9 20, or 9 20 30 40")],
    examples: [ex("(slice p53 0 50)", "the first 50 positions of p53, as indices", { fixture: "case1" })],
    hints: [
      hint("seqql.arg.count", "slice takes start-stop pairs: (slice 1brs/A 9 20)"),
      hint("seqql.type", "slice takes indices, whole numbers from 0: for labels like 45A use (range 1brs/A 45A)"),
      hint("seqql.slice.range", "indices run from 0 to the length; the summary gives each run's slice"),
      hint("seqql.slice.clipped", "indices run from 0 to the length; the summary gives each run's slice"),
    ],
  }),
  sym("seq.modifier.expand", "expand", "modifier", "(expand Q n :left n :right n)", {
    doc: "widen each run",
    args: [
      q(),
      pos("n", "int", "positions added on both sides", { optional: true }),
      opt("left", "int", "positions added before"),
      opt("right", "int", "positions added after"),
    ],
    examples: [
      ex('(and (expand (sel "1brs/A[10-20]") 5) (not (sel "1brs/A[10-20]")))', "around a site, not the site", {
        skill: 10,
        fixture: "case7",
      }),
    ],
    hints: [hint("seqql.arg.count", "expand takes one positional after the selection, the width: (expand Q 5)")],
  }),
  sym("seq.modifier.onto", "onto", "modifier", "(onto Q R :kind K :keep b)", {
    doc: "the same places in another sequence, codons too; onto an alignment: its columns, :kind row its rows",
    args: [
      q("the positions carried"),
      pos("R", "ref", "a unit, molecule, row or alignment", { optional: true }),
      opt(
        "kind",
        SYSTEM_KINDS,
        "only for an alignment's rows (row) or columns (columns, the default onto an alignment); without R, every reachable system of the kind; none across a codon mapping",
      ),
      opt("keep", "bool", "keep the source runs too", { default: false }),
    ],
    examples: [
      ex("(onto (input selection) msa)", "the alignment columns of the selection", skill(5)),
      ex('(onto (sel "protein[5-10]") gene)', "the codons of residues 5-10, no :kind needed", { fixture: "case14" }),
    ],
    hints: [
      hint(
        "seqql.onto.unmapped",
        "a group with no image is dropped; a row reaches other sequences through its binding",
      ),
      hint(
        "seqql.type",
        ":kind is a kind of sequence (unit, molecule, row, columns), and codon is not one: across a codon mapping leave :kind out, (onto Q gene) gives the codons",
        "codon",
      ),
    ],
  }),
  sym("seq.modifier.group", "group", "modifier", "(group Q :by G)", {
    doc: "regroup: none (one group), run or system",
    args: [
      q(),
      arg("by", { enum: ["none", "run", "system"] }, "none: one group; run: each maximal run; system: one per system"),
    ],
    examples: [ex("(group (letters C :in p53) :by none)", "p53's cysteines as one group", { fixture: "case1" })],
    hints: [hint("seqql.arg.missing", "group needs :by: (group Q :by none)")],
  }),
  sym("seq.modifier.partners", "partners", "modifier", "(partners Q :via A :side S :where T)", {
    doc: "the ends of pairs (annotation A) whose other end is in Q",
    args: [
      q("the partner selection"),
      arg("via", "ref", "a pairs annotation"),
      opt(
        "side",
        { enum: ["any", "a", "b"] },
        "which ends may be returned: a or b, the pairs' columns (any: both); not a chain, which goes in Q",
        { default: "any" },
      ),
      where("item", "the pairs to use"),
    ],
    examples: [
      ex("(and 1brs/A (partners 1brs/D :via contacts))", "residues of A in contact with D", {
        skill: 6,
        fixture: "case7",
      }),
    ],
    hints: [hint("seqql.ref.kind", ":via takes a pairs annotation (contacts, hbonds)")],
  }),
  sym("seq.modifier.each", "each", "modifier", "(each Q E)", {
    doc: "run E per group of Q",
    args: [q("the groups run on"), pos("E", "query", "the query run per group, (current) the group")],
    examples: [
      ex("(each (items :of regions) (and (current) (letters C)))", "the cysteines of each region, one group each", {
        fixture: "case1",
      }),
    ],
  }),
  sym("seq.filter.pick", "pick", "filter", "(pick Q :where T)", {
    doc: "keep groups passing T",
    args: [q(), { ...arg("where", "test", "the test each group passes"), binds: "group" }],
    examples: [ex("(pick (items :of domains) :where (matches (label) kinase))", "kinase domains", skill(3))],
    hints: [hint("seqql.arg.missing", "pick needs a test: (pick Q :where (> (length) 50))")],
  }),
  sym("seq.filter.first", "first", "filter", "(first Q n)", {
    doc: "the first n groups",
    args: [q(), pos("n", "int", "how many", { optional: true, default: 1 })],
    examples: [ex("(first (items :of regions) 2)", "the first two regions", { fixture: "case1" })],
  }),
  sym("seq.filter.ids", "ids", "filter", "(ids Q id...)", {
    doc: "keep items with these ids",
    args: [q("the item groups"), rest("id", "str", "item ids")],
    examples: [ex("(ids (items :of domains) PF00069 PF07714)", "two domains by id")],
    hints: [hint("seqql.arg.missing", "ids needs at least one id: (ids (items :of domains) PF00069)")],
  }),

  // Properties, by slot (TESTS).
  prop("seq.position.index", "int", {
    doc: "its index (0-based, as slice takes it)",
    examples: [
      ex("(positions :in 1brs/A :where (< (index) 10))", "the first ten positions of 1brs/A", { fixture: "case7" }),
    ],
  }),
  prop("seq.position.letter", "str", {
    doc: 'the letter ("" where letters are not held)',
    examples: [ex("(positions :in p53 :where (= (letter) W))", "tryptophans of p53", { fixture: "case1" })],
  }),
  prop("seq.position.label", "str", {
    doc: 'its label in the default numbering or :numbering N ("" when none)',
    signature: "(label :numbering N)",
    args: [opt("numbering", "str", "the numbering; default the system's")],
    examples: [ex('(positions :in 1brs/A :where (matches (label) "[A-Z]$"))', "positions with an insertion code")],
  }),
  // A full name is snake_case, as every name a document holds (CM-15); the text keeps its alias, `value-of`.
  prop("seq.position.value_of", "num", {
    alias: "value-of",
    doc: "the value of values annotation A here; NaN when none",
    signature: "(value-of A)",
    args: [pos("A", "ref", "a values annotation")],
    examples: [ex("(positions :in model :where (> (value-of plddt) 70))", "confident positions", { fixture: "case8" })],
    hints: [hint("seqql.ref.kind", "value-of takes a values annotation (cons, plddt)")],
  }),
  prop("seq.system.ref", "str", {
    doc: "the system's ref",
    examples: [ex("(systems :where (matches (ref) ^1brs/))", "the units of 1brs", { fixture: "case7" })],
  }),
  prop("seq.system.kind", "str", {
    doc: "the system's kind: unit, molecule, row or columns",
    examples: [ex("(pick (all) :where (= (kind) row))", "every alignment row", { fixture: "case3" })],
  }),
  prop("seq.system.label", "str", {
    doc: "the system's label",
    examples: [ex("(rows :of fam :where (matches (label) HUMAN))", "the human members' rows")],
  }),
  prop("seq.system.length", "int", {
    doc: "positions in the system",
    examples: [ex("(systems :where (> (length) 100))", "sequences longer than 100", { fixture: "case1" })],
  }),
  prop("seq.system.origin", "str", {
    doc: "the system's origin",
    examples: [
      ex("(systems :where (= (origin) sequence))", "units that are their molecule's own sequence", {
        fixture: "case1",
      }),
    ],
  }),
  prop("seq.system.score", "num", {
    doc: "a unit's score",
    examples: [ex("(systems :where (>= (score) 0.5))", "units scoring at least 0.5")],
  }),
  prop("seq.system.alphabet", "str", {
    doc: "the system's alphabet: protein, dna, rna or generic",
    examples: [ex("(systems :where (= (alphabet) dna))", "DNA sequences", { fixture: "case14" })],
  }),
  prop("seq.system.field", "json", {
    doc: "another field: :of a row table, :ns app fields",
    signature: "(field name :of K :ns N)",
    args: [
      pos("name", "str", "the field: molecule, complex, alignment, row, or a column of :of or :ns"),
      opt("of", "ref", "a rows key table"),
      opt("ns", "str", "the namespace of an app's custom fields"),
    ],
    examples: [
      ex("(systems :where (= (field auth_asym_id :ns mvs) A))", "units of chain A, by the app's field", {
        fixture: "case7",
      }),
    ],
    hints: [
      hint("seqql.field.unknown", "a system's own fields: molecule, complex, alignment, row; others need :of or :ns"),
    ],
  }),
  prop("seq.item.id", "str", {
    doc: "the item's id",
    examples: [ex("(pick (items :of domains) :where (in (id) PF00069 PF07714))", "two domains by id")],
  }),
  prop("seq.item.label", "str", {
    doc: "the item's label",
    examples: [
      ex("(items :of regions :where (matches (label) binding))", "p53's DNA-binding region", { fixture: "case1" }),
    ],
  }),
  prop("seq.item.category", "str", {
    doc: "the item's category",
    examples: [ex("(partners 1brs/D :via contacts :where (= (category) hbond))", "ends hydrogen-bonded to 1brs/D")],
  }),
  prop("seq.item.value", "json", {
    doc: "the item's value",
    examples: [ex("(pick (items :of domains) :where (> (value) 0.8))", "domains valued above 0.8")],
  }),
  prop("seq.item.type", "str", {
    doc: "its annotation's semantic type",
    examples: [ex("(items :where (= (type) uniprot:Domain))", "every domain item")],
  }),
  prop("seq.item.annotation", "str", {
    doc: "its annotation's ref",
    examples: [ex("(items :where (= (annotation) variants))", "the items of variants")],
  }),
  prop("seq.item.length", "int", {
    doc: "positions in the item",
    examples: [ex("(items :of regions :where (> (length) 50))", "regions longer than 50", { fixture: "case1" })],
  }),
  prop("seq.item.field", "json", {
    doc: "any other column",
    signature: "(field name)",
    args: [pos("name", "str", "the column or field: tooltip, fragments, evalue")],
    examples: [ex("(items :of hits :where (< (field evalue) 0.001))", "hits with an e-value below 0.001")],
    hints: [hint("seqql.field.unknown", "an annotation's own columns are listed when the name is unknown")],
  }),
  prop("seq.group.length", "int", {
    doc: "positions in the group",
    examples: [ex("(pick (letters C) :where (>= (length) 2))", "runs of two or more cysteines", { fixture: "case14" })],
  }),
  prop("seq.group.runs", "int", {
    doc: "the group's run count",
    examples: [ex("(pick (items :of domains) :where (> (runs) 1))", "domains in more than one piece")],
  }),
  prop("seq.group.score", "num", {
    doc: "a hit's edit count; NaN otherwise",
    examples: [ex("(pick (search GXGXXG :edits 1) :where (= (score) 0))", "exact hits only")],
  }),
  prop("seq.group.overlaps", "bool", {
    doc: "shares a position with Q",
    signature: "(overlaps Q)",
    args: [q("computed once, not per group")],
    examples: [
      ex(
        "(pick sites :where (overlaps (pick domains :where (matches (label) SH2))))",
        "sites in an SH2 domain (a bare name selects all of the annotation)",
        skill(7),
      ),
      ex('(pick (items :of variants) :where (overlaps (sel "domains[#PF00069]")))', "variants in a domain, by its id"),
    ],
  }),
  prop("seq.group.within", "bool", {
    doc: "within n positions of Q, in the same system",
    signature: "(within Q n)",
    args: [q("computed once, not per group"), pos("n", "int", "the distance in positions")],
    examples: [ex('(pick (items :of sites) :where (within (sel "1brs/A[40-90]") 5))', "sites near a stretch")],
  }),
  prop("seq.group.count", "int", {
    doc: "how many of Q's positions it holds",
    signature: "(count Q)",
    args: [q("computed once, not per group")],
    examples: [
      ex("(pick (items :of domains) :where (>= (count (letters C)) 4))", "domains with four or more cysteines"),
    ],
  }),
];

export type BaseSymbolName = (typeof BASE_SYMBOLS)[number]["name"];
