/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// QL1: the base table: every v1 symbol, both names, complete metadata.
import { describe, expect, it } from "vitest";
import { BASE_SYMBOLS, createLanguage, SeqqlIssueCodes, type SymbolDef, type Type } from "../../src/query/index.js";

// Every code seqspec:seqql#errors names: a hint's code must be a registered one.
const CODES = new Set(SeqqlIssueCodes.map((c) => c.code));

/** The heads of an S-expression's lists, strings skipped: enough to check an example names real symbols. */
const heads = (text: string): string[] =>
  [...text.replace(/"(?:[^"\\]|\\.)*"/g, '""').matchAll(/\(([^\s()"]+)/g)].map((m) => m[1] ?? "");
const keywords = (signature: string) => [...signature.matchAll(/:([a-z][a-z0-9-]*)/g)].map((m) => m[1]).sort();

function validDefault(type: Type, v: unknown): boolean {
  if (typeof type === "object") return type.enum.includes(v as string);
  if (type === "int") return Number.isInteger(v) && (v as number) >= 0;
  if (type === "num") return typeof v === "number";
  if (type === "bool") return typeof v === "boolean";
  return typeof v === "string";
}

const { language } = createLanguage();
const byName = (name: string) => language.symbols.get(name) as SymbolDef;

describe("the base table", () => {
  it("has 62 symbols with 53 short names, each kind counted as the spec counts it", () => {
    expect(BASE_SYMBOLS).toHaveLength(62);
    expect(new Set(BASE_SYMBOLS.map((s) => s.name)).size).toBe(62);
    expect(new Set(BASE_SYMBOLS.map((s) => s.alias)).size).toBe(53);
    const count = (pred: (s: SymbolDef) => boolean) => BASE_SYMBOLS.filter(pred).length;
    expect(count((s) => s.name.startsWith("core."))).toBe(11);
    expect(
      count((s) => ["seq.combinator.intersect", "seq.combinator.merge", "seq.modifier.complement"].includes(s.name)),
    ).toBe(3);
    expect(count((s) => s.kind === "generator" || s.kind === "macro")).toBe(12);
    expect(count((s) => s.kind === "modifier" && s.name !== "seq.modifier.complement")).toBe(7);
    expect(count((s) => s.kind === "filter")).toBe(3);
    expect(count((s) => s.kind === "property")).toBe(26);
  });

  it("names each symbol in its namespace, its kind matching the full name", () => {
    for (const s of BASE_SYMBOLS) {
      const ns = s.name.split(".").slice(0, 2).join(".");
      const expected: Record<string, string> = {
        "seq.generator": "generator",
        "seq.modifier": "modifier",
        "seq.combinator": "combinator",
        "seq.filter": "filter",
        "seq.macro": "macro",
      };
      if (s.name.startsWith("core.")) expect(s.kind, s.name).toBe("core");
      else if (expected[ns]) expect(s.kind, s.name).toBe(expected[ns]);
      else {
        expect(s.kind, s.name).toBe("property");
        expect(s.reads, s.name).toBe(s.name.split(".")[1]);
      }
      if (s.kind !== "property") expect(s.reads, s.name).toBeUndefined();
    }
  });

  it("gives every symbol complete metadata", () => {
    for (const s of BASE_SYMBOLS) {
      const at = s.name;
      expect(s.alias, at).toBeTruthy();
      expect(s.signature.startsWith(`(${s.alias}`) && s.signature.endsWith(")"), at).toBe(true);
      expect(s.doc.length > 0 && !s.doc.includes("\n"), at).toBe(true);
      expect(s.examples.length >= 1 && s.examples.length <= 2, at).toBe(true);
      for (const e of s.examples) {
        expect(e.says.length, at).toBeGreaterThan(0);
        // the example uses the symbol, and every head in it is a short or a full name of the table
        expect(heads(e.text), `${at}: ${e.text}`).toContain(s.alias);
        for (const h of heads(e.text))
          expect(language.aliases.has(h) || language.symbols.has(h), `${at}: ${h}`).toBe(true);
        if (e.fixture !== undefined) expect(e.fixture, at).toMatch(/^case\d+$/);
      }
      for (const h of s.hints) {
        expect(CODES.has(h.code), `${at}: ${h.code}`).toBe(true);
        expect(h.text.length, at).toBeGreaterThan(0);
      }
      expect(s.result, at).toBeTruthy();
      // the signature names exactly the keywords, and positionals come first, a rest one last
      const kws = s.args.filter((a) => !a.positional);
      expect(keywords(s.signature), at).toEqual(kws.map((a) => a.name).sort());
      const firstKw = s.args.findIndex((a) => !a.positional);
      if (firstKw >= 0)
        expect(
          s.args.slice(firstKw).every((a) => !a.positional),
          at,
        ).toBe(true);
      const positionals = s.args.filter((a) => a.positional);
      positionals.forEach((a, i) => {
        if (a.rest) expect(i, at).toBe(positionals.length - 1);
        if (a.optional)
          expect(
            positionals.slice(i).every((b) => b.optional || b.rest),
            at,
          ).toBe(true);
      });
      for (const a of s.args) {
        expect(a.doc.length, `${at} ${a.name}`).toBeGreaterThan(0);
        if (a.default !== undefined) expect(validDefault(a.type, a.default), `${at} ${a.name}`).toBe(true);
        expect(a.binds !== undefined, `${at} ${a.name}`).toBe(a.name === "where");
        if (a.binds) expect(a.type, at).toBe("test");
      }
    }
  });

  it("binds the slots of the binder table", () => {
    const binds = BASE_SYMBOLS.flatMap((s) => s.args.filter((a) => a.binds).map((a) => `${s.alias} ${a.binds}`));
    expect(binds.sort()).toEqual(
      ["systems system", "rows system", "items item", "partners item", "positions position", "pick group"].sort(),
    );
  });

  it("gives selections the query type, tests bool, and properties their results", () => {
    for (const s of BASE_SYMBOLS) {
      if (["generator", "macro", "modifier", "combinator", "filter"].includes(s.kind))
        expect(s.result, s.name).toBe("query");
      if (s.kind === "core") expect(s.result, s.name).toBe("bool");
    }
    expect(byName("seq.position.index").result).toBe("int");
    expect(byName("seq.position.value_of").result).toBe("num");
    expect(byName("seq.item.field").result).toBe("json");
    expect(byName("seq.group.overlaps").result).toBe("bool");
  });

  it("carries A.3's SELECTIONS and COMBINE lines, in its order", () => {
    const lines = BASE_SYMBOLS.filter((s) => s.kind !== "core" && s.kind !== "property").map(
      (s) => `${s.signature}  ${s.doc}`,
    );
    expect(lines).toEqual([...A3_COMBINE.slice(0, 3), ...A3_SELECTIONS, ...A3_COMBINE.slice(3)]);
  });

  it("marks A.3's EXAMPLES for the skill, and no others", () => {
    const skill = BASE_SYMBOLS.flatMap((s) => s.examples.filter((e) => e.skill).map((e) => `${e.text}  ${e.says}`));
    expect(skill.sort()).toEqual([...A3_EXAMPLES].sort());
  });

  it("orders properties by slot as A.3's TESTS lists them", () => {
    const bySlot = (slot: string) => BASE_SYMBOLS.filter((s) => s.reads === slot).map((s) => s.alias);
    expect(bySlot("position")).toEqual(["index", "letter", "label", "value-of"]);
    expect(bySlot("system")).toEqual(["ref", "kind", "label", "length", "origin", "score", "alphabet", "field"]);
    expect(bySlot("item")).toEqual(["id", "label", "category", "value", "type", "annotation", "length", "field"]);
    expect(bySlot("group")).toEqual(["length", "runs", "score", "overlaps", "within", "count"]);
  });
});

// seqspec:seqql#overview Appendix A.3, without Stories' `chain` macro.
const A3_SELECTIONS = [
  "(all)  every unit and alignment row, one group each",
  "(empty)  nothing",
  "(sel S)  the shorthand",
  "(systems :of R :kind K :where T)  whole sequences of kind unit, molecule, row or columns, one group each",
  "(rows :of A :where T)  an alignment's rows, one group each",
  "(items :of A :where T)  items of annotation A (all annotations if no :of), one group each",
  "(values :of A :min x :max x :above x :below x)  positions where values annotation A passes",
  "(letters L :in Q)  positions with one of the letters L",
  "(positions :in Q :where T)  positions passing a test",
  "(search P :in Q :edits n :classes C :regex b :ends e)  hits of a pattern, one group each; plain letters (ATG) match themselves, :classes only reads codes such as X or N",
  "(input name)  a selection the app passes",
  "(current)  inside each, the group being run on",
];
const A3_COMBINE = [
  "(and Q Q...)  each group of the first, intersected with the others",
  "(or Q Q...)  all groups, duplicates dropped",
  "(not Q :in U)  the positions of U (default all) not in Q",
  "(range Q part... :numbering N)  narrow to labels, inclusive",
  "(slice Q start stop...)  narrow to indices, from 0, end excluded",
  "(expand Q n :left n :right n)  widen each run",
  "(onto Q R :kind K :keep b)  the same places in another sequence, codons too; onto an alignment: its columns, :kind row its rows",
  "(group Q :by G)  regroup: none (one group), run or system",
  "(partners Q :via A :side S :where T)  the ends of pairs (annotation A) whose other end is in Q",
  "(each Q E)  run E per group of Q",
  "(pick Q :where T)  keep groups passing T",
  "(first Q n)  the first n groups",
  "(ids Q id...)  keep items with these ids",
];
const A3_EXAMPLES = [
  '(sel "1brs/A[10-20]")  residues 10-20 of 1brs/A',
  '(sel "p53@0:50")  the first 50 positions, as indices',
  "(pick (items :of domains) :where (matches (label) kinase))  kinase domains",
  "(values :of cons :above 0.8)  columns with conservation above 0.8",
  "(onto (input selection) msa)  the alignment columns of the selection",
  "(and 1brs/A (partners 1brs/D :via contacts))  residues of A in contact with D",
  "(pick sites :where (overlaps (pick domains :where (matches (label) SH2))))  sites in an SH2 domain (a bare name selects all of the annotation)",
  "(and 1brs/A (letters C))  cysteines of 1brs/A",
  "(rows :of msa :where (>= (field identity :of traits) 0.9))  rows at least 90% identical",
  '(and (expand (sel "1brs/A[10-20]") 5) (not (sel "1brs/A[10-20]")))  around a site, not the site',
  "(search GXGXXG :in 1brs/A :classes protein)  a motif",
];
