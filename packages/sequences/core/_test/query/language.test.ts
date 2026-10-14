/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// QL1: the registry (seqspec:seqql#extension), the alias index (the binder's candidates), the builder
// (seqspec:seqql#overview) and the JSON form (seqspec:seqql#names).
import { describe, expect, it } from "vitest";
import {
  BASE_SYMBOLS,
  createLanguage,
  defineSymbol,
  type Expression,
  isExpression,
  type Language,
  Q,
  type Slot,
  type SymbolDef,
  type SymbolKind,
  type SymbolSpec,
} from "../../src/query/index.ts";

const ext = (name: string, alias: string | null, kind: SymbolKind, reads?: Slot, more?: Partial<SymbolSpec>) =>
  defineSymbol({
    name,
    alias,
    kind,
    ...(reads ? { reads } : {}),
    args: [],
    result: kind === "property" ? "str" : "query",
    signature: `(${alias ?? name})`,
    doc: "an extension's symbol",
    examples: [{ text: `(${alias ?? name})`, says: "it" }],
    hints: [],
    ...more,
  });

const names = (language: Language, alias: string) => language.aliases.get(alias)?.map((s) => s.name);
const SLOT_CHAINS: Record<string, readonly Slot[]> = {
  positions: ["position", "system"],
  systems: ["system"],
  rows: ["system"],
  items: ["item"],
  partners: ["item"],
  pick: ["group", "item", "system"],
};

/**
 * seqspec:seqql#binding-and-checking step 3 over the index's candidates, as far as QL1 can: the first slot of the
 * binder's chain with a candidate whose signature takes the keywords given. The binder (QL2a) does this on the tree.
 */
function resolveProperty(language: Language, binder: string, alias: string, given: readonly string[] = []) {
  const candidates = language.aliases.get(alias) ?? [];
  for (const slot of SLOT_CHAINS[binder] ?? []) {
    const fit = candidates.find(
      (s) => s.reads === slot && given.every((k) => s.args.some((a) => !a.positional && a.name === k)),
    );
    if (fit) return fit.name;
  }
  return null;
}

describe("the alias index", () => {
  const { language, issues } = createLanguage();

  it("holds the base table with no issues", () => {
    expect(issues).toEqual([]);
    expect(language.symbols.size).toBe(62);
    expect(language.aliases.size).toBe(53);
  });

  it("lists every candidate of a shared short name, in table order", () => {
    expect(names(language, "label")).toEqual(["seq.position.label", "seq.system.label", "seq.item.label"]);
    expect(names(language, "length")).toEqual(["seq.system.length", "seq.item.length", "seq.group.length"]);
    expect(names(language, "score")).toEqual(["seq.system.score", "seq.group.score"]);
    expect(names(language, "field")).toEqual(["seq.system.field", "seq.item.field"]);
    expect(names(language, "and")).toEqual(["core.logic.and", "seq.combinator.intersect"]);
    expect(names(language, "or")).toEqual(["core.logic.or", "seq.combinator.merge"]);
    expect(names(language, "not")).toEqual(["core.logic.not", "seq.modifier.complement"]);
    const shared = [...language.aliases].filter(([, list]) => list.length > 1).map(([a]) => a);
    expect(shared.sort()).toEqual(["and", "field", "label", "length", "not", "or", "score"]);
  });

  it("tells shared candidates apart by operand kind (and, or, not) or by slot (properties), as the binder needs", () => {
    for (const [alias, list] of language.aliases) {
      if (list.length < 2) continue;
      if (["and", "or", "not"].includes(alias)) {
        expect(list.map((s) => s.result)).toEqual(["bool", "query"]);
      } else {
        expect(list.every((s) => s.kind === "property")).toBe(true);
        expect(new Set(list.map((s) => s.reads)).size).toBe(list.length);
      }
    }
  });

  it("resolves properties along the slot chain from its candidates", () => {
    expect(resolveProperty(language, "pick", "label")).toBe("seq.item.label");
    expect(resolveProperty(language, "pick", "length")).toBe("seq.group.length");
    expect(resolveProperty(language, "pick", "kind")).toBe("seq.system.kind");
    expect(resolveProperty(language, "pick", "field")).toBe("seq.item.field");
    expect(resolveProperty(language, "pick", "field", ["of"])).toBe("seq.system.field");
    expect(resolveProperty(language, "positions", "label", ["numbering"])).toBe("seq.position.label");
    expect(resolveProperty(language, "positions", "length")).toBe("seq.system.length");
    expect(resolveProperty(language, "rows", "label")).toBe("seq.system.label");
    expect(resolveProperty(language, "partners", "category")).toBe("seq.item.category");
    expect(resolveProperty(language, "items", "letter")).toBeNull();
  });
});

describe("registration", () => {
  const molql = ext("molstar.select", "molql", "generator");
  const other = ext("other.molql", "molql", "modifier");

  it("adds an extension's symbols and short names, leaving the base as it was", () => {
    const chain = ext("stories.chain", "chain", "macro");
    const { language, issues } = createLanguage([{ symbols: [molql] }, { symbols: [chain] }]);
    expect(issues).toEqual([]);
    expect(language.symbols.get("molstar.select")).toBe(molql);
    expect(names(language, "molql")).toEqual(["molstar.select"]);
    expect(names(language, "chain")).toEqual(["stories.chain"]);
    // no module state: a language made after it has none of the extension's
    const plain = createLanguage().language;
    expect(plain.symbols.has("molstar.select")).toBe(false);
    expect(plain.aliases.has("molql")).toBe(false);
  });

  it("keeps a symbol with no short name reachable by its full name only", () => {
    const bare = ext("molstar.raw", null, "generator");
    const { language, issues } = createLanguage([{ symbols: [bare] }]);
    expect(issues).toEqual([]);
    expect(language.symbols.get("molstar.raw")).toBe(bare);
    expect([...language.aliases.values()].flat()).not.toContain(bare);
  });

  it("refuses a collision for both symbols, whatever the load order", () => {
    const orders = [[[molql], [other]], [[other], [molql]], [[other, molql]], [[molql, other]]];
    for (const order of orders) {
      const { language, issues } = createLanguage(order.map((symbols) => ({ symbols })));
      expect(language.aliases.has("molql")).toBe(false);
      expect(language.symbols.get("molstar.select")).toBe(molql);
      expect(language.symbols.get("other.molql")).toBe(other);
      expect(issues.map((i) => [i.code, i.message])).toEqual([
        [
          "seqql.alias.collision",
          'the short name "molql" is defined by molstar.select and other.molql, which the binder cannot tell apart; use the full names',
        ],
      ]);
    }
  });

  it("allows one short name for properties of different slots, and refuses it within one slot", () => {
    const a = ext("a.hydro", "hydro", "property", "item");
    const b = ext("b.hydro", "hydro", "property", "item");
    const c = ext("c.hydro", "hydro", "property", "system");
    const permutations = [
      [a, b, c],
      [c, b, a],
      [b, c, a],
    ];
    for (const symbols of permutations) {
      const { language, issues } = createLanguage(symbols.map((s) => ({ symbols: [s] })));
      expect(names(language, "hydro")).toEqual(["c.hydro"]);
      expect(issues.map((i) => i.code)).toEqual(["seqql.alias.collision"]);
      expect(issues[0]?.message).toContain("a.hydro and b.hydro");
    }
    const { language, issues } = createLanguage([{ symbols: [a, c] }]);
    expect(issues).toEqual([]);
    expect(names(language, "hydro")).toEqual(["a.hydro", "c.hydro"]);
    // a property and a selection with one short name: the binder has no slot to tell them by
    const gen = ext("d.hydro", "hydro", "generator");
    const both = createLanguage([{ symbols: [gen, c] }]);
    expect(both.language.aliases.has("hydro")).toBe(false);
    expect(both.issues.map((i) => i.code)).toEqual(["seqql.alias.collision"]);
  });

  it("refuses base short names and words the reader reads otherwise as aliases, keeping the symbols", () => {
    for (const alias of [
      "label",
      "and",
      "true",
      "false",
      "12",
      "-1.5e3",
      ":of",
      "$x",
      "a.b",
      "a b",
      "f(x)",
      "x,y",
      "",
    ]) {
      const s = ext("x.thing", alias, "property", "item");
      const { language, issues } = createLanguage([{ symbols: [s] }]);
      expect(
        issues.map((i) => i.code),
        alias,
      ).toEqual(["seqql.alias.reserved"]);
      expect(language.symbols.get("x.thing")).toBe(s);
      if (alias === "label") expect(names(language, "label")).toHaveLength(3);
    }
    for (const alias of ["molql", "chain", "hydro-2", "x:y", "n@1"]) {
      const { issues } = createLanguage([{ symbols: [ext("x.thing", alias, "generator")] }]);
      expect(issues, alias).toEqual([]);
    }
  });

  it("refuses full names outside an extension namespace, and a second definition of one name", () => {
    for (const name of ["seq.generator.mine", "core.x", "svs.x", "seqv.x", "plain", ".x", "x."]) {
      const { language, issues } = createLanguage([{ symbols: [ext(name, null, "generator")] }]);
      expect(
        issues.map((i) => i.code),
        name,
      ).toEqual(["seqql.symbol.namespace"]);
      expect(language.symbols.size).toBe(62);
    }
    const first = ext("molstar.select", null, "generator");
    const second = ext("molstar.select", "molql", "generator");
    const { language, issues } = createLanguage([{ symbols: [first] }, { symbols: [second] }]);
    expect(issues.map((i) => i.code)).toEqual(["seqql.symbol.duplicate"]);
    expect(language.symbols.get("molstar.select")).toBe(first);
    expect(language.aliases.has("molql")).toBe(false);
  });

  it("reports registration issues with no place", () => {
    const { issues } = createLanguage([{ symbols: [molql, other] }]);
    expect(issues[0]).toMatchObject({ severity: "error", path: "", at: null, open: null, signature: null });
  });
});

describe("the builder and the JSON form", () => {
  const kinases = Q.seq.filter.pick(Q.seq.generator.items({ of: "domains" }), {
    where: Q.core.str.matches(Q.seq.item.label(), "kinase"),
  });

  it("builds the stored form, character for character", () => {
    expect(JSON.stringify(kinases)).toBe(
      '{"head":"seq.filter.pick","args":{"0":{"head":"seq.generator.items","args":{"of":"domains"}},"where":{"head":"core.str.matches","args":{"0":{"head":"seq.item.label"},"1":"kinase"}}}}',
    );
  });

  it("builds the document expression", () => {
    expect(kinases).toEqual({
      head: "seq.filter.pick",
      args: {
        "0": { head: "seq.generator.items", args: { of: "domains" } },
        where: { head: "core.str.matches", args: { "0": { head: "seq.item.label" }, "1": "kinase" } },
      },
    });
  });

  it("keeps sel as an application with its string", () => {
    expect(JSON.stringify(Q.seq.macro.sel("1brs/A[10-20]"))).toBe(
      '{"head":"seq.macro.sel","args":{"0":"1brs/A[10-20]"}}',
    );
  });

  it("puts a rest positional's values in an array at its index", () => {
    expect(Q.seq.modifier.range("1brs/A", "10-20", 45, { numbering: "uniprot" })).toEqual({
      head: "seq.modifier.range",
      args: { "0": "1brs/A", "1": ["10-20", 45], numbering: "uniprot" },
    });
    expect(Q.seq.modifier.slice("1brs/A", 9, 20)).toEqual({
      head: "seq.modifier.slice",
      args: { "0": "1brs/A", "1": [9, 20] },
    });
    expect(Q.core.logic.and(Q.core.rel.gr(Q.seq.item.length(), 30), Q.core.rel.in(Q.seq.item.id(), "a", "b"))).toEqual({
      head: "core.logic.and",
      args: {
        "0": { head: "core.rel.gr", args: { "0": { head: "seq.item.length" }, "1": 30 } },
        "1": [{ head: "core.rel.in", args: { "0": { head: "seq.item.id" }, "1": ["a", "b"] } }],
      },
    });
  });

  it("omits defaults and absent keywords, and writes keywords in signature order", () => {
    expect(Q.seq.generator.systems({ kind: "unit" })).toEqual({ head: "seq.generator.systems" });
    expect(Q.seq.filter.first("domains", 1)).toEqual({ head: "seq.filter.first", args: { "0": "domains" } });
    expect(
      Q.seq.generator.search("GAATTC", { edits: 0, regex: false, ends: "local_best_last", in: undefined }),
    ).toEqual({
      head: "seq.generator.search",
      args: { "0": "GAATTC" },
    });
    const onto = Q.seq.modifier.onto("q", "msa", { keep: true, kind: "row" });
    expect(Object.keys(onto.args ?? {})).toEqual(["0", "1", "kind", "keep"]);
    expect(Q.seq.generator.all()).toEqual({ head: "seq.generator.all" });
  });

  it("throws on a keyword the symbol lacks or too many positionals (a programming error, not a query issue)", () => {
    expect(() => Q.seq.generator.items({ in: "domains" })).toThrow("no keyword :in");
    expect(() => Q.seq.generator.all("x")).toThrow("too many positional");
    expect(() => Q.seq.generator.items({ of: "a" }, "x")).toThrow("too many positional");
  });

  it("reaches every base symbol by its full name", () => {
    for (const s of BASE_SYMBOLS) {
      const node = s.name.split(".").reduce<unknown>((n, part) => (n as Record<string, unknown>)[part], Q);
      expect(node, s.name).toBe(s);
    }
    expect(Q.seq.position.value_of("plddt")).toEqual({ head: "seq.position.value_of", args: { "0": "plddt" } });
  });

  it("round-trips through JSON, and recognizes the JSON shape", () => {
    const exprs: Expression[] = [
      kinases,
      Q.seq.modifier.range("1brs/A", "10-20", 45),
      Q.seq.combinator.intersect("1brs/A", Q.seq.modifier.complement("domains")),
      Q.seq.modifier.onto(Q.seq.generator.input("selection"), "msa", { keep: true }),
      "1brs/A",
      12,
      true,
    ];
    for (const e of exprs) {
      const back: unknown = JSON.parse(JSON.stringify(e));
      expect(back).toEqual(e);
      expect(isExpression(back)).toBe(true);
    }
    for (const bad of [
      null,
      undefined,
      Number.NaN,
      Number.POSITIVE_INFINITY,
      [],
      {},
      { head: 1 },
      { head: "x", extra: 1 },
    ]) {
      expect(isExpression(bad), String(JSON.stringify(bad))).toBe(false);
    }
    expect(isExpression({ head: "x", args: [] })).toBe(false);
    expect(isExpression({ head: "x", args: { "0": [["nested"]] } })).toBe(false);
    expect(isExpression({ head: "x", args: { "0": null } })).toBe(false);
    expect(isExpression({ head: "x", args: { "0": [1, { head: "y" }] } })).toBe(true);
  });

  it("makes symbols frozen functions carrying their definition", () => {
    const s: SymbolDef = Q.seq.item.label;
    expect(typeof s).toBe("function");
    expect(s.name).toBe("seq.item.label");
    expect(Object.isFrozen(s)).toBe(true);
    expect(Object.keys(s)).toContain("name");
  });
});
