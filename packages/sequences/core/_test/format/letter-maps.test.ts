/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Letter maps (LM1; seqspec:format#layer-types, `letter_map`; seqspec:format#extensions; seqspec:format#letter-maps):
// the kernel's named tables, the 128-entry table a draw reads, the registry kind (an extension's map, namespaced), the
// validator's instructive errors for a name and for a dict's entries, and resolution handing viewers the table whether
// named or inline.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import {
  compile,
  KernelRegistry,
  letterMapTable,
  register,
  resolveView,
  type SeqViewSpecDoc,
  SvsRegistrationError,
  supportedAll,
  validate,
} from "../../src/format/index.js";
import { KernelLetterMaps } from "../../src/format/palettes/letter-maps.js";
import { fixtureLoaders, registries } from "./load-fixtures.js";

const reg = registries();
const kernel = (name: string) => KernelLetterMaps.find((m) => m.name === name)?.map ?? {};
const shown = (name: string, s: string) => [...s].map((c) => kernel(name)[c] ?? c).join("");

function doc(map: unknown, type: "letters" | "rows" = "letters"): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  const s = b.data.sequence({ ref: "s", alphabet: "dna", value: "ATGCNRY-" });
  const v = b.view({ ref: "v", label: "V" }).block({ sequence: s });
  v.track({ ref: "t" }).layer({ type, ref: "l", letter_map: map } as never);
  return b.getDocument();
}
const errs = (d: SeqViewSpecDoc) => validate(d, reg).filter((i) => i.severity === "error");

describe("the kernel's maps", () => {
  it("complement: IUPAC's, ambiguity codes and lower case included, unlisted letters kept", () => {
    expect(shown("complement", "ATGCNRYKMBVDHSW-.")).toBe("TACGNYRMKVBHDSW-.");
    expect(shown("complement", "atgcryu")).toBe("tacgyra");
  });
  it("complement_rna reads A as U; lower lowers", () => {
    expect(shown("complement_rna", "ATUGC")).toBe("UAACG");
    expect(shown("lower", "MKTAYI*-")).toBe("mktayi*-");
  });
  it("a complement of a complement is the original, but for U (which reads as T)", () => {
    for (const c of "ACGTRYKMBVDHSWN") expect(shown("complement", shown("complement", c))).toBe(c);
  });
});

describe("letterMapTable", () => {
  it("a flat 128-entry table; null when nothing changes", () => {
    const t = letterMapTable({ A: "T", c: "g" });
    expect(t?.length).toBe(128);
    expect([t?.[65], t?.[99], t?.[66]]).toEqual([84, 103, 66]);
    expect(letterMapTable(null)).toBeNull();
    expect(letterMapTable({})).toBeNull();
    expect(letterMapTable({ A: "A" })).toBeNull();
  });
  it("ignores entries a draw cannot index (the validator reports them)", () => {
    expect(letterMapTable({ AB: "T", A: "é", "\u0001": "A" })).toBeNull();
  });
});

describe("the validator", () => {
  it("accepts a name, a dict and null on letters and rows", () => {
    for (const type of ["letters", "rows"] as const)
      for (const m of ["complement", "lower", { A: "T", T: " " }, {}, null]) expect(errs(doc(m, type))).toEqual([]);
  });
  it("an unknown name lists the registered ones", () => {
    const [e, ...more] = errs(doc("complment"));
    expect(more).toEqual([]);
    expect([e.code, e.path]).toEqual([
      "svs.params.type",
      "/root/children/1/children/0/children/1/children/0/params/letter_map",
    ]);
    expect(e.message).toContain('"complment"');
    expect(e.message).toContain("complement, complement_rna, lower");
  });
  it("a dict entry says which and what is valid", () => {
    const d = errs(doc({ AB: "T", A: "TT", C: 5 }));
    expect(d.map((i) => i.path.split("/").slice(-2).join("/"))).toEqual([
      "letter_map/AB",
      "letter_map/A",
      "letter_map/C",
    ]);
    expect(d[0].message).toContain("key of one printable character");
    expect(d[1].message).toContain("one printable character (the letter shown)");
  });
  it("another type is a type error", () => {
    expect(errs(doc(5))[0].message).toContain("a letter map name");
    expect(errs(doc(["A"]))).toHaveLength(1);
  });
});

describe("the registry kind", () => {
  it("an extension's map is namespaced, and usable by name", () => {
    const r = registries();
    register(r, { letterMaps: [{ name: "org.example:purine", map: { A: "R", G: "R", C: "Y", T: "Y" } }] });
    expect(errs(doc("org.example:purine"))).not.toEqual([]);
    expect(validate(doc("org.example:purine"), r).filter((i) => i.severity === "error")).toEqual([]);
    expect(r.get("letterMaps", "complement")).toBeDefined();
  });
  it("a name that is not namespaced, or one the kernel has, is a conflict", () => {
    const r = registries();
    expect(() => register(r, { letterMaps: [{ name: "purine", map: {} }] })).toThrow(SvsRegistrationError);
    expect(() => register(r, { letterMaps: [{ name: "complement", map: {} }] })).toThrow(/already registered/);
    expect(KernelRegistry.letterMaps).toBe(KernelLetterMaps);
  });
});

describe("resolution", () => {
  async function resolved(d: SeqViewSpecDoc) {
    const { set } = await compile(d, reg, fixtureLoaders());
    const { view } = resolveView(d, set, "v", reg, supportedAll(reg));
    return view.layers.params[0].letter_map;
  }
  it("a named map is handed on as its table; a dict and null as written", async () => {
    expect(await resolved(doc("complement"))).toEqual(kernel("complement"));
    expect(await resolved(doc({ A: "T" }))).toEqual({ A: "T" });
    expect(await resolved(doc(null))).toBeNull();
    expect(await resolved(doc("lower", "rows"))).toEqual(kernel("lower"));
  });
  it("an unregistered name (reported statically) resolves to none", async () => {
    expect(await resolved(doc("nope"))).toBeNull();
  });

  /** Each layer's resolved `letter_map` by ref: a protein, its gene (`of`), a protein alignment and a DNA one. */
  async function byRef(map: unknown, r = reg) {
    const b = createSeqViewSpec();
    const g = b.data.sequence({ ref: "g", value: "ATGGCTTCCAAG" });
    const p = b.data.sequence({ ref: "p", value: "MKTAYIAKQR" });
    b.data.mapping({ ref: "m", from: p, to: g, method: "codon", cds: [[0, 12]] } as never);
    const pa = b.data.alignment({ ref: "pa", value: ["MKVLA", "MRVLG"] });
    const da = b.data.alignment({ ref: "da", value: ["ACGTA", "ACGAA"] });
    const v = b.view({ ref: "v", label: "V" });
    const prot = v.section({ ref: "v.p" }).block({ sequence: p });
    prot.track({ ref: "pt" }).layer({ type: "letters", ref: "p.letters", letter_map: map } as never);
    const gene = v.section({ ref: "v.g" }).block({ sequence: g });
    gene.track({ ref: "gt" }).layer({ type: "letters", ref: "g.letters", letter_map: map } as never);
    gene.track({ ref: "ot" }).layer({ type: "letters", ref: "g.of", of: "p", letter_map: map } as never);
    v.section({ ref: "v.pa" })
      .block({ alignment: pa })
      .track({ ref: "pr" })
      .layer({ type: "rows", ref: "pa.rows", letter_map: map } as never);
    v.section({ ref: "v.da" })
      .block({ alignment: da })
      .track({ ref: "dr" })
      .layer({ type: "rows", ref: "da.rows", letter_map: map } as never);
    const d = b.getDocument();
    const { set } = await compile(d, r, fixtureLoaders());
    const { view } = resolveView(d, set, "v", r, supportedAll(r));
    const out: Record<string, unknown> = {};
    for (let k = 0; k < view.layers.count; k++)
      out[set.strings.get(view.layers.ref[k])] = view.layers.params[k].letter_map;
    return out;
  }
  it("a map listing alphabets leaves another alphabet's letters as stored (a protein is not complemented)", async () => {
    const c = kernel("complement");
    expect(await byRef("complement")).toEqual({
      "p.letters": null,
      "g.letters": c,
      "g.of": null,
      "pa.rows": null,
      "da.rows": c,
    });
    // `lower` and a dict list none: every layer takes them.
    expect(Object.values(await byRef("lower")).every((m) => m === kernel("lower"))).toBe(true);
    expect(Object.values(await byRef({ K: "k" }))).toEqual(Array(5).fill({ K: "k" }));
    // An extension's map for proteins.
    const r = registries();
    register(r, { letterMaps: [{ name: "org.example:charge", map: { K: "+", D: "-" }, alphabets: ["protein"] }] });
    const charge = await byRef("org.example:charge", r);
    expect([charge["p.letters"], charge["g.of"], charge["pa.rows"], charge["g.letters"], charge["da.rows"]]).toEqual([
      { K: "+", D: "-" },
      { K: "+", D: "-" },
      { K: "+", D: "-" },
      null,
      null,
    ]);
  });
  it("emphasis takes the map of its track's letters or rows layer before it (its letters are the track's)", async () => {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: "s", alphabet: "dna", value: "ATGCNRY" });
    const f = b.data.annotation({ ref: "f", geometry: "ranges", target: s, value: { start: [1], stop: [3] } });
    const v = b.view({ ref: "v", label: "V" }).block({ sequence: s });
    v.track({ ref: "t1" }).letters({ ref: "l1", letter_map: "complement" }).emphasis({ ref: "e1", data: f });
    v.track({ ref: "t2" }).emphasis({ ref: "e0", data: f }).letters({ ref: "l2", letter_map: "lower" });
    v.track({ ref: "t3" })
      .letters({ ref: "l3", letter_map: { A: "x" } })
      .letters({ ref: "l4" })
      .emphasis({ ref: "e3", data: f });
    const d = b.getDocument();
    const { set } = await compile(d, reg, fixtureLoaders());
    const { view } = resolveView(d, set, "v", reg, supportedAll(reg));
    const at = (ref: string) => view.layers.params[view.layers.ref.indexOf(set.strings.find(ref))].letter_map;
    expect([at("e1"), at("e0"), at("e3")]).toEqual([kernel("complement"), undefined, undefined]);
  });
});
