/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// NP1b: computed alignments (`realign`, seqspec:model#alignments-and-trees) and the `codon_rows` mapping provider
// (seqspec:format#mappings): inputs read in the set, the aligned rows under the node's key, a mapping per row to its
// input, caching and aborting.

import { AlignerSizeError } from "@molstar/seq-core/algorithms";
import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { type Alignment, type CoordinateSystem, type Issue, KeyKind, type SeqSet } from "@molstar/seq-core/model";
import { beforeEach, describe, expect, it } from "vitest";
import { type AlignerEntry, compile, register, type SeqViewSpecDoc, validate } from "../../src/format/index.js";
import { clearNodeProviderCache } from "../../src/format/load/node-providers.js";
import { fixtureLoaders, registries } from "./load-fixtures.js";

const codes = (issues: Issue[]) => [...new Set(issues.map((i) => i.code))].sort();
const errors = (issues: Issue[]) => issues.filter((i) => i.severity === "error");
const rows = (a: Alignment) =>
  Array.from({ length: a.rowCount }, (_, r) => {
    const slot = a.rows.slotOf(r as never);
    return slot < 0
      ? null
      : String.fromCharCode(...a.rows.residues.subarray(a.rows.residueStart[slot], a.rows.residueStart[slot + 1]));
  });
const named = (set: SeqSet, key: string) =>
  set.alignments.find((a) => a && set.strings.get(a.key) === key) as Alignment;
/** A coordinate system by the name a document gives it. */
function nameOf(set: SeqSet, c: CoordinateSystem | undefined): string {
  if (!c) return "";
  if (c.kind === KeyKind.molecule) return set.strings.get(set.molecules[c.id]?.key ?? -1);
  if (c.kind === KeyKind.unit) return set.strings.get(set.units.key[c.id]);
  if (c.kind === KeyKind.alignment && c.row >= 0) return `${set.strings.get(set.alignments[c.id]?.key ?? -1)}/${c.row}`;
  return String(c.kind);
}
const mappingKeys = (set: SeqSet) => set.mappings.map((m) => set.strings.get(m?.key ?? -1)).sort();

const S1 = "MKVLAGWTQRSTEDNLK";
const S2 = "MKVLAGWTRRSTEDNLK";
const S3 = "KVLAGWTQRSTEDNLK";

type Build = (b: ReturnType<typeof createSeqViewSpec>) => unknown;
function doc(build: Build): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: "Realign" });
  build(b);
  const v = b.view({ ref: "v", label: "V", style: "msa" });
  v.block({ alignment: "re" });
  v.track({ label: "Rows" }).rows({});
  return b.getDocument();
}
/** Three sequences and `re` realigning them (or whatever `inputs` names). */
function threeSequences(re: Record<string, unknown> = {}, more: Build = () => {}): SeqViewSpecDoc {
  return doc((b) => {
    b.data.sequence({ ref: "s1", value: S1, label: "first" });
    b.data.sequence({ ref: "s2", value: S2 });
    b.data.sequence({ ref: "s3", value: S3, label: "third" });
    more(b);
    b.data.alignment({ ref: "re", provider: "realign", inputs: ["s1", "s2", "s3"], ...re } as never);
  });
}
const load = (d: SeqViewSpecDoc, signal?: AbortSignal, reg = registries()) => compile(d, reg, fixtureLoaders(), signal);

beforeEach(clearNodeProviderCache);

describe("realign", () => {
  it("aligns sequences into rows that keep their inputs' order, letters and labels, each row bound to its input", async () => {
    const { set, issues } = await load(threeSequences());
    expect(errors(issues)).toEqual([]);
    const a = named(set, "re");
    expect(a.rowCount).toBe(3);
    expect(rows(a)).toEqual([S1, S2, S3]); // each row's own residues
    expect([0, 1, 2].map((r) => set.strings.get(a.labelOf(r as never)))).toEqual(["first", "s2", "third"]);
    expect(a.columnCount).toBe(17);
    expect(mappingKeys(set)).toEqual(["re.row-0", "re.row-1", "re.row-2"]);
    expect(set.mappings.map((m) => [m?.method, nameOf(set, m?.from), nameOf(set, m?.to)])).toEqual([
      ["identity", "re/0", "s1"],
      ["identity", "re/1", "s2"],
      ["identity", "re/2", "s3"],
    ]);
  });
  it("puts a clipped sequence's first letter under the others' second", async () => {
    const a = named((await load(threeSequences())).set, "re");
    const onFirstColumn = [0, 1, 2].map((r) => a.rows.colToPos(r as never, 0 as never) >= 0);
    expect(onFirstColumn).toEqual([true, true, false]); // KVLAG... has no M: its first column is a gap
    expect(a.rows.colToPos(2 as never, 1 as never)).toBe(0);
  });
  it("reads units by their range, and rows of another alignment", async () => {
    const d = doc((b) => {
      const s = b.data.sequence({ ref: "s1", value: `GG${S1}GG` });
      b.data.unit({ ref: "u1", sequence: s, range: [2, 2 + S1.length] });
      b.data.sequence({ ref: "s2", value: S2 });
      b.data.alignment({ ref: "old", value: ["MKVL--AGWTQ", "MKVLEEAGWTR", "MKVL--AGWTQ"], labels: ["a", "b", "c"] });
      b.data.alignment({ ref: "re", provider: "realign", inputs: ["u1", "s2", "old/1"] } as never);
    });
    const { set, issues } = await load(d);
    expect(errors(issues)).toEqual([]);
    expect(rows(named(set, "re"))).toEqual([S1, S2, "MKVLEEAGWTR"]);
    expect(set.mappings.map((m) => nameOf(set, m?.to))).toEqual(["u1", "s2", "old/1"]);
  });
  it("realigns a whole alignment: its rows ungapped, in order, with its labels and reference row, and a tree over them", async () => {
    const d = doc((b) => {
      b.data.alignment({
        ref: "old",
        value: ["MKV--LAGWTQ", "MK-VLLAGWTR", "MKVL--AGWTQ", "MKVLAGWTQ--"],
        labels: ["a", "b", "c", "d"],
        reference_row: 2,
      });
      const re = b.data.alignment({ ref: "re", provider: "realign", inputs: "old" } as never);
      re.tree({ provider: "nj" } as never);
    });
    const { set, issues } = await load(d);
    expect(errors(issues)).toEqual([]);
    const a = named(set, "re");
    expect(rows(a)).toEqual(rows(named(set, "old"))); // the same residues, realigned
    expect([0, 1, 2, 3].map((r) => set.strings.get(a.labelOf(r as never)))).toEqual(["a", "b", "c", "d"]);
    expect(a.referenceRow).toBe(2);
    expect(a.tree?.leafOrder.length).toBe(4);
    expect(mappingKeys(set)).toEqual(["re.row-0", "re.row-1", "re.row-2", "re.row-3"]);
  });
});

describe("realign: a missing piece or a bad request", () => {
  it("is validated: inputs with a provider, a provider with inputs", () => {
    const reg = registries();
    const noInputs = doc((b) => {
      b.data.sequence({ ref: "s1", value: S1 });
      b.data.alignment({ ref: "re", provider: "realign" } as never);
    });
    expect(codes(validate(noInputs, reg))).toEqual(["svs.params.missing"]);
    const noProvider = doc((b) => {
      b.data.sequence({ ref: "s1", value: S1 });
      b.data.alignment({ ref: "re", value: ["MK"], inputs: ["s1"] } as never);
    });
    expect(codes(validate(noProvider, reg))).toContain("svs.params.exclusive");
    expect(codes(validate(threeSequences({ options: { aligner: 3 } }), reg))).toEqual(["svs.params.type"]);
  });
  it("leaves an unregistered provider, or one that computes another kind, an unavailable alignment", async () => {
    for (const provider of ["nope", "nj", "codon_rows"]) {
      const { set, issues } = await load(threeSequences({ provider }));
      expect(codes(issues).filter((c) => c !== "svs.layer.unplaced")).toEqual(["svs.provider.unavailable"]);
      expect(set.alignments.filter(Boolean)).toHaveLength(0);
    }
  });
  it("names what cannot be aligned: a ligand, a name that is nothing, the alignment itself", async () => {
    const lig = doc((b) => {
      b.data.sequence({ ref: "s1", value: S1 });
      const l = b.data.ligand({ ref: "l1", value: ["C", "C", "O"] });
      b.data.unit({ ref: "ul", ligand: l });
      b.data.alignment({ ref: "re", provider: "realign", inputs: ["s1", "ul"] } as never);
    });
    expect(codes((await load(lig)).issues)).toContain("svs.ligand.unsupported");
    const none = threeSequences({ inputs: ["s1", "s2", "s3", "s2/4"] });
    expect(codes(errors(validate(none, registries())))).toEqual(["svs.ref.wrong_kind"]);
    const loop = doc((b) => {
      b.data.sequence({ ref: "s1", value: S1 });
      b.data.alignment({ ref: "re", provider: "realign", inputs: ["s1", "re/0"] } as never);
    });
    expect(codes(errors((await load(loop)).issues)).length).toBeGreaterThan(0);
  });
  it("reports options the aligner refuses as svs.params.type, and an aligner that aligns two only", async () => {
    const bad = await load(threeSequences({ options: { aligner_options: { gap_open: -1 } } }));
    expect(bad.issues.filter((i) => i.code === "svs.params.type").map((i) => i.path)).toEqual([
      "/root/children/0/children/3/params/options/aligner_options",
    ]);
    const two = await load(threeSequences({ options: { aligner: "needleman_wunsch" } }));
    expect(two.issues.map((i) => i.message).join()).toContain("aligns two sequences only");
    expect(named(two.set, "re")).toBeUndefined();
  });
  it("is svs.provider.size when the aligner's full matrix is too large", async () => {
    const reg = registries();
    const stub: AlignerEntry = {
      name: "x:big",
      options: (reg.get("aligners", "biomsa") as AlignerEntry).options,
      align: async () => null,
      alignMany: async () => {
        throw new AlignerSizeError("x:big", 9000, 9000);
      },
    };
    register(reg, { aligners: [stub] });
    const { set, issues } = await load(threeSequences({ options: { aligner: "x:big" } }), undefined, reg);
    expect(codes(issues).filter((c) => c !== "svs.layer.unplaced")).toEqual(["svs.provider.size"]);
    expect(named(set, "re")).toBeUndefined();
  });
});

describe("realign: caching and aborting", () => {
  it("asks the aligner once for the same inputs and options, and again when a sequence or an option changes", async () => {
    const reg = registries();
    const biomsa = reg.get("aligners", "biomsa") as AlignerEntry;
    let runs = 0;
    register(reg, {
      aligners: [
        {
          ...biomsa,
          name: "x:counted",
          alignMany: (...a: Parameters<NonNullable<AlignerEntry["alignMany"]>>) => {
            runs++;
            return (biomsa.alignMany as NonNullable<typeof biomsa.alignMany>)(...a);
          },
        },
      ],
    });
    const run = (re = {}, more: Build = () => {}) =>
      load(threeSequences({ options: { aligner: "x:counted" }, ...re }, more), undefined, reg);
    const first = await run();
    await run();
    expect(runs).toBe(1);
    expect(rows(named((await run()).set, "re"))).toEqual(rows(named(first.set, "re")));
    await run({ options: { aligner: "x:counted", aligner_options: { gap_open: 500 } } });
    expect(runs).toBe(2);
  });
  it("rejects a compile with the signal's reason when a newer request aborts a long run, and keeps nothing", async () => {
    const AA = "ACDEFGHIKLMNPQRSTVWY";
    let seed = 9;
    const rnd = () => {
      seed = (Math.imul(seed, 1103515245) + 12345) & 0x7fffffff;
      return seed / 0x7fffffff;
    };
    const base = Array.from({ length: 300 }, () => AA[Math.floor(rnd() * 20)]).join("");
    const many = Array.from({ length: 120 }, () =>
      [...base].map((c) => (rnd() < 0.2 ? AA[Math.floor(rnd() * 20)] : c)).join(""),
    );
    const d = doc((b) => {
      for (const [k, v] of many.entries()) b.data.sequence({ ref: `q${k}`, value: v });
      b.data.alignment({ ref: "re", provider: "realign", inputs: many.map((_, k) => `q${k}`) } as never);
    });
    const c = new AbortController();
    const running = load(d, c.signal);
    setTimeout(() => c.abort(new Error("newer request")), 1);
    await expect(running).rejects.toThrow("newer request");
  });
});

describe("codon_rows", () => {
  const PROTEIN = ["MAS", "MKV"];
  const CDNA = ["ATGGCTTCC", "ATGAAAGTT"];
  function pair(mapping: Record<string, unknown> = {}, cdna = CDNA): SeqViewSpecDoc {
    const b = createSeqViewSpec({ title: "Coding rows" });
    b.data.alignment({ ref: "prot", value: PROTEIN, labels: ["p1", "p2"], alphabet: "protein" });
    b.data.alignment({ ref: "dna", value: cdna, labels: ["d1", "d2", "d3"].slice(0, cdna.length), alphabet: "dna" });
    b.data.mapping({
      ref: "cr",
      from: "prot",
      to: "dna",
      method: "provider",
      provider: "codon_rows",
      ...mapping,
    } as never);
    const v = b.view({ ref: "v", label: "V" });
    v.block({ alignment: "prot" });
    v.track({ label: "Rows" }).rows({});
    return b.getDocument();
  }
  it("pairs row k with row k by a codon mapping each, every status exact", async () => {
    const { set, issues } = await load(pair());
    expect(errors(issues)).toEqual([]);
    expect(mappingKeys(set)).toEqual(["cr.row-0", "cr.row-1"]);
    expect(set.mappings.map((m) => [m?.method, nameOf(set, m?.from), nameOf(set, m?.to), m?.span])).toEqual([
      ["codon", "prot/0", "dna/0", 3],
      ["codon", "prot/1", "dna/1", 3],
    ]);
  });
  it("takes the pairs and the genetic code as options", async () => {
    const d = pair({ options: { rows: [[1, 2]], genetic_code: 11 } }, [...CDNA, "ATGAAAGTT"]);
    const { set, issues } = await load(d);
    expect(errors(issues)).toEqual([]);
    expect(set.mappings.map((m) => [set.strings.get(m?.key ?? -1), nameOf(set, m?.to), m?.codon?.geneticCode])).toEqual(
      [["cr.row-1", "dna/2", 11]],
    );
  });
  it("says what to write when the rows cannot be paired", async () => {
    const counts = await load(pair({}, [...CDNA, "ATGAAAGTT"]));
    expect(counts.issues.map((i) => [i.code, i.path, i.message.includes("give rows")])).toEqual([
      ["svs.mapping.codon", "/root/children/0/children/2/params/options/rows", true],
    ]);
    const range = await load(pair({ options: { rows: [[0, 5]] } }));
    expect(range.issues.map((i) => i.code)).toEqual(["svs.mapping.codon"]);
    expect(range.set.mappings.filter(Boolean)).toHaveLength(0);
    const notAln = await load(pair({ to: "prot/1" }));
    expect(codes(notAln.issues)).toContain("svs.ref.derived_unknown");
  });
  it("is checked like a codon mapping: a protein row that the codons do not translate is a conflict, not an error", async () => {
    const { set, issues } = await load(pair({}, ["ATGGCTTCC", "ATGAAAAAA"]));
    expect(errors(issues)).toEqual([]);
    expect(set.mappings.filter(Boolean)).toHaveLength(2);
  });
});
