/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The builder's own tests (moved from the format's compile.test.ts in B1), its copied helpers, and its dependency
// rule: no runtime @molstar/seq-* module (seqspec:format#builders).

import {
  createRegistries,
  fromJSON,
  fromMsgPack,
  KernelRegistry,
  register,
  type SeqViewSpecDoc,
  toJSON,
  toMsgPack,
  validate,
} from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { basePairs, createSeqViewSpec, translate, VARIANT_ROWS } from "../../src/builder/index.ts";

const nodeRefs = (doc: SeqViewSpecDoc) => {
  const out: string[] = [];
  const walk = (n: SeqViewSpecDoc["root"]) => {
    if (n.ref) out.push(n.ref);
    n.children?.forEach(walk);
  };
  walk(doc.root);
  return out;
};

describe("generated refs", () => {
  it("avoid explicit refs and number repeats in building order; equal sources are one node", () => {
    const b = createSeqViewSpec();
    const src = { host_file: "x.fasta", format: "fasta" };
    b.data.sequence({ ref: "sequence", value: "MK" });
    b.data.sequence({ source: src, record: "a" });
    b.data.sequence({ source: src, record: "b" });
    const v = b.view({ ref: "v", label: "V" });
    v.track({ label: "A b!" }).letters();
    v.track({ label: "A b" }).letters();
    v.track({}).layer({ type: "org.example:dots", fallback: "points" });
    const doc = b.getDocument();
    expect(nodeRefs(doc)).toEqual([
      "sequence",
      "sequence-2.src",
      "sequence-2",
      "sequence-3",
      "v",
      "v.main",
      "v.a-b",
      "v.letters",
      "v.a-b-2",
      "v.letters-2",
      "v.org-example-dots",
      "v.org-example-dots-2",
    ]);
  });
});

describe("the sugar", () => {
  it("offset and gapped pairs become segments and a CIGAR; columnar values are dictionary-coded", () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", value: "MKVL" });
    const c = b.data.sequence({ ref: "c", value: "AMKL" });
    b.data.mapping({ ref: "m1", from: a, to: c, method: "offset", offset: 1 });
    b.data.mapping({ ref: "m2", from: a, to: c, method: "pairwise", gapped: ["-MKVL", "AMK-L"] });
    b.data.annotation({
      ref: "t",
      geometry: "ranges",
      target: a,
      columnar: true,
      value: [
        { start: 0, label: "x" },
        { start: 1, label: "x" },
        { start: 2, label: "y" },
      ],
    });
    const data = b.getDocument().root.children?.[0].children ?? [];
    expect(data[2].params).toMatchObject({ method: "segments", segments: [[0, 1, null]] });
    expect(data[3].params).toMatchObject({ method: "pairwise", cigar: "1D2M1I1M" });
    expect(data[4].params?.value).toEqual({ start: [0, 1, 2], label: { dict: ["x", "y"], codes: [0, 0, 1] } });
  });

  it("columnar dictionary coding is one pass: null-only, mixed and unique columns stay plain", () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", value: "MKVL" });
    const rows = [
      { start: 0, none: null, uniq: "p", mix: "x", nul: "u", n: 1 },
      { start: 1, none: null, uniq: "q", mix: 2, nul: null, n: 2 },
      { start: 2, none: null, uniq: "r", mix: "x", nul: "u", n: 3 },
      { start: 3, none: null, uniq: null, mix: "x", nul: "u", n: 4 },
    ];
    b.data.annotation({ ref: "t", geometry: "ranges", target: a, columnar: true, value: rows });
    const cols = b.getDocument().root.children?.[0].children?.[1].params?.value as Record<string, unknown>;
    expect(cols.none).toEqual([null, null, null, null]);
    expect(cols.uniq).toEqual(["p", "q", "r", null]); // every string distinct: nothing to gain
    expect(cols.mix).toEqual(["x", 2, "x", "x"]);
    expect(cols.nul).toEqual({ dict: ["u"], codes: [0, null, 0, 0] });
    expect(cols.n).toEqual([1, 2, 3, 4]);
  });

  it("columnar coding of 200,000 rows with 100,000 distinct strings does not scan a dictionary per row", () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", value: "MK" });
    const value = Array.from({ length: 200_000 }, (_, i) => ({ start: i, label: `f${i % 100_000}` }));
    b.data.annotation({ ref: "t", geometry: "ranges", target: a, columnar: true, value });
    const table = b.getDocument().root.children?.[0].children?.[1].params?.value as { label: { dict: string[] } };
    expect(table.label.dict.length).toBe(100_000);
  });

  it("split: a field's name is { by, order: null }; alt on variants is UniProt's grid; objects as given", () => {
    const b = createSeqViewSpec();
    const t = b.view({ ref: "v", label: "V" }).track({ ref: "t" });
    t.ranges({ split: "category" })
      .variants({ split: "alt" })
      .points({ split: { by: "x", order: ["b"] } });
    const layers = b.getDocument().root.children?.[1].children?.[0].children?.[0].children ?? [];
    expect(layers.map((l) => l.params?.split)).toEqual([
      { by: "category", order: null },
      { by: "alt", order: VARIANT_ROWS },
      { by: "x", order: ["b"] },
    ]);
    expect(VARIANT_ROWS.join("")).toBe("GAVLISTCMDNEQRKHFYWP-*");
  });
});

describe("a section's style", () => {
  it("an inline preset is written as given and round-trips through JSON and MessagePack", () => {
    const reg = createRegistries();
    register(reg, [KernelRegistry]);
    const b = createSeqViewSpec();
    b.data.sequence({ ref: "s", value: "MEEPQSDPSV" });
    const style = { preset: "msa", column_heads: "own_row", data_floor: 0, row_height: 18 } as const;
    b.view({ ref: "v", label: "V", style }).track({}).letters();
    b.view({ ref: "w", label: "W", style: "compact" }).track({}).letters();
    const doc = b.getDocument();
    const section = (d: SeqViewSpecDoc, view: number) => d.root.children?.[view + 1].children?.[0].params;
    expect(section(doc, 0)?.style).toEqual(style);
    expect(section(doc, 1)?.style).toBe("compact");
    expect(validate(doc, reg)).toEqual([]);
    for (const read of [fromJSON(toJSON(doc)), fromMsgPack(toMsgPack(doc))]) {
      expect(read.issues).toEqual([]);
      expect(read.doc).toEqual(doc);
      expect(section(read.doc as SeqViewSpecDoc, 0)?.style).toEqual(style);
    }
  });
});

describe("typed arrays", () => {
  const value = Float64Array.of(1, Number.NaN, 2, 3);
  const start = Int32Array.of(0, 2);
  const codes = Uint8Array.of(1, 0);

  it("table columns keep lossless typed arrays as given, typed codes too; elsewhere they become arrays", () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", value: "MKVL" });
    b.data.annotation({ ref: "v", geometry: "values", target: a, value: { value } });
    const ranges = { start, stop: Uint16Array.of(1, 3), label: { dict: ["x", "y"], codes } };
    const v = b.view({ ref: "v1", label: "V" });
    v.track({}).ranges({ annotation: { geometry: "ranges", target: a, value: ranges } });
    v.track({}).values({ data: "v", domain: Float64Array.of(0, 3) as never });
    const doc = b.getDocument();
    const [, vals] = doc.root.children?.[0].children ?? [];
    expect(vals.params?.value).toEqual({ value });
    expect((vals.params?.value as { value: unknown } | undefined)?.value).toBe(value);
    const [own, vl] = doc.root.children?.[1].children?.[0].children?.flatMap((t) => t.children ?? []) ?? [];
    const cols = own.children?.[0].params?.value as typeof ranges;
    expect([cols.start, cols.stop.constructor, cols.label.codes]).toEqual([start, Uint16Array, codes]);
    expect(cols.label.codes).toBe(codes);
    expect(vl.params?.domain).toEqual([0, 3]);
  });

  it("a Float32Array column is the producer's choice, kept; `float32` makes fractional columns Float32", () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", value: "MKVL" });
    const score = Float32Array.of(0.15, 0.3);
    b.data.annotation({ ref: "f", geometry: "values", target: a, value: { value: score } });
    const value = { index: [0, 1, 2], value: [0.15, null, 2], label: ["x", "y", "z"] };
    // `float32` is advanced and untyped: passed through a variable, not a literal.
    const g32 = { ref: "g", geometry: "values", target: a, value, float32: true } as const;
    b.data.annotation(g32);
    b.view({ ref: "v", label: "V" })
      .track({})
      .values({ data: "f", domain: Float32Array.of(0.15, 1) as never });
    const doc = b.getDocument();
    const [, f, g] = doc.root.children?.[0].children ?? [];
    expect((f.params?.value as { value: unknown } | undefined)?.value).toBe(score);
    const cols = g.params?.value as Record<string, unknown>;
    expect(g.params?.float32).toBeUndefined();
    expect([cols.index, cols.value, cols.label]).toEqual([
      [0, 1, 2],
      Float32Array.of(0.15, Number.NaN, 2),
      value.label,
    ]);
    expect(doc.root.children?.[1].children?.[0].children?.[0].children?.[0].params?.domain).toEqual([0.15, 1]);
    const bad = { geometry: "values", target: a, value: { value: new Uint8ClampedArray(2) } } as const;
    expect(() => b.data.annotation({ ref: "h", ...bad })).toThrow(
      'annotation "h": column "value" is a Uint8ClampedArray; give an array, or Int8-32, Uint8-32, Float32 or Float64.',
    );
    const t = b.view({ ref: "w", label: "W" }).track({});
    const codes = {
      geometry: "ranges",
      target: a,
      value: { start, label: { dict: ["x"], codes: new BigInt64Array(2) } },
    };
    expect(() => t.ranges({ annotation: codes as never })).toThrow(/column "label" is a BigInt64Array/);
  });
});

describe("a grid's typed rows", () => {
  it("are kept as given by the builder and written as JSON arrays only when serialized", () => {
    const b = createSeqViewSpec();
    b.data.sequence({ ref: "s", value: "MK" });
    const row = Float32Array.of(0.25, Number.NaN, 2, 3.5);
    b.data.annotation({ ref: "g", geometry: "grid", value: { target: ["s"], b_target: ["s"], value: [row] } });
    const doc = b.getDocument();
    const cols = doc.root.children?.[0].children?.[1].params?.value as { value: unknown[] };
    expect(cols.value[0]).toBe(row);
    expect(JSON.parse(toJSON(doc)).root.children[0].children[1].params.value.value).toEqual([[0.25, null, 2, 3.5]]);
  });
});

describe("base pairs", () => {
  it("brackets and letter pairs pair; unbalanced WUSS throws", () => {
    const b = createSeqViewSpec();
    const t = b.data.sequence({ ref: "t", value: "GGAACCUU" });
    basePairs(b.data, { ref: "bp", target: t, structure: "<<A..a>>" });
    const bp = b.getDocument().root.children?.[0].children?.[1];
    expect(bp?.params?.value).toEqual({ start: [0, 1, 2], b_start: [7, 6, 5], category: Array(3).fill("base_pair") });
    expect(() => basePairs(b.data, { target: t, structure: "((.)" })).toThrow(
      "basePairs: WUSS: the pair opened at column 0 is not closed.",
    );
    expect(() => basePairs(b.data, { target: t, structure: ".)" })).toThrow('WUSS: ")" at column 1 closes nothing.');
  });
});

describe("translate", () => {
  // case 14's gene: its CDS spliced from two exons, ending in TAA.
  const GENE =
    "GAATTCGCCACCATGGCTTCCAAGGAACTGACCCCACGTGACGGCCGTAAGTCTTTCCTAACTTCTGTTTCCACAGATTGGTACTTCAACAGCCAGGTGGAGTAAGCGGCCGCATAAAGT";
  const CDS: [number, number][] = [
    [12, 46],
    [76, 105],
  ];
  const nodes = (doc: SeqViewSpecDoc) => doc.root.children?.[0].children ?? [];
  const revcomp = (s: string) =>
    [...s]
      .reverse()
      .map((c) => ({ A: "T", T: "A", C: "G", G: "C" })[c] ?? c)
      .join("");
  const protein = (dna: string, opts: Omit<Parameters<typeof translate>[1], "from" | "ref"> = {}) => {
    const b = createSeqViewSpec();
    translate(b.data, { ref: "p", from: b.data.sequence({ ref: "g", value: dna }), ...opts });
    return nodes(b.getDocument())[1].params?.value;
  };

  it("writes case 14's protein and its codon mapping as ordinary nodes", () => {
    const b = createSeqViewSpec();
    const gene = b.data.sequence({ ref: "gene", label: "Gene fragment", alphabet: "dna", value: GENE });
    const p = translate(b.data, { ref: "protein", label: "Protein", from: gene, cds: CDS });
    expect(p.ref).toBe("protein");
    expect(nodes(b.getDocument()).slice(1)).toEqual([
      {
        kind: "sequence",
        ref: "protein",
        params: {
          value: "MASKELTPRDGHWYFNSQVE",
          alphabet: "protein",
          label: "Protein",
          provenance: {
            source: "@molstar/seq-core/builder",
            method: "translation",
            parameters: { from: "gene", cds: CDS, strand: "+", genetic_code: 1 },
          },
        },
      },
      { kind: "mapping", ref: "protein.codon", params: { from: "protein", to: "gene", method: "codon", cds: CDS } },
    ]);
  });

  it("reads to the first stop without cds, the reverse complement with it, starts and ambiguity codes by the table", () => {
    // Unspliced, read on through the intron to the TAA of exon 2.
    expect(protein(GENE, { start: 12 })).toBe("MASKELTPRDGRKSFLTSVSTDWYFNSQVE");
    const minus = CDS.map(([s, t]) => [120 - t, 120 - s] as [number, number]);
    expect(protein(revcomp(GENE), { cds: minus, strand: "-" })).toBe("MASKELTPRDGHWYFNSQVE");
    expect(protein("GTGGTGTGAAANTAA", { cds: [[0, 15]], genetic_code: 11 })).toBe("MV*X");
    expect(protein("GTGGTGTGAAANTAA", { cds: [[0, 15]] })).toBe("VV*X");
    expect(protein("ATGAGATGATAA", { cds: [[0, 12]], genetic_code: 2 })).toBe("M*W"); // AGA a stop, TGA W
    expect(protein("ATGAGATGATAA", { genetic_code: 2 })).toBe("M");
  });

  it("throws on a sequence it does not hold by value, and on input the mapping would reject", () => {
    const b = createSeqViewSpec();
    const read = b.data.sequence({ source: { host_file: "g.fasta", format: "fasta" } });
    expect(() => translate(b.data, { from: read })).toThrow(/by value/);
    const g = b.data.sequence({ value: GENE });
    expect(() => translate(b.data, { from: g, cds: [[12, 47]] })).toThrow(/multiple of 3/);
    expect(() => translate(b.data, { from: g, strand: "-" })).toThrow(/cds/);
    expect(() => translate(b.data, { from: g, genetic_code: 7 })).toThrow(/genetic code 7/);
  });

  it("throws where the codon mapping would: ranges out of order or overlapping, a bad alphabet, an empty translation", () => {
    const b = createSeqViewSpec();
    const g = b.data.sequence({ ref: "g", value: "ATGGCCTAAATGGCCTAA" });
    const t = (opts: Partial<Parameters<typeof translate>[1]>) => () => translate(b.data, { from: g, ...opts });
    // On +, ranges ascend and do not overlap; on -, they descend (seqspec:model#codon-maps).
    expect(
      t({
        cds: [
          [6, 9],
          [0, 6],
        ],
      }),
    ).toThrow("range 1 overlaps range 0 or is out of ascending order.");
    expect(
      t({
        cds: [
          [0, 6],
          [3, 9],
        ],
      }),
    ).toThrow("range 1 overlaps range 0 or is out of ascending order.");
    expect(
      t({
        cds: [
          [0, 6],
          [6, 12],
        ],
      }),
    ).not.toThrow();
    expect(
      t({
        cds: [
          [0, 6],
          [6, 12],
        ],
        strand: "-",
      }),
    ).toThrow("out of descending (strand -) order.");
    expect(
      t({
        cds: [
          [6, 12],
          [0, 6],
        ],
        strand: "-",
      }),
    ).not.toThrow();
    expect(t({ cds: [[0, 1.5]] })).toThrow(/outside the sequence/);
    expect(t({ start: 0.5 })).toThrow(/outside the sequence/);
    // A first codon that is a stop leaves nothing.
    const stop = b.data.sequence({ value: "TAAGGG" });
    expect(() => translate(b.data, { from: stop })).toThrow(/translation is empty/);
    expect(() => translate(b.data, { from: stop, cds: [[0, 3]] })).toThrow(/translation is empty/);
    // The alphabet: given, it must be dna or rna and hold its letters; else the model's detection decides.
    const prot = b.data.sequence({ value: "MKVLAAGIVG", alphabet: "protein" });
    expect(() => translate(b.data, { from: prot })).toThrow(/must be DNA or RNA; its alphabet is protein/);
    const guess = b.data.sequence({ value: "MKVLAAGIVGMK" });
    expect(() => translate(b.data, { from: guess })).toThrow(/must be DNA or RNA/);
    const rna = b.data.sequence({ value: "AUGGCC", alphabet: "rna" });
    expect(() => translate(b.data, { from: rna })).not.toThrow();
    const wrong = b.data.sequence({ value: "ATGGCC", alphabet: "rna" });
    expect(() => translate(b.data, { from: wrong })).toThrow(/not one of rna's/);
    const mixed = b.data.sequence({ value: "ATGUCC" });
    expect(() => translate(b.data, { from: mixed })).toThrow(/not one of dna's/);
    const lower = b.data.sequence({ value: "atggccatgaaatgcctgcatr" });
    expect(() => translate(b.data, { from: lower })).not.toThrow();
  });
});

describe("presets", () => {
  it("multiRegion: a block per region with its range, a sequence's own unit or a unit; valid as written", () => {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: "s", value: "MKVLAAGIVGLLLAQRSTVW" });
    const u = b.data.unit({ ref: "u", sequence: s });
    const v = b.view({ ref: "v", label: "V" });
    b.preset.multiRegion(
      v,
      [
        [s, 0, 5],
        [u, 10, 15],
      ],
      { direction: -1 },
    );
    v.track({ label: "L" }).letters();
    const doc = b.getDocument();
    const blocks = doc.root.children?.[1].children?.[0].children?.filter((c) => c.kind === "block");
    expect(blocks?.map((c) => c.params)).toEqual([
      { sequence: "s", range: [0, 5], direction: -1 },
      { unit: "u", range: [10, 15], direction: -1 },
    ]);
    const reg = createRegistries();
    register(reg, [KernelRegistry]);
    expect(validate(doc, reg)).toEqual([]);
  });

  it("reverseComplement: direction -1 and a complement map on the section's own letters; exonsOnly: positions", () => {
    const b = createSeqViewSpec();
    const g = b.data.sequence({ ref: "g", alphabet: "dna", value: "GAATTCGCCACCATGGCTTCCAAGGAA" });
    const p = translate(b.data, { ref: "p", from: g, cds: [[12, 27]] });
    const exons = b.data.annotation({
      ref: "exons",
      geometry: "ranges",
      target: g,
      value: { start: [0, 12], stop: [6, 27] },
    });
    const v = b.view({ ref: "v", label: "V" });
    const s = v.section({ ref: "v.rc" }).block({ sequence: g });
    s.track({ label: "Gene" }).letters();
    s.track({ label: "Own map" }).letters({ letter_map: "lower" });
    s.track({ label: "Protein" }).letters({ of: p });
    b.preset.reverseComplement(s);
    const e = b.preset.exonsOnly(v.section({ ref: "v.ex" }).block({ sequence: g }), exons);
    e.track().letters();
    // Before its tracks: the map is written on the layers the section has when the document is made.
    const whole = b.preset.reverseComplement(v.section({ ref: "v.main2" }).block({ sequence: g }), {
      map: "complement_rna",
    });
    whole.track({ ref: "t" }).layer({ type: "rows", ref: "r" });
    const doc = b.getDocument();
    const find = (kind: string, ref: string) => {
      const out: SeqViewSpecDoc["root"][] = [];
      const walk = (n: SeqViewSpecDoc["root"]) => {
        if (n.kind === kind && n.ref === ref) out.push(n);
        n.children?.forEach(walk);
      };
      walk(doc.root);
      return out[0];
    };
    expect(find("section", "v.rc")?.params).toEqual({ direction: -1 });
    const maps = find("section", "v.rc")?.children?.flatMap((t) => (t.kind === "track" ? (t.children ?? []) : []));
    expect(maps?.map((l) => l.params?.letter_map)).toEqual(["complement", "lower", undefined]);
    expect(find("section", "v.ex")?.params).toEqual({ positions: "exons" });
    expect(find("layer", "r")?.params).toEqual({ type: "rows", letter_map: "complement_rna" });
    expect(find("section", "v.main2")?.params).toEqual({ direction: -1 });
    const reg = createRegistries();
    register(reg, [KernelRegistry]);
    expect(validate(doc, reg)).toEqual([]);
  });

  it("plasmidMap and circularGenome: a circular section and a linear one following it; valid as written", () => {
    const b = createSeqViewSpec();
    const p = b.data.sequence({ ref: "p", alphabet: "dna", cyclic: true, value: "ACGT".repeat(50) });
    const at = (start: number[], stop: number[], more = {}) => ({ start, stop, ...more });
    const f = b.data.annotation({
      ref: "f",
      geometry: "ranges",
      target: p,
      value: at([10, 150], [60, 20], { strand: ["+", "-"], label: ["a", "b"] }),
    } as never);
    const pp = b.data.annotation({
      ref: "pp",
      geometry: "pairs",
      target: p,
      value: { start: [5], b_target: ["p"], b_start: [105] },
    } as never);
    const sites = b.data.annotation({
      ref: "sites",
      geometry: "ranges",
      target: p,
      value: { boundary: [40, 90], label: ["A", "B"] },
    } as never);
    const map = b.preset.plasmidMap(b.view({ ref: "m", label: "M" }), p, {
      ref: "m.ring",
      features: f,
      sites,
      pairs: pp,
    });
    const genome = b.preset.circularGenome(b.view({ ref: "g", label: "G" }), p, { features: f, detail: false });
    expect(genome.detail).toBe(null);
    const doc = b.getDocument();
    const sections = (view: string) =>
      doc.root.children?.find((c) => c.ref === view)?.children?.filter((c) => c.kind === "section") ?? [];
    const [ring, detail] = sections("m");
    expect([ring.params, detail.params]).toEqual([
      { label: "Plasmid map", style: "circular" },
      { label: "Plasmid map, linear", follow: "m.ring" },
    ]);
    const layers = (sec: (typeof ring)[]) =>
      sec.flatMap((s) => s.children?.flatMap((t) => (t.kind === "track" ? (t.children ?? []) : [])) ?? []);
    // On the ring the sites are ticks on the features' track (CG3); linearly a track of their own.
    expect(layers([ring]).map((l) => l.params?.type)).toEqual(["ruler", "ranges", "ranges", "pairs"]);
    expect(layers([ring])[1].params).toMatchObject({ labels: true, direction_field: "strand" });
    expect(layers([ring])[2].params).toMatchObject({
      data: "sites",
      labels: true,
      display: "union",
      color: { theme: "rules", params: { rules: [] } },
    });
    const tracks = (sec: typeof ring) => sec.children?.filter((t) => t.kind === "track").length;
    expect([tracks(ring), tracks(detail)]).toEqual([3, 5]);
    expect(map.detail).not.toBe(null);
    const [g] = sections("g");
    expect(layers([g]).map((l) => [l.params?.type, (l.params?.filter as { values?: string[] })?.values])).toEqual([
      ["ruler", undefined],
      ["ranges", ["+"]],
      ["ranges", ["-"]],
    ]);
    const reg = createRegistries();
    register(reg, [KernelRegistry]);
    expect(validate(doc, reg)).toEqual([]);
  });

  it("sashimi: a pairs track per strand, above and below their line, stroked and labeled by the count", () => {
    const b = createSeqViewSpec();
    const g = b.data.sequence({ ref: "g", alphabet: "dna", value: "ACGT".repeat(30) });
    const j = b.data.annotation({
      ref: "j",
      geometry: "pairs",
      target: g,
      value: { start: [10, 40], b_target: ["g", "g"], b_start: [30, 90], count: [12, 3], strand: ["+", "-"] },
    } as never);
    const s = b.view({ ref: "v", label: "V" }).block({ sequence: g });
    b.preset.sashimi(s, j, { range: [1, 6] });
    b.preset.sashimi(s, j, { label: "Unstranded", strand: null, count: "reads" });
    const doc = b.getDocument();
    const tracks = doc.root.children?.[1].children?.[0].children?.filter((c) => c.kind === "track") ?? [];
    const width = (field: string, range: number[]) => ({ theme: "field_scale", params: { field, range } });
    const common = { type: "pairs", data: "j", display: "arcs", labels: true };
    expect(tracks.map((t) => [t.params, t.children?.map((l) => l.params)])).toEqual([
      [
        { label: "Junctions (+)", height: 3 },
        [{ ...common, side: "above", width: width("count", [1, 6]), filter: { field: "strand", values: ["+"] } }],
      ],
      [
        { label: "Junctions (-)", height: 3 },
        [{ ...common, side: "below", width: width("count", [1, 6]), filter: { field: "strand", values: ["-"] } }],
      ],
      [{ label: "Unstranded", height: 3 }, [{ ...common, side: "above", width: width("reads", [1, 8]) }]],
    ]);
    const reg = createRegistries();
    register(reg, [KernelRegistry]);
    expect(validate(doc, reg)).toEqual([]);
  });
});

describe("sources written in place (CB-08)", () => {
  it("equal params in another key order are one source; inline values by the value itself", () => {
    const b = createSeqViewSpec();
    b.data.sequence({ ref: "a", source: { uri: "https://example.org/a.fa", format: "fasta" } });
    b.data.sequence({ ref: "b", source: { format: "fasta", uri: "https://example.org/a.fa" } });
    b.data.sequence({ ref: "c", source: { value: ">c\nMK", format: "fasta" } });
    b.data.sequence({ ref: "d", source: { format: "fasta", value: ">c\nMK" } });
    b.data.sequence({ ref: "e", source: { value: ">e\nMK", format: "fasta" } });
    const sources = (b.getDocument().root.children?.[0]?.children ?? []).filter((n) => n.kind === "source");
    expect(sources.map((s) => s.ref)).toEqual(["a.src", "c.src", "e.src"]);
  });

  it("nested params in another key order, and a field left undefined, are one source too (LA-41)", () => {
    const b = createSeqViewSpec();
    const uri = "https://example.org/a.tsv";
    b.data.sequence({ ref: "a", source: { uri, format: "tsv", declared: { x: 1, y: [2, { p: 1, q: 2 }] } } } as never);
    b.data.sequence({ ref: "b", source: { declared: { y: [2, { q: 2, p: 1 }], x: 1 }, format: "tsv", uri } } as never);
    b.data.sequence({
      ref: "c",
      source: { uri, format: "tsv", declared: { x: 1, y: [2, { p: 1, q: 2 }] }, max_bytes: undefined },
    } as never);
    const sources = (b.getDocument().root.children?.[0]?.children ?? []).filter((n) => n.kind === "source");
    expect(sources.map((s) => s.ref)).toEqual(["a.src"]);
  });
});
