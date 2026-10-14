/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The builder has no runtime dependency (seqspec:format#builders), so it keeps small copies of what it needs; each is
// named as the twin of its original, and this file keeps every pair equal: the format and version, the namespace
// syntax, the float32 printer, the CIGAR and WUSS readers, NCBI's genetic codes and codon reading, the codon mapping's
// order rule and alphabet, and the kernel layer types.
// The helpers are not the builder's public API, so they are read from their file.

import { cigarFromRows } from "@molstar/seq-core/algorithms";
import { createSeqViewSpec, translate } from "@molstar/seq-core/builder";
import { float32Number, GC_PRT, GeneticCode, Molecules, SeqSets, Stockholm } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import {
  cigarFromRows as builderCigar,
  SvsFormat as builderFormat,
  isNamespaced as builderIsNamespaced,
  SvsVersion as builderVersion,
  wussPairs as builderWuss,
  codonOf,
  f32,
  GENETIC_CODES,
} from "../../src/builder/helpers.js";
import { isNamespaced, KERNEL_LAYER_TYPE_NAMES, SvsFormat, SvsVersion } from "../../src/format/index.js";
import { KernelLayerTypes } from "../../src/format/tables/kernel.js";

/** A small deterministic generator (an LCG over 32 bits). */
function random(seed: number) {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 2 ** 32;
  };
}

describe("the builder's twins", () => {
  it("the format and version of a document", () => {
    expect([builderFormat, builderVersion]).toEqual([SvsFormat, SvsVersion]);
    expect(createSeqViewSpec().getDocument()).toMatchObject({ format: SvsFormat, version: SvsVersion });
  });

  it("isNamespaced: the same names, short names, reverse-DNS namespaces, case and odd characters", () => {
    const names = [
      "uniprot:Domain",
      "org.example:dots",
      "org.example.sub-1:Name_2.x",
      "svs:thing",
      "ns:1abc",
      "Org:x",
      "org.Example:x",
      "org..x:y",
      ".org:x",
      "org.:x",
      "org:",
      ":x",
      "org:-x",
      "org:x y",
      "org:x:y",
      "1org:x",
      "plain",
      "",
      "a:b",
      "a_b-c:D-e",
    ];
    for (const n of [...names, 7, null, undefined, {}]) expect(builderIsNamespaced(n)).toBe(isNamespaced(n));
  });

  it("f32: the shortest decimals of a float32, on every kind of float32", () => {
    const next = random(7);
    const bits = new Uint32Array(1);
    const f = new Float32Array(bits.buffer);
    const xs = [
      0,
      -0,
      0.15,
      0.3,
      0.1,
      1,
      16777216,
      3.4028234663852886e38,
      1.401298464324817e-45,
      2 ** -96,
      2 ** 87,
      2 ** 90,
      1 / 3,
    ];
    for (let i = 0; i < 40_000; i++) {
      bits[0] = Math.floor(next() * 2 ** 32);
      if (Number.isFinite(f[0])) xs.push(f[0]);
    }
    for (let e = -149; e < 128; e++) xs.push(Math.fround(2 ** e), Math.fround(-(2 ** e) * 1.5));
    xs.push(Number.NaN, Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY);
    expect(xs.filter((x) => !Object.is(f32(x), float32Number(x)))).toEqual([]);
  });

  it("cigarFromRows: runs, gaps in both rows skipped, unequal rows", () => {
    const next = random(11);
    const row = (n: number) => Array.from({ length: n }, () => "ACG-"[Math.floor(next() * 4)]).join("");
    const pairs = [
      ["", ""],
      ["-", "-"],
      ["MKV", "MKV"],
      ["-MKVL", "AMK-L"],
      ["--A-", "B--C"],
    ];
    for (let i = 0; i < 300; i++) {
      const n = Math.floor(next() * 30);
      pairs.push([row(n), row(n)]);
    }
    for (const [a, b] of pairs) expect(builderCigar(a, b)).toBe(cigarFromRows(a, b));
    expect(() => builderCigar("A", "AA")).toThrow(RangeError);
    expect(() => cigarFromRows("A", "AA")).toThrow(RangeError);
  });

  it("wussPairs: the same partners, and unbalanced input is an error in both", () => {
    const next = random(13);
    const alphabet = "<>()[]{}Aa.:_,-~";
    const strings = ["<<A..a>>", "((.))", "..", "", "AAaa", "<(>)", "{[(<>)]}", "Aa<>aA", "(Aa)"];
    for (let i = 0; i < 400; i++) {
      strings.push(
        Array.from({ length: Math.floor(next() * 24) }, () => alphabet[Math.floor(next() * alphabet.length)]).join(""),
      );
    }
    for (const ss of strings) {
      const issues: Parameters<typeof Stockholm.pairs>[1] = [];
      const model = Stockholm.pairs(ss, issues, "");
      if (model) expect(Array.from(builderWuss(ss))).toEqual(Array.from(model));
      else {
        expect(issues.length).toBe(1);
        expect(() => builderWuss(ss)).toThrow(issues[0].message.replace(/^/, "basePairs: "));
      }
    }
  });

  it("GENETIC_CODES are the model's, and codonOf reads every codon of every table as the model does", () => {
    expect(GENETIC_CODES).toEqual(GC_PRT);
    const bases = "TCAGUN";
    const code = (c: string) => c.charCodeAt(0);
    for (const id of Object.keys(GC_PRT).map(Number)) {
      for (const a of bases) {
        for (const b of bases) {
          for (const c of bases) {
            const { residue, start } = codonOf(GENETIC_CODES[id], a + b + c);
            expect(residue).toBe(String.fromCharCode(GeneticCode.translate(id, code(a), code(b), code(c))));
            expect(start).toBe(GeneticCode.isStart(id, code(a), code(b), code(c)));
          }
        }
      }
    }
  });

  it("the layer type names are the kernel's, in the table's order", () => {
    const track = createSeqViewSpec().view({ label: "V" }).track({}) as unknown as Record<string, unknown>;
    expect(Object.keys(track).filter((k) => k !== "layer")).toEqual([...KERNEL_LAYER_TYPE_NAMES]);
    expect(KernelLayerTypes.map((t) => t.name)).toEqual([...KERNEL_LAYER_TYPE_NAMES]);
  });
});

describe("translate and the codon mapping", () => {
  const GENE =
    "GAATTCGCCACCATGGCTTCCAAGGAACTGACCCCACGTGACGGCCGTAAGTCTTTCCTAACTTCTGTTTCCACAGATTGGTACTTCAACAGCCAGGTGGAGTAAGCGGCCGCATAAAGT";
  interface Case {
    dna?: string;
    alphabet?: "dna" | "rna" | "protein";
    cds?: [number, number][];
    start?: number;
    strand?: "+" | "-";
    code?: number;
    empty?: boolean;
  }
  /** Ranges as "12-46 76-105". */
  const r = (ranges: string) => ranges.split(" ").map((x) => x.split("-").map(Number) as [number, number]);
  const CASES: Case[] = [
    { cds: r("12-46 76-105") },
    { cds: r("76-105 12-46") }, // out of order on +
    { cds: r("12-46 40-69") }, // overlap
    { cds: r("12-46 46-49") }, // touching ranges are fine
    { cds: r("12-46 76-104") }, // not a multiple of 3
    { cds: r("12-46 100-121") }, // past the sequence
    { cds: r("12-46 76-105"), strand: "-" }, // ascending on -
    { cds: r("76-105 12-46"), strand: "-" },
    { cds: r("76-105 40-46"), strand: "-" },
    { cds: r("76-105 40-49"), strand: "-" },
    { cds: r("0-6"), code: 7 },
    { cds: r("0-6"), code: 11 },
    { start: 12 },
    { start: 90 },
    { start: 200 },
    { start: -3 },
    { cds: r("12.5-45.5") },
    { cds: r("20-20") },
    { cds: r("0-3"), empty: true, dna: "TAAGGG" },
    { start: 0, empty: true, dna: "TAAGGGCCC" },
    { cds: r("0-6"), dna: "AUGGCC", alphabet: "rna" },
    { cds: r("0-6"), dna: "ATGGCC", alphabet: "rna" },
    { cds: r("0-6"), dna: "AUGGCC", alphabet: "dna" },
    { cds: r("0-6"), dna: "ATGGCC", alphabet: "protein" },
    { cds: r("0-6"), dna: "ATGGCC" },
    { cds: r("0-6"), dna: "AUGGCC" },
    { cds: r("0-6"), dna: "AUGGCT" },
    { cds: r("0-6"), dna: "MKVLAA" },
    { cds: r("0-6"), dna: "MKVLAA", alphabet: "dna" },
    { cds: r("0-12"), dna: "atggccnaaatg" },
    { cds: r("0-12"), dna: "atgrykmswatg" }, // under 90% plain: detected as protein
    { cds: r("0-12"), dna: "ATGRYKMSWATG", alphabet: "dna" },
  ];

  it("throws exactly where the model's codon mapping (or its molecule) reports an issue", () => {
    for (const c of CASES) {
      const dna = c.dna ?? GENE;
      const code = c.code ?? 1;
      const b = createSeqViewSpec();
      const from = b.data.sequence({ ref: "g", value: dna, ...(c.alphabet && { alphabet: c.alphabet }) });
      let protein: string | null = null;
      try {
        const p = translate(b.data, {
          from,
          ...(c.cds && { cds: c.cds }),
          ...(c.start !== undefined && { start: c.start }),
          ...(c.strand && { strand: c.strand }),
          genetic_code: code,
        });
        protein = p.node.params.value as string;
      } catch {}
      const detected = Molecules.detectAlphabet(Uint8Array.from(dna.toUpperCase(), (ch) => ch.charCodeAt(0)));
      const { issues } = SeqSets.build({
        molecules: [
          { key: "g", alphabet: c.alphabet ?? detected, residues: dna, ownUnit: true },
          { key: "p", alphabet: "protein", residues: protein ?? (c.empty ? "" : "MKVLAAGIVGMK"), ownUnit: true },
        ],
        mappings: [
          {
            key: "p.codon",
            from: "p",
            to: "g",
            method: "codon",
            ...(c.cds ? { cds: c.cds } : { start: c.start ?? 0 }),
            ...(c.strand && { strand: c.strand }),
            geneticCode: code,
          },
        ],
      });
      const rejected = issues.some((i) => i.severity === "error");
      expect(rejected, JSON.stringify({ case: c, protein, issues })).toBe(protein === null);
    }
  });
});
