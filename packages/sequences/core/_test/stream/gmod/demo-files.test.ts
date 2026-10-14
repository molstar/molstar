/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The synthetic genome files (demo-files.ts) read by the four readers as the File formats story reads them: each
// manifest, a region's data, and what the text they were written from says is there.

import { describe, expect, it } from "vitest";
import { bigBedProvider, bigWigProvider, twoBitProvider, vcfProvider } from "../../../src/stream/gmod/index.ts";
import { req } from "./cases.ts";
import { DEMO_SYSTEMS, demoGenome, memoryReader } from "./demo-files.ts";

const signal = new AbortController().signal;
const g = demoGenome();
const systems = DEMO_SYSTEMS.map(([key, n]) => ({ key, extent: [0, n] }));
/** chrS1's letters as the FASTA text has them. */
const chrS1 = g.text.fasta.split(">")[1].split("\n").slice(1).join("");

describe("the synthetic genome files", () => {
  it("2bit: the sequences, an N run and soft-masked repeats", async () => {
    const p = twoBitProvider(memoryReader(g.twobit));
    const m = await p.fetch(req(null), signal);
    expect(m.kind === "manifest" && m.manifest.systems.map((s) => ({ key: s.key, extent: s.extent }))).toEqual(systems);
    const r = await p.fetch(req("chrS1", 390, 3_010), signal);
    expect(r.kind === "letters" && r.letters).toBe(chrS1.slice(390, 3_010));
    expect(chrS1.slice(395, 405)).toMatch(/^N{5}[ACGT]{5}$/);
    expect(chrS1.slice(3_000, 3_010)).toMatch(/^[acgt]{10}$/);
  });

  it("BigWig: GC percent in 50-base windows, its range in the manifest", async () => {
    const p = bigWigProvider(memoryReader(g.bigwig));
    const m = await p.fetch(req(null), signal);
    if (m.kind !== "manifest") throw new Error(m.kind);
    expect(m.manifest.systems.map((s) => s.key)).toEqual(["chrS1", "chrS2"]);
    expect(m.manifest.levels).toEqual([]);
    const [lo, hi] = m.manifest.domain as [number, number];
    expect(lo).toBeGreaterThanOrEqual(0);
    expect(hi).toBeLessThanOrEqual(100);
    const r = await p.fetch(req("chrS1", 6_000, 6_100), signal);
    if (r.kind !== "items") throw new Error(r.kind);
    const value = r.items.value as Float64Array;
    const gc = (s: string) =>
      Math.round((100 * [...s.toUpperCase()].filter((c) => c === "G" || c === "C").length) / 50);
    expect([value[0], value[50]]).toEqual([gc(chrS1.slice(6_000, 6_050)), gc(chrS1.slice(6_050, 6_100))]);
  });

  it("BigBed: a row per exon, a gene's exons under one id", async () => {
    const p = bigBedProvider(memoryReader(g.bigbed), { category: "strand", vocabulary: ["+", "-"] });
    const m = await p.fetch(req(null), signal);
    expect(m.kind === "manifest" && m.manifest.vocabulary).toEqual({ category: ["+", "-"] });
    const r = await p.fetch(req("chrS1", 0, 10_000), signal);
    if (r.kind !== "items") throw new Error(r.kind);
    const t = r.items as Record<string, unknown[]>;
    expect(t.label.slice(0, 5)).toEqual(["GENA", "GENA", "GENA", "GENA", "GENB"]);
    expect([t.start.slice(0, 4), t.stop.slice(0, 4)]).toEqual([
      [2_000, 3_800, 5_200, 7_000],
      [2_300, 4_100, 5_600, 7_500],
    ]);
    expect(new Set(t.id.slice(0, 4)).size).toBe(1);
    expect(t.category[4]).toBe("-");
    const s2 = await p.fetch(req("chrS2", 0, 6_000), signal);
    expect(s2.kind === "items" && (s2.items as Record<string, unknown[]>).label).toEqual(["GENE", "GENE"]);
  });

  it("tabix VCF: the header's lengths, a region's records by the index, read in BGZF blocks", async () => {
    const read = memoryReader(g.vcf);
    const p = vcfProvider(read, g.tbi, { info: ["CLNSIG"] });
    const m = await p.fetch(req(null), signal);
    expect(m.kind === "manifest" && m.manifest.systems.map((s) => ({ key: s.key, extent: s.extent }))).toEqual(systems);
    const lines = g.text.vcf.split("\n").filter((l) => l.startsWith("chrS1\t"));
    const inRegion = lines.filter((l) => {
      const pos = Number(l.split("\t")[1]);
      return pos > 10_000 && pos <= 14_000;
    });
    const r = await p.fetch(req("chrS1", 10_000, 14_000), signal);
    if (r.kind !== "items") throw new Error(r.kind);
    const t = r.items as Record<string, unknown[]>;
    expect(t.name).toEqual(inRegion.map((l) => l.split("\t")[2]));
    expect(t.CLNSIG).toEqual(inRegion.map((l) => l.split("\t")[7].slice("CLNSIG=".length)));
    expect(new Set(t.category)).toEqual(
      new Set(["SNV", "insertion", "deletion"].filter((c) => t.category.includes(c))),
    );
    expect(t.category).toContain("SNV");
    // The open read the header's blocks; the region one run of blocks (and its last block's header first).
    expect(read.reads.length).toBeGreaterThanOrEqual(2);
  });
});
