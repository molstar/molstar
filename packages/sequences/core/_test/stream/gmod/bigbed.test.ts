/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// BigBed against the recorded slices of UCSC's gencodeV50.bb (bigGenePred: the knownGene track). The expected genes are
// the REST API's (getData/track, knownGene) over the same regions, recorded in fixtures/expected.json.

import { readFileSync } from "node:fs";
import { StreamNotServed } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { CASES, makers, req } from "./cases.js";
import { sliceReader } from "./slices.js";

interface Gene {
  name: string;
  start: number;
  end: number;
  strand: string;
  exons: [number, number][];
}
const expected = JSON.parse(readFileSync(new URL("./fixtures/expected.json", import.meta.url), "utf8"))
  .bigbed as Record<string, Gene[]>;
const signal = new AbortController().signal;
const provider = (maxItems?: number) => makers.bigbed(sliceReader("bigbed"), maxItems);

describe("BigBed", () => {
  it("lists the chromosomes; a level of counts, then the file's zoom levels", async () => {
    const p = await provider().fetch(req(null), signal);
    if (p.kind !== "manifest") throw new Error("not a manifest");
    expect(p.manifest.systems.find((s) => s.key === "chr21")?.extent).toEqual([0, 46_709_983]);
    // One level counted from items, then the file's zoom levels (ST4b).
    expect(p.manifest.levels.slice(0, 3)).toEqual([256, 44_748, 178_992]);
    expect(p.manifest.levels.length).toBe(9);
    expect(p.manifest.vocabulary).toEqual({});
  });

  for (const name of ["chr1 genes", "chr21 genes"])
    it(`${name}: the genes, their exons as rows of one id, are the API's`, async () => {
      const r = CASES.bigbed.find(([n]) => n === name)?.[1] ?? req(null);
      const p = await provider().fetch(r, signal);
      if (p.kind !== "items") throw new Error("not items");
      const t = p.items as Record<string, (string | number)[]>;
      const got = new Map<string, Gene>();
      t.id.forEach((id, j) => {
        const g = got.get(id as string) ?? {
          name: t.name[j] as string,
          start: 0,
          end: 0,
          strand: t.strand[j] as string,
          exons: [],
        };
        g.exons.push([t.start[j] as number, t.stop[j] as number]);
        got.set(id as string, g);
      });
      for (const [id, g] of got) [g.start, g.end] = [Number(id.split("@")[1].split("-")[0]), Number(id.split("-")[1])];
      const sorted = (gs: Gene[]) =>
        gs
          .map((g) => ({ ...g, start: +g.start, end: +g.end }))
          .sort((a, b) => a.name.localeCompare(b.name) || a.start - b.start);
      expect(sorted([...got.values()])).toEqual(sorted(expected[name]));
      expect(t.label).toEqual(t.name);
      expect(Object.keys(t)).toEqual(
        expect.arrayContaining(["start", "stop", "id", "label", "name", "strand", "geneName", "geneType"]),
      );
      expect(Object.keys(t)).not.toContain("blockSizes");
    });

  it("answers a level with the items overlapping each bin", async () => {
    const [, r] = CASES.bigbed[3];
    const p = await provider().fetch(r, signal);
    if (p.kind !== "bins") throw new Error("not bins");
    const count = p.bins.count as Uint32Array;
    expect(p.level).toBe(0);
    expect(count.length).toBe(256);
    const want = Array.from({ length: 256 }, (_, j) => {
      const [a, b] = [r.start + 256 * j, r.start + 256 * (j + 1)];
      return expected["chr21 level 0"].filter((g) => g.start < b && g.end > a).length;
    });
    expect([...count]).toEqual(want);
    expect(want.some((c) => c > 0)).toBe(true);
  });

  it("answers a data request holding more than maxItems with the counts of the finest level that fits (the zoom gate)", async () => {
    const [, r] = CASES.bigbed[1];
    const p = await provider(100).fetch(r, signal); // 334 genes
    if (p.kind !== "bins") throw new Error("not bins");
    expect([p.level, p.start, p.stop, (p.bins.count as Uint32Array).length]).toEqual([0, 10_000, 75_536, 256]);
    const whole = await provider(1000).fetch(r, signal);
    expect(whole.kind).toBe("items");
  });

  it("fails clearly for what it does not serve", async () => {
    const p = provider();
    await expect(p.fetch(req("chr21", 0, 256, 9), signal)).rejects.toBeInstanceOf(StreamNotServed);
    await expect(p.fetch(req("chrZ", 0, 256), signal)).rejects.toThrow("BigBed: no chromosome chrZ");
  });

  it("fails naming the range when a region is not in the recorded slices", async () => {
    await expect(provider().fetch(req("chr2", 50_000_000, 50_065_536), signal)).rejects.toThrow(
      /bigbed: bytes \d+\.\.\d+ are not in the recorded slices/,
    );
  });
});

describe("BigBed under a host's naming (ST4b)", () => {
  it("labels items by a named field and takes another as their category, its vocabulary the manifest's", async () => {
    const { bigBedProvider } = await import("../../../src/stream/gmod/index.js");
    const types = ["coding", "nonCoding", "pseudo"];
    const p = bigBedProvider(sliceReader("bigbed"), {
      label: "geneName",
      category: "transcriptClass",
      vocabulary: types,
    });
    const m = await p.fetch(req(null), signal);
    if (m.kind !== "manifest") throw new Error("not a manifest");
    expect(m.manifest.vocabulary).toEqual({ category: types });
    const [, r] = CASES.bigbed[1];
    const got = await p.fetch(r, signal);
    if (got.kind !== "items") throw new Error("not items");
    const t = got.items as Record<string, string[]>;
    expect(t.label).toEqual(t.geneName);
    expect(t.category).toEqual(t.transcriptClass);
    expect(new Set(t.category).size).toBeGreaterThan(1);
  });
});
