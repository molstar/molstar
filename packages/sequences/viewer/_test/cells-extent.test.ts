/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// GM-02 and GM-03: a draw's edge-rate work is bounded by what its block shows, not its system's length
// (seqspec:viewer#performance): case 20's protein letters placed on chr2 through its codons build their classes over
// the protein's parts, never the 242 Mb domain (7.7 s and 242 MB on the first frame before); a ruler numbering the
// genome by that protein builds its tick lists over the protein's reach (about 16 s before), and two pieces of one
// unit each over its own.

import { createSeqViewSpec, translate } from "@molstar/seq-core/builder";
import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { Coords, Numbering } from "@molstar/seq-core/model";
import { describe, expect, it, vi } from "vitest";
import { CASES } from "../../core/_test/format/cases.js";
import { CODING_CDS } from "../../core/_test/format/stream-provider.js";
import { Dirty } from "../src/draw/composer.js";
import { composed, recording, textsOf } from "./composer-fixtures.js";
import { lit } from "./overlay-fixtures.js";

/** The length of every class pyramid built (`cellClasses`), a count rather than a clock (LC-34). */
const built = vi.hoisted(() => [] as number[]);
vi.mock("../src/reductions.js", async (original) => {
  const m = await original<typeof import("../src/reductions.js")>();
  return {
    ...m,
    cellClasses: (...a: Parameters<typeof m.cellClasses>) => {
      built.push(a[3]);
      return m.cellClasses(...a);
    },
  };
});

describe("GM-02: cell classes over what the block shows", () => {
  it("case 20's genome section builds its classes over the protein's parts, its buffers small", async () => {
    const t = await composed(CASES.case20(), { section: "coding.genome" });
    const m0 = process.memoryUsage().arrayBuffers;
    built.length = 0;
    t.frame(recording(), Dirty.all);
    const mb = (process.memoryUsage().arrayBuffers - m0) / 1e6;
    expect(mb).toBeLessThan(32);
    // The protein's codons span about 1.2 Mb of chr2's 242: classes over 2^21 positions at most.
    expect(built.length).toBeGreaterThan(0);
    expect(Math.max(...built)).toBeLessThanOrEqual(2 ** 21);
    t.ctx.dispose();
  }, 60_000);
});

const CODING = "ATGGCTTCCAAGGAACTGACCCCACGTGACGGCCATTGGTACTTCAACAGCCAGGTGGAG";
/** A gene of `n` positions whose 161 codons start at `at` (not aligned to any power of two). */
function codingDoc(n: number, at: number) {
  const cds = `${CODING.repeat(8)}TAA`;
  const value = "ACGT".repeat(n / 4).split("");
  value.splice(at, cds.length, ...cds);
  const b = createSeqViewSpec();
  const gene = b.data.sequence({ ref: "gene", alphabet: "dna", value: value.join("") });
  const protein = translate(b.data, { ref: "protein", from: gene, cds: [[at, at + cds.length]] });
  b.view({ ref: "v", label: "V" }).section({ ref: "v.main" }).block({ sequence: gene }).track({ label: "" }).letters({
    of: protein,
  });
  return { doc: b.getDocument(), stop: at + cds.length };
}

describe("GM-02: the classes from an aligned offset", () => {
  for (const cell of [0.006, 0.02, 0.1])
    it(`zoomed out (${cell} px a position): the protein's columns colored, and only those`, async () => {
      const at = 70_003;
      const { doc, stop } = codingDoc(200_000, at);
      const t = await lit(doc, { width: 1200, height: 300 });
      const x0 = Math.floor(at * cell - 600);
      t.vp?.setView(cell, Math.max(0, x0), "api");
      const be = t.frame(recording());
      const ax = t.ctx.layout.axis("v.main");
      if (!ax) throw new Error("no axis");
      const [left, right] = [ax.x(0, at as never), ax.x(0, stop as never)];
      const r = be.rectItems;
      const xs: number[][] = [];
      for (let i = 0; i < r.n; i++) if (r.h[i] > 4) xs.push([r.x[i], r.x[i] + r.w[i]]);
      expect(xs.length).toBeGreaterThan(0);
      for (const [a, e] of xs) {
        expect(a).toBeGreaterThanOrEqual(Math.floor(left) - 1);
        expect(e).toBeLessThanOrEqual(Math.ceil(right) + 1);
      }
      t.ctx.dispose();
    });
});

describe("GM-02: an extent across a high power of two (LC-25)", () => {
  // Under three times the extent wherever it lies: an offset aligned to a power of two past the whole range fell to 0.
  for (const [what, n, at] of [
    ["a protein across 2^20", 1_100_000, 2 ** 20 - 200],
    ["a protein across 2^18", 400_000, 2 ** 18 - 240],
  ] as const)
    it(`${what}: classes under three times the protein's extent, its columns colored and only those`, async () => {
      const { doc, stop } = codingDoc(n, at);
      const t = await lit(doc, { width: 1200, height: 300 });
      const cell = 0.006;
      built.length = 0;
      t.vp?.setView(cell, Math.max(0, Math.floor(at * cell - 600)), "api");
      const be = t.frame(recording());
      expect(built.length).toBeGreaterThan(0);
      for (const length of built) expect(length).toBeLessThan(3 * (stop - at));
      const ax = t.ctx.layout.axis("v.main");
      if (!ax) throw new Error("no axis");
      const [left, right] = [ax.x(0, at as never), ax.x(0, stop as never)];
      const r = be.rectItems;
      let cells = 0;
      for (let i = 0; i < r.n; i++)
        if (r.h[i] > 4) {
          cells++;
          expect(r.x[i]).toBeGreaterThanOrEqual(Math.floor(left) - 1);
          expect(r.x[i] + r.w[i]).toBeLessThanOrEqual(Math.ceil(right) + 1);
        }
      expect(cells).toBeGreaterThan(0);
      t.ctx.dispose();
    });
});

type Node = { kind: string; ref?: string; params?: Record<string, unknown>; children?: Node[] };
/** Case 20 with its genome section's ruler numbering the genome by protein P1 (`of`). */
function rulerOf(): SeqViewSpecDoc {
  const d = CASES.case20() as unknown as { root: Node };
  const walk = (n: Node) => {
    if (n.kind === "section" && n.ref === "coding.genome") {
      const ruler = n.children?.find((c) => c.kind === "track") as Node;
      const layer = ruler.children?.find((c) => c.kind === "layer") as Node;
      expect(layer.params?.type ?? "ruler").toBe("ruler");
      layer.params = { ...layer.params, of: "p1" };
    }
    for (const c of n.children ?? []) walk(c);
  };
  walk(d.root);
  return d as unknown as SeqViewSpecDoc;
}

describe("GM-03: the ruler's tick lists over what it numbers", () => {
  it("a ruler `of` the protein on chr2 draws its first frame in well under a second, labeling its residues", async () => {
    const t = await composed(rulerOf(), { section: "coding.genome" });
    const be = recording();
    // Positions mapped to the protein (the tick lists' and the labels'), a count rather than a clock (LC-34): the
    // protein's reach, not chr2's 242 M.
    const maps = [vi.spyOn(Coords, "toAt"), vi.spyOn(Coords, "forEachFromAt")];
    t.frame(be, Dirty.all);
    const n = maps[0].mock.calls.length + maps[1].mock.calls.length;
    for (const m of maps) m.mockRestore();
    expect(n).toBeGreaterThan(0);
    expect(n).toBeLessThan(2 ** 23);
    // P1's residue numbers (prefixed with its label), at the focus's zoom.
    // Zoomed to the first exon's start: P1's residue numbers (residue 1 at the CDS's first codon).
    t.vp?.zoomTo(0, CODING_CDS[0][0] - 20, CODING_CDS[0][0] + 280);
    const near = recording();
    t.frame(near, Dirty.all);
    expect(
      textsOf(near)
        .filter((s) => /^\d+$/.test(s))
        .map(Number),
    ).toEqual([25, 50, 75]);
    t.ctx.dispose();
  }, 60_000);
});

describe("GM-03: two pieces of one unit far apart (LC-33)", () => {
  it("each piece's ticks come from its own positions, not the hull across the gap", async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", alphabet: "dna", value: "ACGT".repeat(500_000) });
    // A numbering with no 0 (read from the tick lists, not counted by one), the rulers' default.
    a.numbering({ name: "auth", method: "offset", first: -100, skip_zero: true, default: true });
    const sec = b.view({ ref: "v", label: "V" }).section({ ref: "v.main" });
    sec.block({ sequence: a, range: [1_000_000, 1_000_300] });
    sec.block({ sequence: a, range: [1_999_000, 1_999_300] });
    sec.track({ label: "" }).ruler();
    const t = await lit(b.getDocument(), { width: 1200, height: 300 });
    const labels = vi.spyOn(Numbering, "labelOf");
    const be = t.frame(recording());
    expect(textsOf(be).length).toBeGreaterThan(0);
    // About 600 positions' labels (twice each): the hull across the gap read a million.
    expect(labels.mock.calls.length).toBeGreaterThan(0);
    expect(labels.mock.calls.length).toBeLessThan(10_000);
    labels.mockRestore();
    t.ctx.dispose();
  });
});
