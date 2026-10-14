/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// VL-08 (seqspec:format#layer-types `secondary`): a strand's arrowhead points to its C-terminal end on screen (a
// reversed piece's left), drawn on the piece holding that end only; a strand split by a filter's break draws its other
// piece as a body.

import { createSeqViewSpec, translate } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { drawn, record } from "./draw-fixtures.js";

/** A strand [10, 20) on a 40-residue block in `direction`; `filtered`: positions 14 and 15 hidden. */
async function strand(direction: 1 | -1, filtered = false) {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: "a", value: "MKVLATGHEK".repeat(4) });
  const ss = b.data.annotation({
    ref: "ss",
    geometry: "ranges",
    target: a,
    value: { start: [10], stop: [20], category: ["strand"] },
  });
  const keep = b.data.annotation({
    ref: "keep",
    geometry: "ranges",
    target: a,
    value: { start: [0, 16], stop: [14, 40] },
  });
  const v = b.view({ ref: "v", label: "V", height: "box" });
  const sec = v.section({ ref: "v.s", direction, ...(filtered ? { positions: keep } : {}) });
  sec.block({ sequence: a });
  sec.track({ label: "S" }).secondary({ data: ss });
  const t = await drawn(b.getDocument(), { width: 600, height: 200 });
  const { rec, frame } = record(t.ctx, "S");
  const p = rec.pathItems;
  const shapes: number[][] = [];
  const tips: number[] = [];
  for (let i = 0; i < p.n; i++) {
    const c1 = i + 1 < p.n ? p.cstart[i + 1] : p.nc;
    const xs: number[] = [];
    const ys: number[] = [];
    for (let c = p.cstart[i]; c < c1; c += 2) {
      xs.push(p.coords[c]);
      ys.push(p.coords[c + 1]);
    }
    shapes.push(xs);
    // The tip: the vertex on the coil line, half way between the head's top and bottom.
    const mid = (Math.min(...ys) + Math.max(...ys)) / 2;
    tips.push(xs[ys.indexOf(mid)] ?? Number.NaN);
  }
  const v0 = frame.v;
  const q = (i: number) => v0.shownOf(0, i);
  return { shapes, tips, left: (i: number) => v0.xd(0, q(i)), right: (i: number) => v0.xd(0, q(i) + 1) };
}

describe("VL-08: a strand's head at its C-terminal end, on screen", () => {
  it("forward: one arrow, its tip on residue 19's right edge", async () => {
    const { shapes, tips, left, right } = await strand(1);
    expect(shapes.length).toBe(1);
    expect(tips).toEqual([right(19)]);
    expect(Math.max(...shapes[0])).toBe(right(19));
    expect(Math.min(...shapes[0])).toBe(left(10));
    // Seven vertices: the body and the head.
    expect(shapes[0].length).toBe(7);
  });

  it("reversed: one arrow, its tip on residue 19's left edge (residue 19 is on the left)", async () => {
    const { shapes, tips, left, right } = await strand(-1);
    expect(left(19)).toBeLessThan(left(10));
    expect(shapes.length).toBe(1);
    expect(tips).toEqual([left(19)]);
    // The tip is the vertex alone at the far left; the body ends at residue 10's right edge.
    expect(Math.min(...shapes[0])).toBe(left(19));
    expect(Math.max(...shapes[0])).toBe(right(10));
    expect(shapes[0].length).toBe(7);
  });

  it("across a filter's break: one piece (the shown positions run on), one head", async () => {
    for (const direction of [1, -1] as const) {
      const { shapes, tips, left, right } = await strand(direction, true);
      expect(shapes.length).toBe(1);
      expect(tips).toEqual([direction > 0 ? right(19) : left(19)]);
    }
  });

  /** A strand on a 40-residue block, `cyclic` from `origin`, `keep` its filter: its shapes' vertex counts and tips. */
  async function strandOn(o: { start: number; stop: number; origin?: number; keep?: [number[], number[]] }) {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", value: "MKVLATGHEK".repeat(4), cyclic: o.origin !== undefined });
    const ss = b.data.annotation({
      ref: "ss",
      geometry: "ranges",
      target: a,
      value: { start: [o.start], stop: [o.stop], category: ["strand"] },
    });
    const keep = o.keep
      ? b.data.annotation({ ref: "keep", geometry: "ranges", target: a, value: { start: o.keep[0], stop: o.keep[1] } })
      : null;
    const sec = b
      .view({ ref: "v", label: "V", height: "box" })
      .section({ ref: "v.s", ...(keep ? { positions: keep } : {}) });
    sec.block({ sequence: a, ...(o.origin !== undefined ? { origin: o.origin } : {}) });
    sec.track({ label: "S" }).secondary({ data: ss });
    const t = await drawn(b.getDocument(), { width: 600, height: 200 });
    const { rec, frame } = record(t.ctx, "S");
    const p = rec.pathItems;
    const sizes: number[] = [];
    const ends: number[] = [];
    for (let i = 0; i < p.n; i++) {
      const c1 = i + 1 < p.n ? p.cstart[i + 1] : p.nc;
      sizes.push((c1 - p.cstart[i]) / 2);
      let hi = Number.NEGATIVE_INFINITY;
      for (let c = p.cstart[i]; c < c1; c += 2) hi = Math.max(hi, p.coords[c]);
      ends.push(hi);
    }
    const v0 = frame.v;
    return { sizes, ends, right: (i: number) => v0.xd(0, v0.shownOf(0, i) + 1) };
  }

  it("across a cyclic origin (38 to 2): one head, on the piece holding residue 1, from either origin (LC-04)", async () => {
    for (const origin of [0, 30]) {
      const { sizes, ends, right } = await strandOn({ start: 38, stop: 2, origin });
      expect(sizes.filter((n) => n === 7).length, `origin ${origin}`).toBe(1);
      expect(ends[sizes.indexOf(7)], `origin ${origin}`).toBe(right(1));
    }
  });

  it("its C-terminal residue filtered out: the head on the last residue shown (LC-05)", async () => {
    const { sizes, ends, right } = await strandOn({
      start: 10,
      stop: 20,
      keep: [
        [0, 20],
        [19, 40],
      ],
    });
    expect(sizes).toEqual([7]);
    expect(ends[0]).toBe(right(18));
  });

  it("split by an intron (a protein's strand placed on its gene): the head on the second exon's piece only", async () => {
    const b = createSeqViewSpec();
    const GENE = `GAATTCGCCACCATGGCTTCCAAGGAACTGACCCCACGTGACGGCCGTAAGTCTTTCCTAACTTCTGTTTCCACAGATTGGTACTTCAACAGCCAGGTGGAGTAAGCGG`;
    const g = b.data.sequence({ ref: "gene", alphabet: "dna", value: GENE });
    const protein = translate(b.data, {
      ref: "protein",
      from: g,
      cds: [
        [12, 46],
        [76, 105],
      ],
    });
    // Residues [8, 14): codons across the intron (exon 1 ends inside residue 11).
    const ss = b.data.annotation({
      ref: "ss",
      geometry: "ranges",
      target: protein,
      value: { start: [8], stop: [14], category: ["strand"] },
    });
    b.view({ ref: "v", label: "V" }).block({ sequence: g }).track({ label: "S" }).secondary({ data: ss });
    const t = await drawn(b.getDocument(), { width: 1200, height: 200 });
    t.ctx.viewports.get("v.main")?.setView(10, 0, "api");
    const { rec } = record(t.ctx, "S");
    const p = rec.pathItems;
    const sizes = Array.from({ length: p.n }, (_, i) => ((i + 1 < p.n ? p.cstart[i + 1] : p.nc) - p.cstart[i]) / 2);
    // Two pieces (one per exon): a body, then the arrow.
    expect(sizes).toEqual([4, 7]);
  });
});
