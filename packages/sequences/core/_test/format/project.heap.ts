/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for the projection core (VW5; seqspec:interaction#channels: after a route's first use nothing allocates),
// on RB0b's stepped `heapScaling` (a step is one hover change): `project` through a codon route split by an intron, its
// inverse image, links across sections, a source used twice and a row through its binding; `close` by the hull search
// (a nucleotide) and by containment (a residue), then the closed extents projected.

import { Extents } from "@molstar/seq-core/model";
import { beforeAll, describe, expect, it } from "vitest";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from "../model/heap.js";
import { constructs, gene, open, rows } from "./project-fixtures.js";

/** Steps in the shorter run; the longer runs four times as many. */
const N = 50_000;
let g: Awaited<ReturnType<typeof open>>;
let c: Awaited<ReturnType<typeof open>>;
let r: Awaited<ReturnType<typeof open>>;
const sink = {
  n: 0,
  add(s: number, t: number) {
    this.n += t - s;
  },
};

beforeAll(async () => {
  [g, c, r] = await Promise.all([open(gene({ chromosome: true })), open(constructs()), open(rows())]);
}, 60_000);

describe("the projection allocates nothing per move", () => {
  it("project: residues onto their codons, nucleotides onto residues, across the link to the chromosome", () => {
    const ext = Extents.create(4);
    const [protein, unit] = [g.system("protein"), g.system("gene/unit")];
    const blocks = [g.block("gene"), g.block("protein"), g.block("chr")];
    const step = (i: number) => {
      ext.count = 0;
      Extents.add(ext, protein, i % 20, (i % 20) + 1);
      Extents.add(ext, unit, i % 120, (i % 120) + 1);
      for (let b = 0; b < blocks.length; b++) g.p.project(ext, blocks[b], sink);
    };
    const result = heapScaling(N, step);
    expect(sink.n).toBeGreaterThan(0);
    expect(isFlat(result), describeScaling(result, N)).toBe(true);
    // The template's control (CF-05): the same step with a leak is caught.
    expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
  });

  it("project: a source used twice, and a row through its binding to its chain", () => {
    const ext = Extents.create(4);
    const [ub, diub] = [c.system("ub"), c.block("diub")];
    const [row, x, a] = [r.system("a/0"), r.block("x"), r.block("a")];
    const rowExt = Extents.create(4);
    const result = heapScaling(N, (i) => {
      ext.count = 0;
      Extents.add(ext, ub, i % 76, (i % 76) + 1);
      c.p.project(ext, diub, sink);
      rowExt.count = 0;
      Extents.add(rowExt, row, i % 4, (i % 4) + 1);
      r.p.project(rowExt, x, sink);
      r.p.project(rowExt, a, sink);
    });
    expect(isFlat(result), describeScaling(result, N)).toBe(true);
  });

  it("close: a nucleotide by the hull search, a residue by containment, then projected", () => {
    const ext = Extents.create(8);
    const [chr, protein, gene, block] = [
      g.system("chr/unit"),
      g.system("protein"),
      g.system("gene/unit"),
      g.block("chr"),
    ];
    let added = 0;
    const result = heapScaling(N, (i) => {
      ext.count = 0;
      Extents.add(ext, i % 2 ? chr : gene, 100 * (i % 2) + (i % 120), 100 * (i % 2) + (i % 120) + 1);
      Extents.add(ext, protein, i % 20, (i % 20) + 1);
      g.p.close(ext);
      added += ext.count - 2;
      g.p.project(ext, block, sink);
    });
    expect(added).toBeGreaterThan(0);
    expect(isFlat(result), describeScaling(result, N)).toBe(true);
  });
});
