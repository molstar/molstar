/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST7's heap guard: reading a sparse composed path (seqspec:model#paths) allocates nothing per step. A DMD-sized protein
// (79 exons over 2.2 Mb) coded on a lazy 200 Mb chromosome, its path read per element both ways (`Paths.at`,
// `Paths.spanAt`) and by ranges either way (`forEachFrom`: a streamed request's runs; `forEachTo`: an item's image), as
// the stream manager's pass and its draws read it. The positive control is the template's: the same step with a leak.

import { describe, expect, it } from "vitest";
import { type CoordinateSystem, Mappings, type PathStep, Paths, SeqSets } from "../../src/model/index.ts";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from "./heap.ts";
import { gene as dmd } from "./support/generators.ts";

const N = 250_000;
const MASK = 0xfffff;
const AT = 100_000_000;
const g = dmd();
const { set } = SeqSets.build({
  molecules: [
    { key: "chr", alphabet: "dna", lazy: { length: 200_000_000 }, ownUnit: true },
    { key: "p", alphabet: "protein", residues: g.protein, ownUnit: true },
  ],
  mappings: [{ key: "c", from: "p", to: "chr", method: "codon", cds: g.exons.map(([s, t]) => [AT + s, AT + t]) }],
});
const sys = (key: string): CoordinateSystem => {
  const s = { kind: -1, id: -1, row: -1 };
  Mappings.systemOf(set, key, s);
  return s;
};
const steps: PathStep[] = [];
Mappings.path(set, sys("chr/unit"), sys("p/unit"), steps);
const down = Mappings.compose(set, sys("chr/unit"), steps);
Mappings.path(set, sys("p/unit"), sys("chr/unit"), steps);
const up = Mappings.compose(set, sys("p/unit"), steps);
const span = g.exons[g.exons.length - 1][1] - g.exons[0][0];

const ctx = { n: 0 };
let acc = 0;
function runSum(s: number, t: number, c: { n: number }) {
  c.n = (c.n + s + t) & MASK;
}
const read = (i: number) => {
  const r = (i * 7919) % 3685;
  acc = (acc + Paths.at(down, AT + g.exons[0][0] + ((i * 104_729) % span))) & MASK;
  acc = (acc + Paths.spanAt(up, r, i % 3)) & MASK;
  Paths.forEachFrom(down, r, r + (i % 40) + 1, runSum, ctx);
  if ((i & 3) === 0) Paths.forEachTo(down, AT + ((i * 31) % span), AT + ((i * 31) % span) + 5000, runSum, ctx);
};

describe("ST7: a sparse path's reads allocate nothing", () => {
  it("holds the CDS, not the chromosome", () => {
    expect([down.length, down.heldTo.length < 12_000, up.heldSpans?.length]).toEqual([200_000_000, true, 3 * 3685]);
  });

  it("Paths.at, spanAt, forEachFrom and forEachTo over a 79-exon path through 200 Mb", () => {
    const r = heapScaling(N, read);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(Number.isFinite(acc + ctx.n)).toBe(true);
  });

  // CM-05: the control kept arrays on a list, whose reallocation in the shorter run can cancel the longer's growth.
  it("the same step with a leak is caught (the positive control)", () => {
    expect(controlCaught(() => heapScaling(N, leaking(read)))).toBe(true);
  });
});
