/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for translate's hot paths (TR1; seqspec:libraries#translate): `forEachAt` at pointer rate and `translate`
// in batches, through cached routes (one map, composed with a codon step either way, copies, complex members), on
// RB0b's stepped guard: runs of `n` and `4n` steps, the longer growing at most 64 KB beyond the shorter, and no GC in
// either.

import { type CoordinateSystem, CoordinateSystemKind, Coords, Extents } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { type SystemDescriptor, TranslateMaps, Translator } from "../../src/translate/index.ts";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from "../model/heap.ts";

/** Steps in the shorter run; the longer runs four times as many. */
const N = 250_000;
const K = CoordinateSystemKind;
const desc = (key: string, stop: number, members: string[] | null = null): SystemDescriptor => ({
  key,
  kind: members ? "complex" : "external",
  start: 0,
  stop,
  refs: {},
  members,
  structure: null,
  numberings: [],
  owner: "test",
});

const t = new Translator();
const [p, g, ub, diub] = [0, 1, 2, 3].map((id): CoordinateSystem => ({ kind: K.molecule, id, row: -1 }));
const [ua, ubu] = [1, 2].map((id): CoordinateSystem => ({ kind: K.unit, id, row: -1 }));
// A 300-residue protein over a gene in three exons, and a host system onto the protein.
t.addMap(Coords.create(p, g, 3, [0, 300, 50, 1, 300, 300, 500, 1, 600, 300, 1000, 1]));
const x = t.register(desc("host:x", 300), (s) => [Coords.create(s, p, 1, [0, 300, 0, 1])]);
// A source used twice.
t.addMap(Coords.create(diub, ub, 1, [0, 76, 0, 1, 76, 76, 0, 1]));
// A structure of two chains onto two units of one complex.
const labels = Int32Array.from({ length: 200 }, (_, i) => i + 1);
t.register(desc("seqstar:A", 200), [], ua);
t.register(desc("seqstar:B", 200), [], ubu);
const seq = t.register(desc("seqstar:1brs", 0, ["seqstar:A", "seqstar:B"]), [], { kind: K.complex, id: 0, row: -1 });
const chain = (key: string, u: CoordinateSystem) =>
  t.register(desc(key, 190), (s) => [TranslateMaps.byLabels(s, labels.subarray(10), u, labels)]);
const [ca, cb] = [chain("molstar:3d/A@1", ua), chain("molstar:3d/B@1", ubu)];
const structure = t.register(desc("molstar:3d", 0, ["molstar:3d/A@1", "molstar:3d/B@1"]));

/** What the steps visited or translated: a guard over steps that do nothing would pass blind (CB-05). */
let visits = 0;
const count = (i: number, c: { n: number }) => {
  c.n = (c.n + i) & 0xfffff;
};
const acc = { n: 0 };

function guard(name: string, n: number, step: (i: number) => void) {
  it(name, () => {
    visits = 0;
    const r = heapScaling(n, step);
    expect(visits, "the steps visit images").toBeGreaterThan(n);
    expect(isFlat(r), describeScaling(r, n)).toBe(true);
    // The template's control (CB-05): the same step with a leak is caught.
    expect(controlCaught(() => heapScaling(n, leaking(step)))).toBe(true);
  });
}

describe("translate heap guards", () => {
  guard("forEachAt through a composed codon route, forward and backward", N, (i) => {
    visits += t.forEachAt(x, i % 300, g, count, acc);
    visits += t.forEachAt(g, 50 + (i % 1200), x, count, acc);
  });

  guard("forEachAt visiting copies, and complex members", N, (i) => {
    visits += t.forEachAt(ub, i % 76, diub, count, acc);
    visits += t.forEachAt(ua, i % 200, structure, count, acc);
  });

  const runs = Extents.create(4);
  Extents.add(runs, ca, 0, 40);
  Extents.add(runs, cb, 100, 150);
  Extents.add(runs, ca, 170, 190);
  const out = Extents.create(4);
  const gene = Extents.create(4);
  Extents.add(gene, g, 40, 900);

  guard("translate complex runs across units, and backward through a codon step", N / 4, () => {
    out.count = 0;
    t.translate(runs, structure, seq, out);
    t.translate(gene, g, x, out);
    visits += out.count;
  });

  it("translated what it should", () => {
    out.count = 0;
    t.translate(runs, structure, seq, out);
    expect(out.count).toBe(3);
  });
});
