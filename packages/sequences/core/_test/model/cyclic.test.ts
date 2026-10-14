/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Cyclic coordinate systems (seqspec:model#cyclic-systems; CY1): a cyclic molecule, its arithmetic (`cyclicDistance`,
// `cyclicRuns`, `cyclicHull` against a brute force), annotation rows across the origin stored as two runs of one item,
// labels that wrap, mapping segments and CDSs split at the origin, and what stays linear (a ligand, a pair, a unit's
// range).

import { describe, expect, it } from "vitest";
import { type Annotation, Coords, type Mapping, type SeqSet, SeqSets } from "../../src/model/index.ts";

const COMPLEMENT: Record<string, string> = { A: "T", T: "A", C: "G", G: "C" };
const revcomp = (s: string) =>
  [...s]
    .reverse()
    .map((c) => COMPLEMENT[c] ?? c)
    .join("");
const codes = (issues: readonly { code: string }[]) => issues.map((i) => i.code);
const rows = (a: Annotation) => Array.from(a.start, (s, r) => [s, a.stop[r], a.item[r]]);

/** A 40-nucleotide cyclic molecule "ring" and its own unit, with the given annotations and mappings. */
function ring(extra: Partial<Parameters<typeof SeqSets.build>[0]> = {}) {
  const residues = "ACGTACGTTA".repeat(4);
  return SeqSets.build({
    ...extra,
    molecules: [
      { key: "ring", alphabet: "dna", residues, cyclic: true, ownUnit: true },
      { key: "line", alphabet: "dna", residues, ownUnit: true },
      ...(extra.molecules ?? []),
    ],
  });
}

describe("cyclic coordinate systems", () => {
  it("is a polymer's flag: a ligand's never", () => {
    const { set, issues } = ring({ molecules: [{ key: "x", kind: "ligand", elements: ["C", "O"], cyclic: true }] });
    expect(issues).toEqual([]);
    expect(set.molecules.map((m) => m?.cyclic)).toEqual([true, false, false]);
  });

  it("measures distance the short way round", () => {
    expect([0, 5, 20, 35, 39].map((b) => Coords.cyclicDistance(2, b, 40))).toEqual([2, 3, 18, 7, 3]);
    expect(Coords.cyclicDistance(39, 0, 40)).toBe(1);
  });

  it("splits a range across the origin into its two runs, an empty one left out", () => {
    const out = new Int32Array(4);
    const runs = (s: number, t: number) => Array.from(out.subarray(0, 2 * Coords.cyclicRuns(s, t, 40, out)));
    expect([runs(5, 9), runs(36, 4), runs(36, 0)]).toEqual([
      [5, 9],
      [36, 40, 0, 4],
      [36, 40],
    ]);
  });

  it("gives the shortest arc holding the runs, the circle less its largest gap (against a brute force)", () => {
    let seed = 7;
    const next = (n: number) => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed % n;
    };
    const out = new Int32Array(4);
    for (let trial = 0; trial < 500; trial++) {
      const L = 2 + next(30);
      const runs: number[] = [];
      for (let k = 1 + next(4); k > 0; k--) {
        const s = next(L);
        runs.push(s, s + 1 + next(L - s));
      }
      const pairs = Array.from({ length: runs.length / 2 }, (_, k) => [runs[2 * k], runs[2 * k + 1]]);
      pairs.sort((a, b) => a[0] - b[0]);
      const covered = new Uint8Array(L);
      for (const [s, t] of pairs) covered.fill(1, s, t);
      // The brute force: every arc [a, a + w) holding every covered position, the shortest, the one at 0 first.
      let best = [0, L];
      for (let w = 1; w <= L && best[1] === L; w++)
        for (let a = 0; a < L; a++) {
          const holds = covered.every((c, i) => !c || (i - a + L) % L < w);
          if (holds) {
            best = [a, w];
            break;
          }
        }
      const n = Coords.cyclicHull(Int32Array.from(pairs.flat()), pairs.length, L, out);
      const got = Array.from(out.subarray(0, 2 * n));
      const width = n === 2 ? got[1] - got[0] + got[3] : got[1] - got[0];
      expect(width, `trial ${trial}: ${pairs} on ${L}`).toBe(best[1]);
      // Every covered position lies in the arc; across the origin only when no linear hull is as short.
      for (let i = 0; i < L; i++)
        if (covered[i]) expect(got.some((x, k) => k % 2 === 0 && i >= x && i < got[k + 1])).toBe(true);
      if (n === 2) expect(got[0] > got[3] && got[1] === L && got[2] === 0).toBe(true);
    }
    expect(Coords.cyclicHull(new Int32Array(0), 0, 40, out)).toBe(0);
  });

  it("stores an annotation row across the origin as two runs of one item, by index and by label", () => {
    const { set, issues } = ring({
      annotations: [
        {
          key: "f",
          geometry: "ranges",
          columns: {
            target: ["ring", "ring/unit", "ring", "ring"],
            start: [36, 30, 2, 38],
            stop: [4, 35, 6, 0],
            id: [null, null, "ori", null],
          },
        },
        { key: "lab", geometry: "ranges", target: "ring", numbering: "ordinal", columns: { beg: [39], end: [2] } },
        { key: "v", geometry: "values", target: "ring", columns: { start: [38], stop: [1], value: [0.5] } },
      ],
    });
    expect(issues).toEqual([]);
    const [f, lab, v] = set.annotations as Annotation[];
    // Row 0 wraps (two runs, item 0); row 1 is plain; row 3, [38, 0), is [38, 40) alone.
    expect(rows(f)).toEqual([
      [36, 40, 0],
      [0, 4, 0],
      [30, 35, 1],
      [2, 6, 2],
      [38, 40, 3],
    ]);
    expect(f.itemCount).toBe(4);
    expect(rows(lab)).toEqual([
      [38, 40, 0],
      [0, 2, 0],
    ]);
    expect([rows(v), Array.from(v.sparse ?? [])]).toEqual([
      [
        [38, 40, 0],
        [0, 1, 0],
      ],
      [0.5, 0.5],
    ]);
  });

  it("keeps a linear molecule's, and a pair's, start > stop an error that says what wraps", () => {
    const { issues } = ring({
      annotations: [
        { key: "f", geometry: "ranges", target: "line", columns: { start: [36], stop: [4] } },
        {
          key: "p",
          geometry: "pairs",
          target: "ring",
          columns: { start: [36], stop: [4], b_start: [10], b_stop: [11] },
        },
      ],
    });
    expect(codes(issues)).toEqual(["model.position.range", "model.position.range"]);
    expect(issues[0].message).toContain("start > stop wraps only a cyclic sequence's origin");
  });

  it("does not yet let a unit's range wrap: a block's range shows the wrap", () => {
    const { issues } = ring({ units: [{ key: "u", molecule: "ring", range: [30, 5] }] });
    expect(codes(issues)).toEqual(["model.unit.range"]);
    expect(issues[0].message).toContain("a block's range");
  });

  it("splits a mapping's segment at a cyclic origin, either side", () => {
    // "frag" is the ring's 12 positions from 34 round to 6.
    const residues = "ACGTACGTTA".repeat(4);
    const frag = residues.slice(34) + residues.slice(0, 6);
    const { set, issues } = ring({
      molecules: [{ key: "frag", alphabet: "dna", residues: frag, ownUnit: true }],
      mappings: [
        { key: "in", from: "frag", to: "ring", method: "segments", segments: [[0, 34, 12]] },
        { key: "out", from: "ring", to: "frag", method: "segments", segments: [[34, 0, 12]] },
        { key: "bad", from: "frag", to: "line", method: "segments", segments: [[0, 34, 12]] },
      ],
    });
    expect(codes(issues)).toEqual(["model.mapping.range"]);
    const [a, b] = set.mappings as Mapping[];
    const segments = (m: Mapping) => Array.from(m.map.at, (at, k) => [at, m.map.length[k], m.map.toLo[k]]);
    expect(segments(a)).toEqual([
      [0, 6, 34],
      [6, 6, 0],
    ]);
    expect(segments(b)).toEqual([
      [0, 6, 6],
      [34, 6, 0],
    ]);
    expect([Array.from(a.status).every((s) => s === 0), Array.from(b.forward.subarray(34, 40))]).toEqual([
      true,
      [0, 1, 2, 3, 4, 5],
    ]);
  });

  it("reads a CDS across a cyclic origin, from two ranges or one start, on both strands", () => {
    const cds = "ATGAAATTTGGG"; // M K F G
    const plus = `${cds.slice(5)}${"C".repeat(10)}${cds.slice(0, 5)}`; // the CDS from 17 round to 7
    const minus = revcomp(plus); // the CDS from 4 down round to 15
    const build = (dna: string, m: object) =>
      SeqSets.build({
        molecules: [
          { key: "p", alphabet: "dna", residues: dna, cyclic: true, ownUnit: true },
          { key: "prot", alphabet: "protein", residues: "MKFG", ownUnit: true },
        ],
        mappings: [{ key: "c", from: "prot", to: "p", method: "codon", ...m } as never],
      });
    const cases: [string, object, number][] = [
      [
        plus,
        {
          cds: [
            [17, 22],
            [0, 7],
          ],
        },
        17,
      ],
      [plus, { start: 17 }, 17],
      [
        minus,
        {
          cds: [
            [0, 5],
            [15, 22],
          ],
          strand: "-",
        },
        4,
      ],
      [minus, { start: 15, strand: "-" }, 4],
    ];
    for (const [dna, m, first] of cases) {
      const { set, issues } = build(dna, m);
      expect(issues, JSON.stringify(m)).toEqual([]);
      const c = set.mappings[0] as Mapping;
      expect([Array.from(c.status), Coords.toAt(c.map, 0, 0), c.map.count]).toEqual([[0, 0, 0, 0], first, 2]);
    }
    // Across the origin only once, and never back over the first range; a linear `to` not at all.
    expect(
      codes(
        build(plus, {
          cds: [
            [17, 22],
            [0, 4],
            [2, 5],
          ],
        }).issues,
      ),
    ).toEqual(["model.mapping.codon"]);
    const linear = (set: SeqSet) => set.molecules[0]?.cyclic;
    const { set, issues } = SeqSets.build({
      molecules: [
        { key: "p", alphabet: "dna", residues: plus, ownUnit: true },
        { key: "prot", alphabet: "protein", residues: "MKFG", ownUnit: true },
      ],
      mappings: [
        {
          key: "c",
          from: "prot",
          to: "p",
          method: "codon",
          cds: [
            [17, 22],
            [0, 7],
          ],
        },
      ],
    });
    expect([linear(set), codes(issues)]).toEqual([false, ["model.mapping.codon"]]);
  });
});
