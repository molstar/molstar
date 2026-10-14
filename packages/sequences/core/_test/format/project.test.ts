/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The projection core (VW5; seqspec:interaction#channels; seqspec:model#projection): one table of coordinate-map cases
// through one code path (identity, a ranged unit's offset, a codon split by an intron on both strands, a gapped row, a
// filtered block), closure by containment over the whole view (a nucleotide to its residue, isoforms in and out of
// phase, a section with no translation, hidden tracks), one-to-many (linear di-ubiquitin), chimeras, complex runs and
// copy sections, rows through their own binding (a 20k-row alignment adds no routes), granularity and boundaries.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { Extents } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { SelectionGranularity } from "../../src/format/index.js";
import { constructs, type Entry, gene, open, rows, seq } from "./project-fixtures.js";

describe("coordinate maps through one code path", () => {
  it("identity: a molecule's and its unit's runs land on the unit's block as they are", async () => {
    const b = createSeqViewSpec();
    const p = b.data.sequence({ ref: "p", value: seq(30) });
    b.view({ ref: "v", label: "V" }).block({ sequence: p, label: "p" }).track({ label: "L" }).letters();
    const t = await open(b.getDocument());
    expect(t.lit("p", [["p", 3, 6]])).toEqual([[3, 6]]);
    expect(
      t.lit("p", [
        ["p/unit", 3, 6],
        ["p/unit", 6, 9],
        ["p", 20, 40],
      ]),
    ).toEqual([
      [3, 9],
      [20, 30],
    ]);
  });

  it("a ranged unit [5, 15): runs clipped to its domain, its indices the molecule's; boundaries light after", async () => {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: "s", value: seq(20) });
    const u = b.data.unit({ ref: "u", sequence: s, range: [5, 15] });
    b.view({ ref: "v", label: "V" }).block({ unit: u, label: "u" }).track({ label: "L" }).letters();
    const t = await open(b.getDocument());
    expect(t.lit("u", [["s", 0, 20]])).toEqual([[5, 15]]);
    expect(t.lit("u", [["u", 7, 9]])).toEqual([[7, 9]]);
    expect(t.lit("u", [["s", 0, 5]])).toEqual([]);
    // A boundary lights the position after it, the last position at the end.
    expect(t.lit("u", [["u", 7, 7]])).toEqual([[7, 8]]);
    expect(t.lit("u", [["u", 15, 15]])).toEqual([[14, 15]]);
  });

  it("a residue to its codon, H12 split by the intron, on strand + and strand -", async () => {
    const plus = await open(gene());
    expect(plus.lit("gene", [["protein", 0, 1]])).toEqual([[12, 15]]);
    expect(plus.lit("gene", [["protein", 11, 12]])).toEqual([
      [45, 46],
      [76, 78],
    ]);
    expect(plus.lit("gene", [["protein/unit", 11, 12]])).toEqual([
      [45, 46],
      [76, 78],
    ]);
    expect(plus.lit("gene", [["protein", 5, 15]])).toEqual([
      [27, 46],
      [76, 87],
    ]);
    // Strand -: the reverse complement, position i at 119 - i.
    const minus = await open(gene({ minus: true }));
    expect(minus.lit("gene", [["protein", 0, 1]])).toEqual([[105, 108]]);
    expect(minus.lit("gene", [["protein", 11, 12]])).toEqual([
      [42, 44],
      [74, 75],
    ]);
  });

  it("a nucleotide to its residue: the inverse image; outside the CDS nothing", async () => {
    const t = await open(gene());
    for (const i of [45, 76, 77]) expect(t.lit("protein", [["gene/unit", i, i + 1]])).toEqual([[11, 12]]);
    expect(t.lit("protein", [["gene/unit", 5, 6]])).toEqual([]);
    expect(t.lit("protein", [["gene/unit", 60, 61]])).toEqual([]);
    expect(t.lit("protein", [["gene/unit", 103, 104]])).toEqual([]); // the stop codon pairs no residue
    expect(t.lit("protein", [["gene/unit", 40, 80]])).toEqual([[9, 13]]);
  });

  it("a gapped row onto its columns: one run per ungapped run", async () => {
    const b = createSeqViewSpec();
    const a = b.data.alignment({ ref: "a", value: ["MK-VL", "MKIVL"] });
    b.view({ ref: "v", label: "V" }).block({ alignment: a, label: "a" }).track({ label: "L" }).letters();
    const t = await open(b.getDocument());
    expect(t.lit("a", [["a/0", 0, 4]])).toEqual([
      [0, 2],
      [3, 5],
    ]);
    expect(t.lit("a", [["a/1", 2, 3]])).toEqual([[2, 3]]);
    expect(t.lit("a", [["a", 1, 3]])).toEqual([[1, 3]]);
  });

  it("a filtered block: the same block-system runs (the filter is the layout's, never the projection's)", async () => {
    const t = await open(gene({ filter: true }));
    expect(t.view.blocks.positions[t.block("gene")]).not.toBeNull();
    expect(t.lit("gene", [["protein", 11, 12]])).toEqual([
      [45, 46],
      [76, 78],
    ]);
    expect(t.lit("gene", [["gene/unit", 50, 60]])).toEqual([[50, 60]]);
  });
});

describe("inserts and boundaries as carets", () => {
  /** WT and a design with GSGS after 20 and WT 29–31 deleted, through a pairwise mapping; a section each. */
  async function design() {
    const wt = seq(40);
    const d = `${wt.slice(0, 20)}GSGS${wt.slice(20, 29)}${wt.slice(32)}`;
    const b = createSeqViewSpec();
    const w = b.data.sequence({ ref: "wt", value: wt });
    const x = b.data.sequence({ ref: "d", value: d });
    b.data.mapping({
      from: w,
      to: x,
      method: "pairwise",
      gapped: [`${wt.slice(0, 20)}----${wt.slice(20)}`, `${wt.slice(0, 20)}GSGS${wt.slice(20, 29)}---${wt.slice(32)}`],
    });
    // A section each: a source's own block is where its data lands within a section (direct placement).
    const v = b.view({ ref: "v", label: "V" });
    v.section({ ref: "v.wt" }).block({ sequence: w, label: "wt" }).track({ label: "L" }).letters();
    v.section({ ref: "v.d" }).block({ sequence: x, label: "d" }).track({ label: "L" }).letters();
    return open(b.getDocument());
  }
  /** A sink that takes carets: its runs and its carets. */
  const sink = () => {
    const out = { runs: [] as number[][], carets: [] as number[] };
    return { out, add: (s: number, t: number) => out.runs.push([s, t]), boundary: (at: number) => out.carets.push(at) };
  };

  it("an insert lights the boundary it sits at on the other side, a deletion its junction; beyond the pairs nothing", async () => {
    const t = await design();
    const into = (entries: Entry[], label: string) => {
      const k = sink();
      t.p.project(t.extents(entries), t.block(label), k);
      return k.out;
    };
    expect(into([["d", 21, 23]], "wt")).toEqual({ runs: [], carets: [20] });
    expect(into([["wt", 29, 32]], "d")).toEqual({ runs: [], carets: [33] });
    expect(into([["wt", 10, 10]], "d")).toEqual({ runs: [], carets: [10] });
    // A sink without carets (a focus, a host's): the position after a boundary, and nothing for an insert.
    expect(t.lit("d", [["wt", 10, 10]])).toEqual([[10, 11]]);
    expect(t.lit("wt", [["d", 21, 23]])).toEqual([]);
  });

  it("a boundary on either side of a map: a caret on the other, inside an unpaired stretch the next pair's", async () => {
    const t = await design();
    const into = (entries: Entry[], label: string) => {
      const k = sink();
      t.p.project(t.extents(entries), t.block(label), k);
      return k.out.carets;
    };
    // The map's `to` side (the design) onto the wild type: the deletion's junction, the graft's edges and its inside.
    expect(into([["d", 33, 33]], "wt")).toEqual([32]);
    expect(into([["d", 20, 20]], "wt")).toEqual([20]);
    expect(into([["d", 24, 24]], "wt")).toEqual([20]);
    expect(into([["d", 22, 22]], "wt")).toEqual([20]);
    // Its `from` side (the wild type) onto the design: inside the deleted stretch, the junction it leaves.
    expect(into([["wt", 30, 30]], "d")).toEqual([33]);
  });

  it("an intron onto its protein, either strand: the caret after the residue whose codon it splits", async () => {
    const sinkCarets = async (minus: boolean, intron: [number, number]) => {
      const t = await open(gene({ minus }));
      const k = sink();
      t.p.project(t.extents([["gene", intron[0], intron[1]]]), t.block("protein"), k);
      return k.out;
    };
    // Plus: exon 1 [12, 46), exon 2 [76, 105); residue 11's codon is 45 | 76, 77. Minus: read downward, 74 | 43, 42.
    expect(await sinkCarets(false, [50, 60])).toEqual({ runs: [], carets: [12] });
    expect(await sinkCarets(true, [50, 60])).toEqual({ runs: [], carets: [12] });
  });
});

describe("closure by containment over the whole view", () => {
  it("a nucleotide closes to the residue whose codon holds it, which lights the whole codon", async () => {
    const t = await open(gene());
    // Through the letters' map (the protein molecule) and the link's (the protein section's unit).
    expect(t.closed([["gene/unit", 77, 78]])).toEqual([
      ["gene/unit", 77, 78],
      ["protein", 11, 12],
      ["protein/unit", 11, 12],
    ]);
    expect(t.lit("gene", [["gene/unit", 77, 78]], "containing")).toEqual([
      [45, 46],
      [76, 78],
    ]);
    expect(t.lit("gene", [["gene/unit", 77, 78]])).toEqual([[77, 78]]);
    // Outside the CDS a nucleotide closes to nothing; a run closes to every residue it touches.
    expect(t.closed([["gene/unit", 103, 104]])).toEqual([["gene/unit", 103, 104]]);
    expect(t.closed([["gene/unit", 44, 77]])).toEqual([
      ["gene/unit", 44, 77],
      ["protein", 10, 12],
      ["protein/unit", 10, 12],
    ]);
    // A residue closes to the residues over the same codon: here the same residue in the protein section's unit.
    expect(t.closed([["protein", 11, 12]])).toEqual([
      ["protein", 11, 12],
      ["protein/unit", 11, 12],
    ]);
    // The protein section lights a hovered nucleotide's residue either way.
    expect(t.lit("protein", [["gene/unit", 77, 78]], "containing")).toEqual([[11, 12]]);
  });

  it("isoforms: one in the same phase is contained, one in another phase is not; one step", async () => {
    const b = createSeqViewSpec();
    const g = b.data.sequence({ ref: "g", alphabet: "dna", value: "ATGGCTTCCAAGGAAT" });
    b.data.sequence({ ref: "p1", alphabet: "protein", value: "MASKE" });
    b.data.sequence({ ref: "p2", alphabet: "protein", value: "ASK" }); // [3, 12), the same phase as p1
    b.data.sequence({ ref: "p3", alphabet: "protein", value: "WLP" }); // [1, 10), another phase
    b.data.mapping({ ref: "m1", from: "p1", to: "g", method: "codon", cds: [[0, 15]] } as never);
    b.data.mapping({ ref: "m2", from: "p2", to: "g", method: "codon", cds: [[3, 12]] } as never);
    b.data.mapping({ ref: "m3", from: "p3", to: "g", method: "codon", cds: [[1, 10]] } as never);
    b.view({ ref: "v", label: "V" }).block({ sequence: g, label: "g" }).track({ label: "P" }).letters({ of: "mapped" });
    const t = await open(b.getDocument());
    // p1's A (3–5) is p2's A: contained; p3's residues overlap it but hold other triplets.
    expect(t.closed([["p1", 1, 2]])).toEqual([
      ["p1", 1, 2],
      ["p2", 0, 1],
    ]);
    expect(t.closed([["p2", 0, 1]])).toEqual([
      ["p1", 1, 2],
      ["p2", 0, 1],
    ]);
    // A nucleotide is contained by a residue of each: p1's M and p3's first (1–3) hold nucleotide 2.
    expect(t.closed([["g/unit", 2, 3]])).toEqual([
      ["g/unit", 2, 3],
      ["p1", 0, 1],
      ["p3", 0, 1],
    ]);
    // One step: p3's residue 0, added for nucleotide 2, is not closed again (it would add p1's A).
    expect(t.lit("g", [["g/unit", 2, 3]], "containing")).toEqual([[0, 4]]);
  });

  it("the view's systems, not the hit section's: a chromosome section with no translation closes through the gene", async () => {
    const t = await open(gene({ chromosome: true }));
    expect(t.closed([["chr/unit", 177, 178]])).toEqual([
      ["chr/unit", 177, 178],
      ["protein", 11, 12],
      ["protein/unit", 11, 12],
    ]);
    expect(t.lit("chr", [["chr/unit", 177, 178]], "containing")).toEqual([
      [145, 146],
      [176, 178],
    ]);
    expect(t.lit("gene", [["chr/unit", 177, 178]], "containing")).toEqual([
      [45, 46],
      [76, 78],
    ]);
    // The same element lights the same everywhere, whichever section it was hovered in.
    expect(t.lit("chr", [["gene/unit", 77, 78]], "containing")).toEqual([
      [145, 146],
      [176, 178],
    ]);
  });

  it("hidden tracks give no closure maps: the nucleotide stays alone; the protein section still links", async () => {
    const t = await open(gene({ hideTrack: true }));
    // The link to the protein section is a span map too, so a nucleotide still closes through it.
    expect(t.closed([["gene/unit", 77, 78]])).toEqual([
      ["gene/unit", 77, 78],
      ["protein/unit", 11, 12],
    ]);
    const b = createSeqViewSpec();
    const g = b.data.sequence({ ref: "g", alphabet: "dna", value: "ATGGCTTCCTAA" });
    b.data.sequence({ ref: "p", alphabet: "protein", value: "MAS" });
    b.data.mapping({ ref: "m", from: "p", to: "g", method: "codon", cds: [[0, 9]] } as never);
    b.view({ ref: "v", label: "V" })
      .block({ sequence: g, label: "g" })
      .track({ label: "P", hidden: true })
      .letters({ of: "p" });
    const h = await open(b.getDocument());
    expect(h.closed([["g/unit", 4, 5]])).toEqual([["g/unit", 4, 5]]);
  });
});

describe("one-to-many, chimeras, complexes and copies", () => {
  it("linear di-ubiquitin: Ub 48 lights both copies; each copy lights Ub 48", async () => {
    const t = await open(constructs());
    expect(t.lit("diub", [["ub", 47, 48]])).toEqual([
      [47, 48],
      [123, 124],
    ]);
    expect(t.lit("ub", [["diub", 123, 124]])).toEqual([[47, 48]]);
    expect(t.lit("ub", [["diub", 70, 80]])).toEqual([
      [0, 4],
      [70, 76],
    ]);
    expect(t.lit("diub", [["ub/unit", 0, 76]])).toEqual([[0, 152]]);
  });

  it("chimeras: each part on its source, the linker on neither, a run across the junction clipped", async () => {
    const t = await open(constructs());
    expect(t.lit("fusion", [["gfp", 10, 11]])).toEqual([[10, 11]]);
    expect(t.lit("fusion", [["p53", 5, 6]])).toEqual([[248, 249]]);
    expect([t.lit("gfp", [["fusion", 240, 241]]), t.lit("p53", [["fusion", 240, 241]])]).toEqual([[], []]);
    expect(t.lit("gfp", [["fusion", 230, 250]])).toEqual([[230, 238]]);
    expect(t.lit("p53", [["fusion", 230, 250]])).toEqual([[0, 7]]);
    expect(t.lit("p53", [["gfp", 0, 238]])).toEqual([]);
  });

  /** Two copies of s (merged at first, then each) and one unit of another molecule. */
  function complex() {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: "s", value: seq(40) });
    const r = b.data.sequence({ ref: "r", value: seq(30, 3) });
    const [a1, a2] = ["a1", "a2"].map((ref) => b.data.unit({ ref, sequence: s }));
    const d = b.data.unit({ ref: "d", sequence: r });
    const v = b.view({ ref: "v", label: "V" });
    for (const u of [a1, a2, d]) v.block({ unit: u, label: u === a1 ? "A (merged)" : u === a2 ? "A2" : "D" });
    v.track({ label: "L" }).letters();
    const w = v.section({ ref: "v.each", copies: "each" });
    w.block({ unit: a1, label: "a1" });
    w.block({ unit: a2, label: "a2" });
    w.block({ unit: d, label: "d" });
    return b.getDocument();
  }

  it("complex runs: runs per member unit, each through its own routes; copies never reach each other", async () => {
    const t = await open(complex());
    const ext: Entry[] = [
      ["a2", 2, 5],
      ["d", 0, 3],
    ];
    expect(t.lit("A (merged)", ext)).toEqual([[2, 5]]); // a2 is a member of the merged block
    expect(t.lit("D", ext)).toEqual([[0, 3]]);
    expect([t.lit("a1", ext), t.lit("a2", ext), t.lit("d", ext)]).toEqual([[], [[2, 5]], [[0, 3]]]);
  });

  it("two copy sections both lit by their molecule", async () => {
    const t = await open(complex());
    expect([t.lit("a1", [["s", 2, 5]]), t.lit("a2", [["s", 2, 5]]), t.lit("d", [["s", 2, 5]])]).toEqual([
      [[2, 5]],
      [[2, 5]],
      [],
    ]);
  });
});

describe("rows go through their own row only", () => {
  it("a bound row through its binding; bound twice or unbound lights nothing elsewhere; its own columns always", async () => {
    const t = await open(rows());
    expect(t.lit("x", [["a/0", 0, 4]])).toEqual([[2, 6]]);
    expect(t.lit("a", [["a/0", 0, 4]])).toEqual([
      [0, 2],
      [3, 5],
    ]);
    expect([t.lit("x", [["a/1", 0, 2]]), t.lit("y", [["a/1", 0, 2]])]).toEqual([[], []]);
    expect(t.lit("x", [["a/2", 0, 4]])).toEqual([]);
    expect(t.lit("a", [["a/2", 0, 4]])).toEqual([[0, 4]]);
    // A column goes through the rows: row 0 holds column 3 at its residue 2, x's 4.
    expect(t.lit("x", [["a", 3, 4]])).toEqual([[4, 5]]);
  });

  it("direct wins: a chain lights its homolog through the alignment only in a section that does not show it", async () => {
    const b = createSeqViewSpec();
    b.data.alignment({ ref: "a", value: ["MKVL", "MRVL"] });
    for (const ref of ["x", "y"]) b.data.sequence({ ref, value: "GGMKVLGG" });
    for (const [row, to] of [
      ["a/0", "x"],
      ["a/1", "y"],
    ])
      b.data.mapping({ ref: `m${to}`, from: row, to, method: "segments", segments: [[0, 2, 4]] } as never);
    const v = b.view({ ref: "v", label: "V" });
    v.block({ sequence: "x", label: "x" });
    v.block({ sequence: "y", label: "y" }).track({ label: "L" }).letters();
    v.section({ ref: "v.y" }).block({ sequence: "y", label: "y alone" });
    const t = await open(b.getDocument());
    // x → a/0 → columns → a/1 → y: a path, but x's molecule is shown in the first section (seqspec:format#placement).
    expect(t.lit("y", [["x", 3, 4]])).toEqual([]);
    expect(t.lit("y alone", [["x", 3, 4]])).toEqual([[3, 4]]);
  });

  it("20,000 rows hovered one by one add no routes beyond their bound unit's", async () => {
    const t = await open(rows(20_000));
    const ext = Extents.create(1);
    const blocks = [t.block("x"), t.block("y"), t.block("a")];
    const sink = { n: 0, add() {} };
    let routes = -1;
    for (let r = 0; r < 20_000; r++) {
      ext.count = 0;
      Extents.add(ext, { kind: 2, id: t.system("a").id, row: r }, 1, 2);
      for (const b of blocks) t.p.project(ext, b, sink);
      if (r === 1) routes = t.p.routes;
    }
    expect(t.p.routes).toBe(routes);
    expect(routes).toBeLessThanOrEqual(3);
  }, 60_000);
});

describe("granularity after closure", () => {
  it("unit, copy-group, column and row widen each entry in its own system; position and range change nothing", async () => {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: "s", value: seq(40) });
    for (const ref of ["a1", "a2"]) b.data.unit({ ref, sequence: s, range: [5, 35] });
    const a = b.data.alignment({ ref: "a", value: ["MK-VL", "MKIVL"] });
    const v = b.view({ ref: "v", label: "V" });
    v.block({ unit: "a1", label: "a1" }).track({ label: "L" }).letters();
    v.section({ ref: "v.a" }).block({ alignment: a, label: "a" });
    const t = await open(b.getDocument());
    const widened = (g: number, entries: Entry[]) => {
      const e = t.extents(entries);
      t.p.widen(e, g);
      return t.keyed(e);
    };
    const G = SelectionGranularity;
    expect(widened(G.position, [["a1", 7, 9]])).toEqual([["a1", 7, 9]]);
    expect(widened(G.range, [["a1", 7, 9]])).toEqual([["a1", 7, 9]]);
    expect(
      widened(G.unit, [
        ["a1", 7, 9],
        ["a/0", 1, 2],
      ]),
    ).toEqual([
      ["a1", 5, 35],
      ["a/0", 0, 4],
    ]);
    expect(widened(G.copyGroup, [["a1", 7, 9]])).toEqual([
      ["a1", 7, 9],
      ["a2", 7, 9],
    ]);
    expect(
      widened(G.column, [
        ["a/0", 1, 3],
        ["a1", 7, 9],
      ]),
    ).toEqual([
      ["a1", 7, 9],
      ["a columns", 1, 2],
      ["a columns", 3, 4],
    ]);
    // CF-15: adjacent columns are one run, not an extent each.
    expect(widened(G.column, [["a/1", 1, 4]])).toEqual([["a columns", 1, 4]]);
    expect(
      widened(G.row, [
        ["a/1", 2, 3],
        ["a", 2, 3],
      ]),
    ).toEqual([
      ["a/1", 0, 5],
      ["a columns", 2, 3],
    ]);
  });
});
