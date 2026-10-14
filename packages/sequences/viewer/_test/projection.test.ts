/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The projector (IN1; seqspec:interaction#channels, seqspec:interaction#selection, seqspec:interaction#focus,
// seqspec:model#projection, seqspec:libraries#translate): loci, `Coordinates` and extents onto every block as a
// MarkSet's runs. Ported first: the old projection.test (case 7's UniProt sequences reach the 1BRS chains the view
// shows, a ref on a shown block marks directly, items, focus) and spans.test (case 14: a residue's codon, a nucleotide
// closing to its residue's codon, His12 split by the intron, items over whole codons, the protein a true translation of
// its CDS), their expected values kept. Then IN1's: the MarkSet equals VW5's core per block, copies, complex runs,
// external systems through a registered translator, ligands in the position form, boundaries, focus
// with links and rows, row rectangles, the lineage guard and `closeLoci`. Ranges are block coordinate system indices.

import {
  type Extents,
  GeneticCode,
  CoordinateSystemKind as K,
  type Loci,
  OrderedSet,
  type SeqLoci,
  type UnitId,
} from "@molstar/seq-core/model";
import { type SystemDescriptor, TranslateMaps, Translator } from "@molstar/seq-core/translate";
import { describe, expect, it } from "vitest";
import { case7, case12, case14 } from "../../core/_test/format/cases.ts";
import { constructs, open as core, gene, rows, runs } from "../../core/_test/format/project-fixtures.ts";
import { MarkBit } from "../src/index.ts";
import { type MarkSet, RowMarks, SectionMarks } from "../src/marks.ts";
import { BlockHulls } from "../src/projection.ts";
import type { Viewport, Viewports } from "../src/viewport.ts";
import { axisOf } from "./axis-fixtures.ts";
import { entries, over, shown } from "./projection-fixtures.ts";

const seqLoci = (set: SeqLoci["set"], unit: number, indices: number[]): Loci =>
  ({
    kind: "seq-loci",
    set,
    elements: [{ unit: unit as UnitId, indices: OrderedSet.ofSortedArray(Int32Array.from(indices)) }],
  }) as Loci;

describe("loci from outside the view mark the shown blocks along mapping paths (ported)", () => {
  it("Complex: barnase marks the merged A, B, C block, barstar the D, E, F one, never the other", async () => {
    const t = await shown(case7(), "complex-one");
    const bn = t.ctx.refs.toLoci("barnase") as SeqLoci;
    // The cause: `barnase` is the UniProt sequence's own unit, which no block of the view shows.
    expect(bn.elements.map((e) => t.set.strings.get(t.set.units.key[e.unit]))).toEqual(["barnase/unit"]);
    // barnase P00648 48–157 is chain A's 0–109 (segments 47 → 0, 110 long).
    expect(t.marked(bn)).toEqual({ "Barnase ×3 (A, B, C)": 110, "Barstar ×3 (D, E, F)": 0 });
    const bs = t.marked(t.ctx.refs.toLoci("barstar") as Loci);
    expect(bs["Barnase ×3 (A, B, C)"]).toBe(0);
    expect(bs["Barstar ×3 (D, E, F)"]).toBeGreaterThan(80);
  });

  it("Each chain: barnase marks A, B and C, barstar D, E and F", async () => {
    const t = await shown(case7(), "complex-each");
    expect(t.marked(t.ctx.refs.toLoci("barnase") as Loci)).toEqual({ A: 110, B: 110, C: 110, D: 0, E: 0, F: 0 });
    const bs = t.marked(t.ctx.refs.toLoci("barstar") as Loci);
    expect([bs.A, bs.B, bs.C]).toEqual([0, 0, 0]);
    expect(bs.D).toBeGreaterThan(80);
    expect([bs.E, bs.F]).toEqual([bs.D, bs.D]);
  });

  it("a ref on a shown block marks directly; copies never reach each other", async () => {
    const t = await shown(case7(), "complex-each");
    const n = t.view.blocks.length[0];
    expect(t.marked(t.ctx.refs.toLoci("1brs/A") as Loci)).toEqual({ A: n, B: 0, C: 0, D: 0, E: 0, F: 0 });
    // A unit merged into a block (B in A, B, C) marks it directly.
    const one = await shown(case7(), "complex-one");
    expect(one.marked(one.ctx.refs.toLoci("1brs/B") as Loci)["Barnase ×3 (A, B, C)"]).toBe(n);
  });

  it("an annotation on the UniProt sequence (its items) marks the chains; the interface marks both sides", async () => {
    const t = await shown(case7(), "complex-each");
    const feats = t.marked(t.ctx.refs.toLoci("feat-P00648") as Loci);
    expect(feats.A).toBeGreaterThan(0);
    expect([feats.A, feats.B]).toEqual([feats.A, feats.C]);
    expect([feats.D, feats.E, feats.F]).toEqual([0, 0, 0]);
    const iface = t.marked(t.ctx.refs.toLoci("iface") as Loci);
    expect(iface.A).toBeGreaterThan(0);
    expect(iface.D).toBeGreaterThan(0);
  });

  it("focus from a UniProt ref lands on the chain it reaches", async () => {
    const t = await shown(case7(), "complex-each");
    expect(t.ctx.focus.fromLoci(t.ctx.refs.toLoci("barnase") as Loci, { transitionMs: 0 })).toBe(true);
    const entries = [...t.ctx.focus.current.value.values()];
    expect(entries.map((e) => [e.block, e.start, e.stop])).toEqual([[0, 0, 110]]);
  });
});

describe("focus over several blocks", () => {
  it("takes every block the target lights, and the viewport frames them all, gaps included", async () => {
    // The Complex view: barnase ×3 and barstar ×3, two blocks side by side; contacts have an end on each.
    const t = await shown(case7(), "complex-one");
    expect(t.ctx.focus.fromLoci(t.ctx.refs.toLoci("contacts") as Loci, { transitionMs: 0 })).toBe(true);
    const [entry] = [...t.ctx.focus.current.value.values()];
    const blocks = Array.from(entry.blocks ?? [], (x) => x);
    expect([entry.block, blocks.length, blocks[0], blocks[3]]).toEqual([0, 6, 0, 1]);
    const ref = t.ctx.viewports.sections.value.find((r) => r.startsWith("complex-one")) as string;
    const vp = t.ctx.viewports.get(ref) as Viewport;
    vp.setSize(900, 400, 1);
    const ax = vp.axis as NonNullable<Viewport["axis"]>;
    const runs = new Int32Array(4);
    /** The view x of block b's position i at the viewport's cell (as the viewport places it). */
    const x = (b: number, i: number, end: boolean) => {
      const n = ax.runToShown(b, i, i + 1, runs);
      const k = end ? runs[2 * n - 1] : runs[0];
      const s = vp.live.current;
      return ax.offset[b] * s.cell + b * ax.gap + k * ax.span[b] * s.cell - s.scrollX;
    };
    const view = () => {
      const w = vp.live.current.width - ax.insetX;
      const [lo, hi] = [x(blocks[0], blocks[1], false), x(blocks[3], blocks[5] - 1, true)];
      return lo >= -0.5 && hi <= w + 0.5 && hi - lo > w / 2;
    };
    vp.zoomToBlocks(entry.blocks as Int32Array);
    expect(view()).toBe(true);
    // One block only (the old behavior) leaves barstar's ends out.
    vp.zoomTo(blocks[0], blocks[1], blocks[2]);
    expect(view()).toBe(false);
    // Narrow: the largest cell that shows it all; narrower than the gap between the blocks: the fit.
    vp.setSize(ax.insetX + 60, 400, 1);
    vp.zoomToBlocks(entry.blocks as Int32Array);
    expect(view()).toBe(true);
    vp.setSize(ax.insetX + ax.gap / 2, 400, 1);
    vp.zoomToBlocks(entry.blocks as Int32Array);
    expect(vp.live.current.cell).toBe(vp.live.current.minCell);
  });
});

describe("a residue's extent on the nucleotide axis is its codon (ported)", () => {
  async function geneView() {
    const t = await shown(case14());
    const [g, p] = [t.unit("gene/unit"), t.unit("protein/unit")];
    const hovered = (unit: number, i: number, closure: "element" | "containing" = "containing") =>
      t.positions("Gene fragment", seqLoci(t.set, unit, [i]), closure);
    return { ...t, gene: g, protein: p, hovered };
  }

  it("hovering Met1 on the protein lights its 3 codon columns; the primary is the codon's first", async () => {
    const t = await geneView();
    expect(t.hovered(t.protein, 0)).toEqual([12, 13, 14]);
    const s = new SectionMarks();
    t.p.project(t.extents(seqLoci(t.set, t.protein, [0])), 0, s.hover);
    s.updatePrimary(axisOf({ set: t.set, view: t.view, issues: [] }), 1);
    expect(s.primary(0)).toBe(12);
  });

  it("hovering any nucleotide of a codon lights its residue, over the whole codon", async () => {
    const t = await geneView();
    for (const i of [12, 13, 14]) expect(t.hovered(t.gene, i)).toEqual([12, 13, 14]);
    // Outside the CDS (the 5′ UTR, the intron, the stop codon) a nucleotide lights alone.
    expect(t.hovered(t.gene, 5)).toEqual([5]);
    expect(t.hovered(t.gene, 60)).toEqual([60]);
    expect(t.hovered(t.gene, 103)).toEqual([103]);
  });

  it("His12, split by the intron, lights all 3 nucleotides from each entry point", async () => {
    const t = await geneView();
    const his = [45, 76, 77];
    expect(t.hovered(t.protein, 11), "the residue").toEqual(his);
    expect(t.hovered(t.gene, 45), "the nucleotide before the intron").toEqual(his);
    expect(t.hovered(t.gene, 76), "the first after it").toEqual(his);
    expect(t.hovered(t.gene, 77), "the second after it").toEqual(his);
    // The intron's edge next to it lights alone.
    expect(t.hovered(t.gene, 46)).toEqual([46]);
  });

  it("a residue's loci and an item on the protein mark whole codons; the domain in two fragments", async () => {
    const t = await geneView();
    // Residues 5–14 (ordinal): [24, 46) and [76, 84), 10 codons.
    expect(t.lit(t.ctx.refs.toLoci("domain") as Loci)).toEqual({
      "Gene fragment": [
        [24, 46],
        [76, 84],
      ],
    });
    expect(t.marked(t.ctx.refs.toLoci("protein") as Loci)).toEqual({ "Gene fragment": 60 });
  });

  it("case 14's data: each residue's codon translates to it; E20's is [99, 102), the stop codon after it", async () => {
    const t = await geneView();
    const dna =
      "GAATTCGCCACCATGGCTTCCAAGGAACTGACCCCACGTGACGGCCGTAAGTCTTTCCTAACTTCTGTTTCCACAGATTGGTACTTCAACAGCCAGGTGGAGTAAGCGGCCGCATAAAGT";
    const protein = "MASKELTPRDGHWYFNSQVE";
    for (let r = 0; r < 20; r++) {
      const [a, b, c] = t.hovered(t.protein, r, "element").map((i) => dna.charCodeAt(i));
      expect(String.fromCharCode(GeneticCode.translate(1, a, b, c)), `residue ${r}`).toBe(protein[r]);
    }
    expect(t.hovered(t.protein, 19, "element")).toEqual([99, 100, 101]);
    expect(
      String.fromCharCode(GeneticCode.translate(1, dna.charCodeAt(102), dna.charCodeAt(103), dna.charCodeAt(104))),
    ).toBe("*");
  });
});

describe("the projector over VW5's core", () => {
  it("per block, a MarkSet holds exactly the core's merged runs, clipped to the block's piece", async () => {
    for (const doc of [gene({ chromosome: true }), gene({ minus: true, filter: true }), constructs(), rows()]) {
      const c = await core(doc);
      const o = over(c.view, c.set);
      const u = c.set.units;
      const sources: [number, number, number, number][] = [];
      for (let k = 0; k < u.count; k++) if (u.alive[k]) sources.push([K.unit, k, -1, u.rangeStop[k]]);
      c.set.alignments.forEach((a, id) => {
        if (!a) return;
        sources.push([K.columns, id, -1, a.columnCount]);
        for (let r = 0; r < a.rows.count; r++) sources.push([K.row, id, r, a.rows.lengthOf(r as never)]);
      });
      for (const [kind, id, row, n] of sources)
        for (let i = 0; i < n; i++) {
          const e = entries([kind, id, row, i, i + 1]);
          for (let s = 0; s < c.view.sections.count; s++) {
            o.p.project(e, s, o.marks);
            for (let b = 0; b < c.view.sections.blockCount[s]; b++) {
              const block = c.view.sections.blockStart[s] + b;
              const [lo, hi] = [c.view.blocks.start[block], c.view.blocks.stop[block]];
              const want = runs(c.p, e, block)
                .map(([x, y]) => [Math.max(x, lo), Math.min(y, hi)])
                .filter(([x, y]) => x < y);
              expect(held(o.marks, b), `${kind}:${id}:${row} ${i} on ${o.label(block)}`).toEqual(want);
            }
          }
        }
    }
  });

  it("one-to-many: a source used twice lights both copies (ub 47 on diub 47 and 123)", async () => {
    const c = await core(constructs());
    const o = over(c.view, c.set);
    const ub = c.system("ub");
    expect(o.lit(entries([ub.kind, ub.id, -1, 47, 48]))).toEqual({
      ub: [[47, 48]],
      diub: [
        [47, 48],
        [123, 124],
      ],
      gfp: [],
      p53: [],
      fusion: [],
    });
  });

  it("complex runs: runs per member unit light a merged block once, and each chain its own", async () => {
    const one = await shown(case7(), "complex-one");
    const [a, b, d] = ["1brs/A", "1brs/B", "1brs/D"].map(one.unit);
    const e = entries([K.unit, a, -1, 0, 10], [K.unit, b, -1, 5, 20], [K.unit, d, -1, 0, 3]);
    expect(one.lit(e)).toEqual({ "Barnase ×3 (A, B, C)": [[0, 20]], "Barstar ×3 (D, E, F)": [[0, 3]] });
    const each = await shown(case7(), "complex-each");
    expect(each.lit(e)).toEqual({ A: [[0, 10]], B: [[5, 20]], C: [], D: [[0, 3]], E: [], F: [] });
  });

  it("duplicates collapse: a nucleotide closed to its residue lights its codon once, with no cap", async () => {
    const t = await shown(case14());
    const e = t.extents(seqLoci(t.set, t.unit("gene/unit"), [77]), "containing");
    expect(e.count).toBeGreaterThan(1);
    expect(t.lit(e)).toEqual({
      "Gene fragment": [
        [45, 46],
        [76, 78],
      ],
    });
    // Every nucleotide closed at once: 120 positions, no cap.
    const all = t.extents(
      seqLoci(
        t.set,
        t.unit("gene/unit"),
        Array.from({ length: 120 }, (_, i) => i),
      ),
      "containing",
    );
    expect(t.marked(all)).toEqual({ "Gene fragment": 120 });
  });

  it("a boundary lights a caret before its position (the end's after the last), no position", async () => {
    const t = await shown(case14());
    const g = t.unit("gene/unit") as UnitId;
    const at = (b: number) => {
      const lit = t.lit({
        kind: "boundary-loci",
        set: t.set,
        elements: [{ unit: g, indices: OrderedSet.ofSingleton(b) }],
      } as Loci);
      return [lit, Array.from(t.marks.caretAt.subarray(0, t.marks.carets))];
    };
    expect(at(20)).toEqual([{ "Gene fragment": [] }, [20]]);
    expect(at(0)).toEqual([{ "Gene fragment": [] }, [0]]);
    expect(at(120)).toEqual([{ "Gene fragment": [] }, [120]]);
  });

  it("a ligand in the position form lights its one shown position for every atom; focus lands on it", async () => {
    const t = await shown(case12());
    const heme = t.ctx.refs.toLoci("dimer.heme") as Loci;
    expect(t.lit(heme)).toEqual({ "Design 7 ×2 (A, B)": [], "dimer.heme": [[0, 1]] });
    const atoms = entries([K.unit, t.unit("dimer.heme"), -1, 5, 9]);
    expect(t.lit(atoms)["dimer.heme"]).toEqual([[0, 1]]);
    const hull = BlockHulls.create();
    expect(t.p.extents(atoms, 0, hull)).toBe(1);
    expect([...hull.data.subarray(0, 3)]).toEqual([1, 0, 1]);
  });

  it("focus with links and rows: a row's residues and a column land on its alignment and its bound chain", async () => {
    const c = await core(rows());
    const o = over(c.view, c.set);
    const row = c.system("a/0");
    const cols = c.system("a");
    const hull = BlockHulls.create();
    const focus = (e: Extents) => [0, 1].map((s) => (o.p.extents(e, s, hull) ? [...hull.data.subarray(0, 3)] : [-1]));
    // a/0 is "MK-VL": residues [1, 3) are columns 1 and 3; its row is bound to x at offset 2.
    expect(focus(entries([row.kind, row.id, 0, 1, 3]))).toEqual([
      [0, 1, 4],
      [0, 3, 5],
    ]);
    expect(focus(entries([cols.kind, cols.id, -1, 3, 5]))).toEqual([
      [0, 3, 5],
      [0, 4, 6],
    ]);
    expect(o.lit(entries([row.kind, row.id, 0, 1, 3]))).toEqual({
      a: [
        [1, 2],
        [3, 4],
      ],
      x: [[3, 5]],
      y: [],
    });
  });

  it("row rectangles: a row's residues on its columns; a band of rows with the same columns stays one", async () => {
    const c = await core(rows());
    const o = over(c.view, c.set);
    const loci = {
      kind: "aligned-loci",
      set: c.set,
      alignment: 0,
      elements: { rows: OrderedSet.ofRange(0, 3), columns: OrderedSet.ofRange(1, 4) },
    } as unknown as Loci;
    const out = new RowMarks();
    o.p.rows(o.extents(loci), 0, MarkBit.selection, out);
    const rects = Array.from({ length: out.count }, (_, r) => [
      out.rowStart[r],
      out.rowStop[r],
      out.start[r],
      out.stop[r],
    ]);
    // Row 0 "MK-VL": K at column 1, V at 3 (its gap at 2 is the column's, lit as a position); rows 1 and 2 [1, 4).
    expect(rects).toEqual([
      [0, 1, 1, 2],
      [0, 1, 3, 4],
      [1, 3, 1, 4],
    ]);
    expect(o.lit(loci).a).toEqual([[1, 4]]);
  });

  it("loci of another set light nothing until carried (lineage)", async () => {
    const t = await shown(case14());
    const other = await shown(case14());
    expect(t.marked(seqLoci(other.set, t.unit("gene/unit"), [3]))).toEqual({ "Gene fragment": 0 });
    expect(t.marked(seqLoci(t.set, t.unit("gene/unit"), [3]))).toEqual({ "Gene fragment": 1 });
  });

  it("closeLoci: a nucleotide closes to itself and its residue, in their own units; element returns it", async () => {
    const t = await shown(case14());
    const loci = seqLoci(t.set, t.unit("gene/unit"), [13]);
    expect(t.p.closeLoci(loci, "element")).toBe(loci);
    const closed = t.p.closeLoci(loci, "containing") as SeqLoci;
    const keyed = closed.elements.map((e) => [
      t.set.strings.get(t.set.units.key[e.unit]),
      OrderedSet.min(e.indices),
      OrderedSet.size(e.indices),
    ]);
    expect(keyed).toEqual([
      ["gene/unit", 13, 1],
      ["protein/unit", 0, 1],
    ]);
  });
});

describe("external systems through a registered translator", () => {
  const descriptor = (key: string, members: string[] | null = null): SystemDescriptor => ({
    key,
    kind: members ? "complex" : "external",
    start: 0,
    stop: 200,
    refs: {},
    members,
    structure: null,
    numberings: [],
    owner: "molstar",
  });

  async function registered() {
    const t = await shown(case7(), "complex-each");
    const tr = new Translator();
    const [a, d] = [t.unit("1brs/A"), t.unit("1brs/D")].map((id) => ({ kind: K.unit, id, row: -1 }));
    // Chain A's residues start 5 before unit A's (a host's own numbering); D's are the same.
    const mapA = (s: { kind: number; id: number; row: number }) => [
      TranslateMaps.bySegments(s, a, Int32Array.of(5, 0, 105)),
    ];
    tr.register(descriptor("3d/A"), mapA);
    tr.register(descriptor("3d/D"), (s) => [TranslateMaps.bySegments(s, d, Int32Array.of(0, 0, 89))]);
    tr.register(descriptor("3d", ["3d/A", "3d/D"]));
    t.p.systems = tr;
    return { ...t, tr, a, mapA };
  }
  const coords = (systems: string[], ...runs: number[]) => ({ systems, runs: Int32Array.from(runs) });

  it("runs per member of a complex go through each member's map, then on the usual route", async () => {
    const t = await registered();
    expect(t.lit(coords(["3d/A", "3d/D"], 0, 10, 20, 1, 0, 3))).toEqual({
      A: [[5, 15]],
      B: [],
      C: [],
      D: [[0, 3]],
      E: [],
      F: [],
    });
    // A run naming the complex itself has no index: nothing.
    expect(t.marked(coords(["3d"], 0, 0, 10))).toEqual({ A: 0, B: 0, C: 0, D: 0, E: 0, F: 0 });
  });

  it("an unregistered key lights nothing meanwhile, and again once it returns", async () => {
    const t = await registered();
    const c = coords(["3d/A"], 0, 10, 12);
    t.tr.unregister("3d/A");
    expect(t.marked(c).A).toBe(0);
    t.tr.register(descriptor("3d/A"), t.mapA);
    expect(t.lit(c).A).toEqual([[5, 7]]);
    expect(t.marked(coords(["nothing"], 0, 0, 5)).A).toBe(0);
  });

  it("a Seq* system registered under its key is read directly; focus takes Coordinates too", async () => {
    const t = await registered();
    t.tr.register({ ...descriptor("seqstar:1brs/B"), kind: "unit", owner: "seqstar" }, [], {
      kind: K.unit,
      id: t.unit("1brs/B"),
      row: -1,
    });
    expect(t.lit(coords(["seqstar:1brs/B"], 0, 2, 4)).B).toEqual([[2, 4]]);
    (t.ctx.viewports as Viewports).projector.systems = t.tr;
    expect(t.ctx.focus.fromLoci(coords(["3d/D"], 0, 7, 9), { transitionMs: 0 })).toBe(true);
    expect([...t.ctx.focus.current.value.values()].map((e) => [e.block, e.start, e.stop])).toEqual([[3, 7, 9]]);
  });

  it("a host's external hover closes like any other: a 3D residue on a chain lights its own position", async () => {
    const t = await registered();
    const e = t.extents(coords(["3d/A"], 0, 30, 31), "containing");
    expect(e.count).toBe(1);
    expect(t.lit(e).A).toEqual([[25, 26]]);
  });
});

/** Block b's runs in a MarkSet, as pairs. */
function held(m: MarkSet, b: number) {
  return Array.from({ length: m.to(b) - m.from(b) }, (_, r) => [m.start[m.from(b) + r], m.stop[m.from(b) + r]]);
}
