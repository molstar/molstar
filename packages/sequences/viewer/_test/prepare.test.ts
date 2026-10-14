/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Prepare (LY2; seqspec:viewer#extension-points, seqspec:sources#levels-and-the-zoom-gate): the kernel types' prepared
// state on fixtures through the graph (the renderers are `PreparedParts` with no draws): features in lanes and their
// segment tables, per-position values by `i − base`, letters `of` other polymers, pairs; and the findings LY2 closes:
// indices, overlapping fragments and one join, a gallop, overlapping polymers, the memory rule for per-position values.

import { createSeqViewSpec, translate } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { case7, case9, case10, case14 } from "../../core/_test/format/cases.js";
import {
  type EmphasisState,
  type Features,
  featureAt,
  type LettersState,
  ofPartAt,
  type PairsState,
  type PerPosition,
} from "../src/prepare.js";
import { segWindow } from "../src/pyramids.js";
import type { RowsState } from "../src/rows.js";
import { layer, opened, preparedRenderers } from "./layout-fixtures.js";
import { allTypes, colorDoc, itemsIn, scaleDoc, transcripts } from "./prepare-fixtures.js";

const P = { renderers: preparedRenderers() };

/** Brute force: items overlapping [from, to) (a boundary: at a position in it). */
const overlaps = (s: number, e: number, from: number, to: number) =>
  s === e ? s >= from && s < to : s < to && e > from;
/** The topmost (last) feature of a segment covering i, by brute force. */
function topAt(f: Features, g: number, i: number): number {
  let at = -1;
  for (let k = f.seg[g]; k < f.seg[g + 1]; k++)
    if ((f.start[k] <= i && f.stop[k] > i) || (f.start[k] === f.stop[k] && f.start[k] === i)) at = k;
  return at;
}

describe("prepare on fixtures", () => {
  it("case 14: the protein's residues over their codons, His12 over its larger part after the intron", async () => {
    const t = await opened(case14(), P);
    const of = (layer(t.ctx, "Protein").state as LettersState).of;
    if (!of) throw new Error("no of");
    // 20 residues as their larger parts, and His12's first base split off by the intron.
    expect([of.count, of.main.filter((m) => m === 1).length]).toEqual([21, 20]);
    const at = (i: number) => of.source.findIndex((x, k) => x === i && of.main[k] === 1);
    expect([of.start[at(0)], of.stop[at(0)], String.fromCharCode(of.code[at(0)])]).toEqual([12, 15, "M"]);
    // His12 (index 11): 45, then [76, 78): drawn over [76, 78); its first base still knows it.
    expect([of.start[at(11)], of.stop[at(11)], String.fromCharCode(of.code[at(11)])]).toEqual([76, 78, "H"]);
    const over = (i: number) => {
      const k = ofPartAt(of, 0, i);
      return k < 0 ? -1 : of.source[k];
    };
    expect([over(45), over(77), over(50)]).toEqual([11, 11, -1]);
    // The domain across the junction in two fragments; the phosphosite one mark over its codon.
    const d = layer(t.ctx, "Domain").state as Features;
    expect(Array.from({ length: d.count }, (_, k) => [d.start[k], d.stop[k]])).toEqual([
      [24, 46],
      [76, 84],
    ]);
    const p = layer(t.ctx, "Phosphosite").state as Features;
    expect([p.count, p.start[0], p.stop[0]]).toEqual([1, 90, 93]);
  });

  it("every type on two chains: lanes, kinds, values, stacks, pairs, emphasis", async () => {
    const t = await opened(allTypes({ length: 300 }), P);
    const f = layer(t.ctx, "Features").state as Features;
    // Lanes never overlap within a block (one fragment an item here).
    for (let k = 1; k < f.count; k++)
      if (f.block[k] === f.block[k - 1] && f.lane[k] === f.lane[k - 1])
        expect(f.start[k]).toBeGreaterThanOrEqual(f.stop[k - 1]);
    expect([f.lanes > 1, f.hasRows, f.split]).toEqual([true, false, null]);
    const ss = layer(t.ctx, "Secondary").state as Features;
    expect(Array.from(ss.kind.slice(0, 3))).toEqual([0, 1, 2]);
    const va = layer(t.ctx, "Variants").state as Features;
    expect(String.fromCharCode(va.code[0])).toBe("A");
    const vals = layer(t.ctx, "Values").state as PerPosition;
    expect([vals.lo, vals.value[0]?.[3], vals.value[1], vals.pyramid[0]?.min.length, vals.pyramid[1]]).toEqual([
      0,
      0.93,
      null,
      8,
      null,
    ]);
    // A `color` layer over the same values: sums and counts only (a `mean` theme's reduction).
    const color = layer(t.ctx, "Color").state as PerPosition;
    expect([color.pyramid[0]?.min.length, color.pyramid[0]?.sum.length]).toEqual([0, 8]);
    const pp = layer(t.ctx, "Contacts").state as PairsState;
    expect(pp.count).toBe(24);
    // The interval index: first keys ascending; each tree node the larger of its children.
    for (let i = 1; i < pp.count; i++) expect(pp.lo[i]).toBeGreaterThanOrEqual(pp.lo[i - 1]);
    for (let i = 1; 2 * i + 1 < pp.tree.length; i++)
      expect(pp.tree[i]).toBe(Math.max(pp.tree[2 * i], pp.tree[2 * i + 1]));
    const em = layer(t.ctx, "Emphasis").state as EmphasisState;
    expect([em.f.count, em.letters.codes.length]).toEqual([40, 2]);
  });

  it("per lane: every prepared feature in view is in its lane's window; `itemsIn` counts them, to a limit", async () => {
    const t = await opened(allTypes({ length: 600 }), P);
    const m = t.ctx.graph.assembled.value?.view;
    if (!m) throw new Error("no view");
    let checked = 0;
    const out = new Int32Array(2);
    for (const p of m.prepared) {
      const f = p.state as Features | null;
      if (!f || !(f.reach instanceof Int32Array)) continue;
      for (let g = 0; g + 1 < f.seg.length; g++)
        for (let from = -5; from < 620; from += 23) {
          segWindow(f.start, f.reach, f.seg[g], f.seg[g + 1], from, from + 40, out);
          for (let k = f.seg[g]; k < f.seg[g + 1]; k++)
            if (overlaps(f.start[k], f.stop[k], from, from + 40)) {
              expect(k >= out[0] && k < out[1]).toBe(true);
              checked++;
            }
        }
      for (let b = 0; b * f.lanes < f.seg.length - 1; b++) {
        let all = 0;
        for (let g = b * f.lanes; g < (b + 1) * f.lanes; g++) {
          segWindow(f.start, f.reach, f.seg[g], f.seg[g + 1], 100, 300, out);
          all += out[1] - out[0];
        }
        expect([itemsIn(f, b, 100, 300, -1, 1e9), itemsIn(f, b, 100, 300, -1, 3)]).toEqual([all, Math.min(3, all)]);
      }
    }
    expect(checked).toBeGreaterThan(100);
  });

  it("dense values: the target's array shared (a subarray from the domain's start), the pyramid's level 0", async () => {
    const t = await opened(scaleDoc({ length: 5000, ranges: 50, variants: 50 }), P);
    const v = layer(t.ctx, "Values");
    const s = v.state as PerPosition;
    const ann = v.m.set.annotations[v.m.resolved.data[v.m.resolved.layers.dataStart[v.k]]];
    const b = s.value.findIndex((x) => x !== null);
    expect([ann?.count, s.value[b]?.buffer === ann?.dense?.[0].buffer, s.annotation[b], s.row[b], s.only[b]]).toEqual([
      0,
      true,
      null,
      null,
      ann?.id,
    ]);
    expect(s.pyramid[b]?.count[0][0]).toBe(4);
  });

  it("color layers: without values no per-position values, its rows say what it covers; with values a value each", async () => {
    const t = await opened(colorDoc(), P);
    const st = layer(t.ctx, "Bands").state as PerPosition;
    expect([st.value, st.pyramid, st.row[0]?.[15], st.row[0]?.[35], st.row[0]?.[25]]).toEqual([
      [null],
      [null],
      0,
      1,
      -1,
    ]);
    const shade = layer(t.ctx, "Shade").state as PerPosition;
    expect([Array.from(shade.value[0] ?? []), shade.row]).toEqual([
      Array.from({ length: 60 }, (_, i) => i / 60),
      [null],
    ]);
    // Letters share the molecule's codes (no copy).
    const letters = layer(t.ctx, "Bands").m.set.molecules[0];
    expect(letters?.code.length).toBe(60);
  });
});

describe("transcripts, fragments and joins (G1d)", () => {
  it("packs lanes by item over its whole extent, joins its fragments across gaps", async () => {
    const t = await opened(transcripts(), P);
    const f = layer(t.ctx, "Genes").state as Features;
    const rows = Array.from({ length: f.count }, (_, k) => [f.start[k], f.lane[k], f.joinFrom[k]]);
    // t1 [10, 140) on lane 0, t2 [40, 110) inside it on lane 1, t3 after t1 on lane 0 again.
    expect(rows).toEqual([
      [10, 0, -1],
      [60, 0, 0],
      [120, 0, 1],
      [150, 0, -1],
      [170, 0, 3],
      [40, 1, -1],
      [90, 1, 5],
    ]);
  });

  /** A transcript as exons and the CDS inside them (one `id`), and an isoform skipping the middle exon. */
  function exonsAndCds() {
    const b = createSeqViewSpec({ title: "Exons and CDS" });
    const a = b.data.sequence({ ref: "a", label: "A", value: "MKTAYIAKQR".repeat(20) });
    const rows = [
      ["tx1", 10, 50], // exon 1, its 5' UTR [10, 30)
      ["tx1", 30, 50], // CDS 1
      ["tx1", 80, 120], // exon 2
      ["tx1", 80, 100], // CDS 2, its 3' UTR [100, 120)
      ["tx1", 150, 170], // exon 3
      ["tx2", 10, 50], // the isoform: exons 1 and 3
      ["tx2", 150, 170],
    ] as const;
    const value = { id: rows.map((r) => r[0]), start: rows.map((r) => r[1]), stop: rows.map((r) => r[2]) };
    const t = b.data.annotation({ ref: "tx", geometry: "ranges", target: a, value });
    b.view({ ref: "v", label: "V" }).block({ sequence: a }).track({ label: "Genes" }).ranges({ data: t });
    return b.getDocument();
  }

  it("an exon under its CDS on one lane: every position hits the one drawn on top, the 3' UTR its exon", async () => {
    const t = await opened(exonsAndCds(), P);
    const f = layer(t.ctx, "Genes").state as Features;
    expect(f.lanes).toBe(2);
    for (let g = 0; g < f.lanes; g++)
      for (let i = 0; i < 200; i++) expect(featureAt(f, 0, g, i), `lane ${g}, ${i}`).toBe(topAt(f, g, i));
    // The 3' UTR of exon 2, after its CDS ends: the exon (the old stack found nothing).
    const k = featureAt(f, 0, 0, 110);
    expect([f.start[k], f.stop[k]]).toEqual([80, 120]);
  });

  it("one hairline across each gap of an item, from the fragment reaching furthest; none for an overlap", async () => {
    const t = await opened(exonsAndCds(), P);
    const f = layer(t.ctx, "Genes").state as Features;
    const joins = [];
    for (let k = 0; k < f.count; k++)
      if (f.joinFrom[k] >= 0) joins.push([f.lane[k], f.stop[f.joinFrom[k]], f.start[k]]);
    expect(joins).toEqual([
      [0, 50, 80],
      [0, 120, 150],
      [1, 50, 150],
    ]);
  });

  it("on a rows track, each row's item joins its own fragments, whatever another row puts between them", async () => {
    const t = await opened(case9(), P);
    const f = layer(t.ctx, "Rows").m.prepared.find((p) => (p.state as Features)?.hasRows)?.state as Features;
    expect(f.count).toBeGreaterThan(0);
    for (let k = 0; k < f.count; k++) {
      const j = f.joinFrom[k];
      if (j < 0) continue;
      expect([f.on[j], f.annotation[j], f.row[j] >= 0]).toEqual([f.on[k], f.annotation[k], true]);
      expect(f.stop[j]).toBeLessThan(f.start[k]);
    }
  });
});

describe("indices are the block coordinate system's; per-position arrays by i − base", () => {
  it("case 10: abl1-kd [241, 493): features and positions in its indices, arrays its length", async () => {
    const doc = case10();
    const t = await opened(doc, P);
    let ranged = 0;
    for (const view of doc.root.children ?? []) {
      if (view.kind !== "view" || !view.ref) continue;
      t.ctx.views.setCurrent(view.ref);
      await t.ctx.graph.update();
      const m = t.ctx.graph.assembled.value?.view;
      if (!m) throw new Error("no view");
      const { sections, blocks, tracks, layers } = m.resolved;
      for (let k = 0; k < layers.count; k++) {
        const st = m.prepared[k]?.state;
        const b0 = sections.blockStart[tracks.section[layers.track[k]]];
        if ((st as Features)?.reach instanceof Int32Array) {
          const f = st as Features;
          for (let i = 0; i < f.count; i++) {
            const b = b0 + f.block[i];
            expect(f.start[i]).toBeGreaterThanOrEqual(blocks.base[b]);
            expect(f.stop[i]).toBeLessThanOrEqual(blocks.base[b] + blocks.length[b]);
            if (blocks.base[b] > 0) ranged++;
          }
        }
        if ((st as PerPosition)?.only instanceof Int32Array)
          (st as PerPosition).value.forEach((v, b) => {
            if (v) expect(v.length).toBe(blocks.length[b0 + b]);
          });
        if ((st as LettersState)?.form instanceof Int8Array)
          (st as LettersState).codes.forEach((c, b) => {
            if (c && !blocks.ligand[b0 + b]) expect(c.length).toBe(blocks.length[b0 + b]);
          });
      }
    }
    expect(ranged).toBeGreaterThan(0);
  });
});

describe("letters `of` overlapping polymers and the cursor", () => {
  function overlapping() {
    const b = createSeqViewSpec();
    const gene = b.data.sequence({ ref: "gene", alphabet: "dna", value: "ATGGCTTCCTAA" });
    b.data.sequence({ ref: "p1", alphabet: "protein", value: "MAS" });
    b.data.sequence({ ref: "p2", alphabet: "protein", value: "AS" });
    b.data.mapping({ ref: "m1", from: "p1", to: "gene", method: "codon", cds: [[0, 9]] } as never);
    b.data.mapping({ ref: "m2", from: "p2", to: "gene", method: "codon", cds: [[3, 9]] } as never);
    const v = b.view({ ref: "v", label: "V" }).block({ sequence: gene });
    v.track({ label: "List" }).letters({ ref: "lst", of: ["p1", "p2"] });
    v.track({ label: "Reversed" }).letters({ ref: "rev", of: ["p2", "p1"] });
    v.track({ label: "Mapped" }).letters({ ref: "mapped", of: "mapped" });
    return b.getDocument();
  }

  it("p1 and p2 over the same codons: the part drawn and hit is the first placement's, at every index", async () => {
    const t = await opened(overlapping(), P);
    const seen: Record<string, string[]> = {};
    for (const label of ["List", "Reversed", "Mapped"]) {
      const l = layer(t.ctx, label);
      const o = (l.state as LettersState).of;
      if (!o) throw new Error("no of");
      expect(o.count).toBe(5);
      seen[label] = [];
      for (let i = 0; i < 12; i++) {
        const k = ofPartAt(o, 0, i);
        // Brute force: the first part covering i, in the parts' (stable) order.
        let want = -1;
        for (let j = o.seg[0]; j < o.seg[1] && want < 0; j++) if (o.start[j] <= i && o.stop[j] > i) want = j;
        expect(k).toBe(want);
        seen[label].push(k < 0 ? "" : l.m.set.strings.get(l.m.set.units.key[o.unit[k]]));
      }
    }
    // Codons 3 to 9 hold both; the first placement's polymer wins there.
    expect(seen.List.slice(3, 9)).toEqual(Array(6).fill(seen.List[0]));
    expect(seen.Reversed.slice(3, 9)).not.toEqual(seen.List.slice(3, 9));
  });

  it("`ofPartAt` galloping from its layer's cursor equals a binary search, for steps, jumps and returns", async () => {
    const b = createSeqViewSpec();
    const length = 30_000;
    const dna = Array.from({ length }, (_, i) => "ACGT"[(i * 7 + (i >> 3)) % 4]).join("");
    const gene = b.data.sequence({ ref: "gene", alphabet: "dna", value: dna });
    const cds: [number, number][] = [];
    for (let s = 30; s + 600 < length; s += 900) cds.push([s, s + 600]);
    const p = translate(b.data, { ref: "p", from: gene, cds });
    b.view({ ref: "v", label: "V" }).block({ sequence: gene }).track({ label: "P" }).letters({ of: p });
    const t = await opened(b.getDocument(), P);
    const o = (layer(t.ctx, "P").state as LettersState).of;
    if (!o) throw new Error("no of");
    expect(o.count).toBeGreaterThan(5000);
    const brute = (i: number) => {
      for (let j = o.seg[0]; j < o.seg[1]; j++) if (o.start[j] <= i && o.stop[j] > i) return j;
      return -1;
    };
    let x = 7;
    const rnd = () => {
      x = (Math.imul(x, 1664525) + 1013904223) >>> 0;
      return x / 2 ** 32;
    };
    const queries: number[] = [];
    for (let i = 0; i < 3000; i++) queries.push(i);
    for (let i = 0; i < 2000; i++) queries.push(Math.floor(rnd() * length));
    for (let i = 0; i < 400; i++) queries.push(5000 + 37 * i);
    for (let i = length - 1; i > length - 400; i--) queries.push(i);
    for (const i of queries) expect(ofPartAt(o, 0, i), `${i}`).toBe(brute(i));
  });
});

describe("split rows", () => {
  it("each value's features a segment table of its own: windows over them hold every overlapping item", async () => {
    const c = await opened(case7(), P);
    const m = c.ctx.graph.assembled.value?.view;
    if (!m) throw new Error("no view");
    const split = m.prepared.map((p) => (p?.state as Features | null)?.split).find((s) => s);
    if (!split) throw new Error("no split layer");
    // Case 7's UniProt features, a row per feature type in data order.
    expect(split.labels).toEqual([
      "Signal",
      "Propeptide",
      "Chain",
      "Active site",
      "Mutagenesis",
      "Beta strand",
      "Helix",
      "Initiator methionine",
    ]);
    const out = new Int32Array(2);
    let checked = 0;
    for (const f of split.rows)
      for (let g = 0; g + 1 < f.seg.length; g++)
        for (let from = 0; from < 800; from += 17) {
          segWindow(f.start, f.reach, f.seg[g], f.seg[g + 1], from, from + 30, out);
          for (let k = f.seg[g]; k < f.seg[g + 1]; k++)
            if (overlaps(f.start[k], f.stop[k], from, from + 30)) {
              expect(k >= out[0] && k < out[1]).toBe(true);
              checked++;
            }
        }
    expect(checked).toBeGreaterThan(10);
    // The compact form is every row's features together.
    const whole = m.prepared.find((p) => (p?.state as Features | null)?.split === split)?.state as Features;
    expect(split.rows.reduce((n, f) => n + f.count, 0)).toBe(whole.count);
  });
});

describe("rows", () => {
  it("case 9: the rows' prepared state has the tree's fields and the header's inset; its features are per row", async () => {
    const t = await opened(case9(), P);
    const r = layer(t.ctx, "Rows");
    const s = r.state as RowsState;
    expect([s.alignment >= 0, s.lazy, s.tree !== null, s.source?.remote !== null]).toEqual([true, true, true, true]);
    const header = t.ctx.registry.get("renderers", "rows")?.header;
    expect(header?.inset(s, {})).toBe(72);
    const f = r.m.prepared.find((p) => (p.state as Features)?.hasRows)?.state as Features;
    // A row's count is its own features (and column data), never the union of every row's.
    const rows = new Set(Array.from(f.on));
    const own = itemsIn(f, 0, 0, 1e6, [...rows][0], 1e9);
    expect(own).toBeLessThan(f.count);
  });
});
