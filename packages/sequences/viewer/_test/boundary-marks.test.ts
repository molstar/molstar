/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Boundaries (an insertion's place) as `ranges`, `points` and `variants` draw and hit them: on the edge between two
// positions, a boundary after the last position (a C-terminal insertion) on that position's trailing edge, a boundary
// kept whole by the per-position types (a mark between positions, not on the position after); on a reversed piece
// the edge on the right; a junction between two pieces of one unit once; before a position a filter hides, after the
// position before it; and an insert placed through a partial map on the boundary after it (seqspec:format#placement).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { Hit, HoverKind } from "../src/index.ts";
import { recording } from "./composer-fixtures.ts";
import { lit } from "./overlay-fixtures.ts";

const S = "v.main";
const N = 40;

type Block = { range?: [number, number]; direction?: 1 | -1 };
function doc(blocks: Block[] = [{}], positions = false) {
  const b = createSeqViewSpec();
  const p = b.data.sequence({ ref: "p", value: "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRV".slice(0, N) });
  const ends = b.data.annotation({
    ref: "ends",
    geometry: "ranges",
    target: p,
    value: { boundary: [0, 20, N], label: ["start", "middle", "end"] },
  });
  // Every position but 20, for the filter.
  const kept = b.data.annotation({
    ref: "kept",
    geometry: "ranges",
    target: p,
    value: { start: [0, 21], stop: [20, N] },
  });
  const sec = b.view({ ref: "v", label: "V" }).section({ ref: S, ...(positions ? { positions: kept } : {}) });
  for (const o of blocks) sec.block({ sequence: p, ...o });
  sec.track({ ref: "t-ranges", label: "Ranges", height: 2 }).ranges({ data: ends });
  sec.track({ ref: "t-points", label: "Points", height: 2 }).points({ data: ends });
  sec.track({ ref: "t-variants", label: "Variants", height: 2 }).variants({ data: ends });
  return b.getDocument();
}

/** The glyphs' centers a frame drew (device px = CSS px at dpr 1), rounded. */
const glyphXs = (be: ReturnType<typeof recording>) =>
  Array.from(be.glyphItems.x.subarray(0, be.glyphItems.n), (x) => Math.round(x));

async function open(blocks?: Block[], positions?: boolean) {
  const t = await lit(doc(blocks, positions), { width: 900, height: 400 });
  // 15 px a position: the whole sequence in view.
  t.vp?.setView(15, 0, "api");
  // The first frame's glyphs (a frame with nothing changed since draws nothing).
  const xs = glyphXs(t.frame(recording()));
  const ax = () => t.ctx.layout.axis(S);
  /** The left edge of position i (CSS px); i = N: the last one's right edge. */
  const edge = (i: number) => {
    const a = ax();
    return i < N ? (a?.x(0, a.shownOf(0, i)) ?? 0) : (a?.x(0, (a.shownOf(0, N - 1) + 1) as never) ?? 0);
  };
  const hit = Hit.create();
  const at = (x: number, y: number) => {
    t.ctx.layout.hitTest(S, x, y, hit);
    return hit.kind === HoverKind.item ? `item ${hit.item}` : hit.kind === HoverKind.seq ? `seq ${hit.index}` : "none";
  };
  /** What a column `x` hits down a track's band, each once. */
  const down = (x: number, track: string) => {
    const l = t.layout(S);
    const e = Array.from({ length: l.entryCount }, (_, k) => k).find((k) => l.entryTrack[k] === t.trackIndex(track));
    const top = l.scrolledTop + l.entryTop[e as number];
    const seen = new Set<string>();
    for (let y = top; y < top + l.entryHeight[e as number]; y += 0.5) seen.add(at(x, y));
    return seen;
  };
  return { ...t, edge, at, xs, down };
}

describe("boundary marks", () => {
  it("ranges: a boundary after the last position draws and hits on its trailing edge", async () => {
    const t = await open();
    const y = t.yOf(S, "t-ranges");
    expect([t.at(t.edge(20) - 2, y), t.at(t.edge(20) + 2, y)]).toEqual(["item 1", "item 1"]);
    expect(t.at(t.edge(N) - 2, y)).toBe("item 2");
  });

  it("points: a boundary's glyph on its edge, not on the position after it; one after the last position too", async () => {
    const t = await open();
    const y = t.yOf(S, "t-points");
    // The edge before 20 is the glyph's center: 3 px either side is it.
    expect([t.at(t.edge(20) - 3, y), t.at(t.edge(20) + 3, y)]).toEqual(["item 1", "item 1"]);
    expect(t.at(t.edge(N) - 3, y)).toBe("item 2");
  });

  it("ranges, points and variants: boundaries drawn on their edges, the end's included", async () => {
    const t = await open();
    // The ranges' three glyphs, the points' three and the variants' three, each on an edge; the ranges' at the block's
    // ends, on the data area's edges here, shifted in by its radius (5 px), so it is not cut.
    expect(t.xs.filter((x) => x === Math.round(t.edge(20))).length).toBe(3);
    for (const [i, inward] of [
      [0, 5],
      [N, -5],
    ]) {
      expect(t.xs.filter((x) => x === Math.round(t.edge(i))).length).toBe(2);
      expect(t.xs.filter((x) => x === Math.round(t.edge(i)) + inward).length).toBe(1);
    }
  });

  it("variants: a boundary hits where its glyph is drawn, at its place in the stack, the end's included", async () => {
    const t = await open();
    // Left of the edge before 20 is position 19's column: only the glyph's geometry finds the boundary there.
    expect(t.down(t.edge(20) - 2, "t-variants").has("item 1")).toBe(true);
    expect(t.down(t.edge(N) - 2, "t-variants").has("item 2")).toBe(true);
  });

  it("a reversed piece: on the right edge of the position after the boundary; the end's on the left", async () => {
    const t = await open([{ direction: -1 }]);
    const a = t.ctx.layout.axis(S);
    const right = (i: number) => Math.round(a?.x(0, ((a?.shownOf(0, i) ?? 0) + 1) as never) ?? 0);
    const left = (i: number) => Math.round(a?.x(0, a?.shownOf(0, i) ?? 0) ?? 0);
    expect(t.xs.filter((x) => x === right(20)).length).toBe(3);
    // The end's (the block's first edge on screen, the data area's): the ranges' glyph shifted in by its radius.
    expect(t.xs.filter((x) => x === left(N - 1)).length).toBe(2);
    expect(t.xs.filter((x) => x === left(N - 1) + 5).length).toBe(1);
  });

  it("two pieces of one unit: the junction's boundary drawn once, by the piece after it", async () => {
    const t = await open([{ range: [0, 20] }, { range: [20, N] }]);
    const a = t.ctx.layout.axis(S);
    const first = Math.round(a?.x(1, a?.shownOf(1, 20) ?? 0) ?? 0);
    const last = Math.round(a?.x(0, ((a?.shownOf(0, 19) ?? 0) + 1) as never) ?? 0);
    expect(first).not.toBe(last);
    expect([t.xs.filter((x) => x === first).length, t.xs.filter((x) => x === last).length]).toEqual([3, 0]);
  });

  it("a position the filter hides: its boundary after the position before it", async () => {
    const t = await open([{}], true);
    const a = t.ctx.layout.axis(S);
    expect(a?.shownOf(0, 20)).toBe(-1);
    const after19 = Math.round(a?.x(0, ((a?.shownOf(0, 19) ?? 0) + 1) as never) ?? 0);
    expect(t.xs.filter((x) => x === after19).length).toBe(3);
  });
});

describe("past a block's end", () => {
  it("the end boundary's glyph hits from the gap after its block, never the last residue", async () => {
    const b = createSeqViewSpec();
    const p = b.data.sequence({ ref: "p", value: "MKTAYIAKQRQISFVKSHFS" });
    const q = b.data.sequence({ ref: "q", value: "RQLEERLGLIEVQAPILSRV" });
    const end = b.data.annotation({
      ref: "end",
      geometry: "ranges",
      target: p,
      value: { boundary: [20], alt: ["GSG"] },
    });
    const sec = b.view({ ref: "v", label: "V" }).section({ ref: S }).block({ sequence: p }).block({ sequence: q });
    sec.track({ ref: "t-ranges", label: "Ranges", height: 2 }).ranges({ data: end });
    sec.track({ ref: "t-points", label: "Points", height: 2 }).points({ data: end });
    sec.track({ ref: "t-variants", label: "Variants", height: 2 }).variants({ data: end });
    const t = await lit(b.getDocument(), { width: 900, height: 300 });
    t.frame();
    const a = t.ctx.layout.axis(S);
    const edge = a ? a.x(0, (a.shownOf(0, 19) + 1) as never) : 0;
    expect(a?.x(1, 0 as never) ?? 0).toBeGreaterThan(edge + 3);
    const hit = Hit.create();
    const down = (x: number, track: string) => {
      const l = t.layout(S);
      const e = Array.from({ length: l.entryCount }, (_, k) => k).find((k) => l.entryTrack[k] === t.trackIndex(track));
      const top = l.scrolledTop + l.entryTop[e as number];
      const seen = new Set<string>();
      for (let y = top; y < top + l.entryHeight[e as number]; y += 0.5) {
        t.ctx.layout.hitTest(S, x, y, hit);
        seen.add(hit.kind === HoverKind.item ? "item" : hit.kind === HoverKind.seq ? "seq" : "none");
      }
      return seen;
    };
    for (const track of ["t-ranges", "t-points", "t-variants"]) {
      expect(down(edge + 2, track).has("item"), track).toBe(true);
      // Never the last residue from past it.
      expect(down(edge + 2, track).has("seq"), track).toBe(false);
    }
  });

  it("layers hit by position (letters of a mapped design, an alignment's rows, a matrix, pairs, a range's bar) never hover from past it", async () => {
    const P = "MKTAYIAKQRQISFVKSHFS";
    const b = createSeqViewSpec();
    const p = b.data.sequence({ ref: "p", value: P });
    const q = b.data.sequence({ ref: "q", value: "RQLEERLGLIEVQAPILSRV" });
    const d = b.data.sequence({ ref: "d", value: P });
    b.data.mapping({ ref: "d.on_p", from: p, to: d, method: "pairwise", gapped: [P, P] });
    const index = Array.from({ length: 20 }, (_, i) => i);
    const m = b.data.annotation({
      ref: "m",
      geometry: "matrix",
      target: p,
      states: ["x", "y"],
      value: { index, value: index.map(() => [1, 2]) },
    } as never);
    const pairs = b.data.annotation({
      ref: "pp",
      geometry: "pairs",
      target: p,
      value: { start: [19], b_target: ["p"], b_start: [2], count: [1] },
    } as never);
    const bar = b.data.annotation({ ref: "bar", geometry: "ranges", target: p, value: { start: [12], stop: [20] } });
    const sec = b.view({ ref: "v", label: "V" }).section({ ref: S }).block({ sequence: p }).block({ sequence: q });
    sec.track({ ref: "t-of", label: "Of", height: 2 }).letters({ of: d });
    sec.track({ ref: "t-matrix", label: "Matrix", height: 2 }).matrix({ data: m });
    sec.track({ ref: "t-pairs", label: "Pairs", height: 2 }).pairs({ data: pairs });
    sec.track({ ref: "t-bar", label: "Bar", height: 2 }).ranges({ data: bar });
    // An alignment's rows, in a view of their own: one block, its cells at their widest leaving room after it.
    const c = createSeqViewSpec();
    const aln = c.data.alignment({
      ref: "aln",
      source: { value: "# STOCKHOLM 1.0\na MKTAYIAKQR\nb MKT-YIAKQR\n//\n", format: "stockholm" },
    });
    c.view({ ref: "v", label: "V" })
      .section({ ref: S })
      .block({ alignment: aln })
      .track({ ref: "t-rows", label: "Rows" })
      .rows();
    const hit = Hit.create();
    const kinds = (t: Awaited<ReturnType<typeof lit>>, sec: string, at: number, track: string) => {
      const a = t.ctx.layout.axis(sec);
      const edge = a ? a.x(0, (a.shownOf(0, at) + 1) as never) : 0;
      const l = t.layout(sec);
      // In the gap before the next block, or (one block, its cells at their widest) in the data area after it.
      expect((a?.axis?.count ?? 0) > 1 ? (a?.x(1, 0 as never) ?? 0) : l.dataRight, track).toBeGreaterThan(edge + 3);
      const seen = new Set<number>();
      for (let e = 0; e < l.entryCount; e++) {
        if (l.entryTrack[e] !== t.trackIndex(track)) continue;
        const top = l.scrolledTop + l.entryTop[e];
        for (let y = top; y < top + l.entryHeight[e]; y += 0.5) {
          t.ctx.layout.hitTest(sec, edge + 2, y, hit);
          seen.add(hit.kind);
          // Inside the last cell, the same rows hover (the test reaches them).
          t.ctx.layout.hitTest(sec, edge - 1, y, hit);
          seen.add(-hit.kind - 1);
        }
      }
      return seen;
    };
    const t = await lit(b.getDocument(), { width: 900, height: 600 });
    t.frame();
    for (const track of ["t-of", "t-matrix", "t-pairs", "t-bar"]) {
      const seen = kinds(t, S, 19, track);
      expect(
        [...seen].filter((k) => k >= 0),
        track,
      ).toEqual([HoverKind.none]);
      expect(
        [...seen].some((k) => k < -1),
        track,
      ).toBe(true);
    }
    const ta = await lit(c.getDocument(), { width: 900, height: 300 });
    ta.frame();
    const seen = kinds(ta, S, 9, "t-rows");
    expect([...seen].filter((k) => k >= 0)).toEqual([HoverKind.none]);
    expect([...seen].some((k) => k < -1)).toBe(true);
  });
});

describe("an insert through a partial map", () => {
  const WT = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRV";
  /** The design: WT with GSGS after 20 and a tag of 3 after the last residue. */
  function mapped() {
    const b = createSeqViewSpec();
    const wt = b.data.sequence({ ref: "wt", value: WT });
    const d = b.data.sequence({ ref: "d", value: `${WT.slice(0, 20)}GSGS${WT.slice(20)}HHH` });
    b.data.mapping({
      ref: "d.on_wt",
      from: wt,
      to: d,
      method: "pairwise",
      gapped: [`${WT.slice(0, 20)}----${WT.slice(20)}---`, `${WT.slice(0, 20)}GSGS${WT.slice(20)}HHH`],
    });
    const own = b.data.annotation({
      ref: "own",
      geometry: "ranges",
      target: d,
      value: { start: [21, 44], stop: [23, 47], label: ["inside the graft", "the tag"] },
    });
    b.view({ ref: "v", label: "V" })
      .section({ ref: S })
      .block({ sequence: wt })
      .track({ ref: "t-own", label: "Own", height: 2 })
      .ranges({ data: own });
    return b.getDocument();
  }

  it("lands on the boundary after it on the other side; beyond the pairs (a tag) draws nothing; counted", async () => {
    const t = await lit(mapped(), { width: 900, height: 300 });
    const be = t.frame(recording());
    const a = t.ctx.layout.axis(S);
    const edge20 = Math.round(a?.x(0, a?.shownOf(0, 20) ?? 0) ?? 0);
    expect(glyphXs(be)).toEqual([edge20]);
    const insert = t.ctx.issues.current.value.filter((i) => i.code === "svs.layer.insert");
    expect(insert.map((i) => i.message.split(" ")[0])).toEqual(["1"]);
  });
});
