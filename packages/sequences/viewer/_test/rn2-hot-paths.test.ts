/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// RN2's hot paths on the recording backend (K4): density composed per device column (stacked in palette order by share,
// majority at its share, one band and no item rects where blocks mix pyramid and item reductions, a width sweep
// flipping the mode once each way); the interval tree against a brute force, with visit counts that stay
// O((k + 1) log n) beside one long pair; the pairs' draw, hit and `tiny` from the end index and the tree, arcs grouped
// by color; the matrix as tiles of at most 4,096 samples, one run on an unfiltered block, nothing rebuilt or scanned on
// a second frame.

import { OrderedSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import {
  Arcs,
  densityColumns,
  Hit,
  keyOf,
  litPairs,
  MarkBit,
  type PairsState,
  PathBatch,
  pairsAt,
  pairsInView,
  pushArc,
} from "../src/draw/api.ts";
import { Dirty, SectionComposer } from "../src/draw/index.ts";
import { MATRIX_TILE, MatrixVisits } from "../src/draw/matrix-tiles.ts";
import { KernelRenderers, type LayerRenderer } from "../src/index.ts";
import { IndexVisits, overlapping, intervalTree as tree } from "../src/intervals.ts";
import { MarkSet } from "../src/marks.ts";
import { Op, RecordingBackend } from "../src/testing/index.ts";
import { entryOf, layerOf, record } from "./draw-fixtures.ts";
import { opened } from "./layout-fixtures.ts";
import {
  bandColumns,
  densityDoc,
  densityItems,
  forced,
  lastDraw,
  pairsDoc,
  rng,
  wideMatrix,
  withDensity,
} from "./rn2-fixtures.ts";

/** The palette's colors of the track's categories, in legend order. */
function legend(ctx: Awaited<ReturnType<typeof opened>>["ctx"], label: string): number[] {
  const { k } = layerOf(ctx, label);
  const m = ctx.graph.assembled.value?.view;
  const out: number[] = [];
  for (const t of m?.prepared[k].themes ?? [])
    t?.legend?.({ item: (_l, c) => out.push(c), ramp: () => {}, size: () => {} });
  return out;
}

describe("density composed per device column", () => {
  /** Per device column of the band, the items' counts by color (k0, k1), brute force: each at its center's column. */
  function brute(
    ctx: Awaited<ReturnType<typeof opened>>["ctx"],
    frame: ReturnType<typeof record>["frame"],
    o: Parameters<typeof densityItems>[0] = {},
  ) {
    const v = frame.v;
    const { k } = layerOf(ctx, "D");
    const f = ctx.graph.assembled.value?.view?.prepared[k].state as {
      count: number;
      block: Int32Array;
      start: Int32Array;
      stop: Int32Array;
      row: Int32Array;
    };
    const counts = new Map<number, [number, number]>();
    const category = densityItems(o).category;
    for (let i = 0; i < f.count; i++) {
      const c = (f.start[i] + Math.max(f.stop[i], f.start[i] + 1) - 1) >> 1;
      const q = v.shownOf(f.block[i], c);
      if (q < v.visLo[f.block[i]] || q >= v.visHi[f.block[i]]) continue;
      const x = v.xd(f.block[i], q);
      if (x < Math.round(v.left * v.dpr) || x >= Math.round(v.right * v.dpr)) continue;
      const e = counts.get(x) ?? [0, 0];
      e[category[f.row[i]] === "k1" ? 1 : 0]++;
      counts.set(x, e);
    }
    return counts;
  }

  it("stacked: a column of two colors shows both, the palette's first at the bottom, heights by share", async () => {
    const t = await opened(densityDoc(), { renderers: withDensity("points") });
    const [k0, k1] = legend(t.ctx, "D");
    const { rec, frame } = record(t.ctx, "D");
    expect(lastDraw.mode).toBe(1);
    const e = entryOf(t.ctx, "D");
    const y1 = Math.round((e.top + e.height) * e.dpr);
    const left = Math.round(e.left * e.dpr);
    const cols = bandColumns(rec, left, Math.round(e.right * e.dpr));
    const counts = brute(t.ctx, frame);
    let both = 0;
    for (const [x, [n0, n1]] of counts) {
      const c = cols[x - left];
      const total = n0 + n1;
      const bar = c.reduce((s, r) => s + r[1], 0);
      expect(c[0][0] + c[0][1], `column ${x} from the bottom`).toBe(y1);
      if (n0 > 0 && n1 > 0) {
        both++;
        expect(c.map((r) => r[2])).toEqual([k0, k1]);
        expect(Math.abs(c[0][1] - (bar * n0) / total)).toBeLessThanOrEqual(1);
      } else expect(c.map((r) => r[2])).toEqual([n0 > 0 ? k0 : k1]);
    }
    expect(both).toBeGreaterThan(50);
  });

  it("majority: the plurality color alone, at its share of the column's bar", async () => {
    const t = await opened(densityDoc({ density: "majority" }), { renderers: withDensity("points") });
    const [k0, k1] = legend(t.ctx, "D");
    const { rec, frame } = record(t.ctx, "D");
    const e = entryOf(t.ctx, "D");
    const left = Math.round(e.left * e.dpr);
    const cols = bandColumns(rec, left, Math.round(e.right * e.dpr));
    const counts = brute(t.ctx, frame);
    const max = Math.max(...[...counts.values()].map(([a, b]) => a + b));
    const h = Math.round((e.top + e.height) * e.dpr) - Math.round(e.top * e.dpr);
    for (const [x, [n0, n1]] of counts) {
      const c = cols[x - left];
      expect(c).toHaveLength(1);
      expect(c[0][2]).toBe(n1 > n0 ? k1 : k0);
      const bar = Math.min(h, Math.max(2, Math.round(h * Math.sqrt((n0 + n1) / max))));
      expect(Math.abs(c[0][1] - (bar * Math.max(n0, n1)) / (n0 + n1))).toBeLessThanOrEqual(1);
    }
  });

  it("coverage from the count pyramid: every column drawn once, stacked, bars from the bottom", async () => {
    const t = await opened(densityDoc({ n: 400_000, type: "ranges", items: 20_000, len: 3000 }), {
      renderers: withDensity("ranges"),
    });
    const { rec } = record(t.ctx, "D");
    expect(lastDraw.mode).toBe(1);
    const e = entryOf(t.ctx, "D");
    const left = Math.round(e.left * e.dpr);
    const cols = bandColumns(rec, left, Math.round(e.right * e.dpr));
    const y1 = Math.round((e.top + e.height) * e.dpr);
    expect(cols.filter((c) => c.length > 0).length).toBeGreaterThan(0.9 * cols.length);
    for (const c of cols) if (c.length > 0) expect(c[0][0] + c[0][1]).toBe(y1);
  });

  it("a layer whose blocks mix pyramid and item reductions draws one band and no item rects", async () => {
    const t = await opened(densityDoc({ n: 200_000, blocks: 2, items: 20_000 }), { renderers: withDensity("points") });
    const { section } = layerOf(t.ctx, "D");
    const ax = t.ctx.layout.axis(section)?.axis;
    if (!ax) throw new Error("no axis");
    // Block 1 read through its items (a filter of all its positions): the same columns, counted item by item.
    const all = OrderedSet.ofRange(ax.base[1], ax.base[1] + ax.length[1]);
    const whole = record(t.ctx, "D").rec;
    forced.positions = [null, all];
    try {
      const { rec } = record(t.ctx, "D");
      const e = entryOf(t.ctx, "D");
      const y1 = Math.round((e.top + e.height) * e.dpr);
      const b = rec.rectItems;
      expect(b.n).toBeGreaterThan(100);
      // One band: every rect a device column wide, each column's rects stacked without gaps from the band's bottom;
      // no item rects.
      for (let i = 0; i < b.n; i++) expect(b.w[i]).toBe(1);
      const cols = bandColumns(rec, 0, Math.round(e.right * e.dpr) + 1);
      const tops = new Map<number, number>();
      cols.forEach((c, x) => {
        let y = y1;
        for (const [ry, rh] of c) {
          expect(ry + rh, `column ${x}`).toBe(y);
          y = ry;
        }
        if (c.length > 0) tops.set(x, y);
      });
      // Both blocks drawn, the item-counted one with the same columns as through its pyramid, give or take a bin.
      const x1 = t.ctx.layout.axis(section)?.blockX0[1] ?? 0;
      const right = [...tops.keys()].filter((x) => x > x1 * e.dpr + 2);
      const before = new Set<number>();
      for (let i = 0; i < whole.rectItems.n; i++)
        if (whole.rectItems.x[i] > x1 * e.dpr + 2) before.add(whole.rectItems.x[i]);
      expect(right.length).toBeGreaterThan(0.9 * before.size);
    } finally {
      forced.positions = null;
    }
  });
});

describe("one mode per layer per frame, with hysteresis", () => {
  it("a jittered width sweep across the threshold flips items and density once each way", async () => {
    const t = await opened(densityDoc({ n: 20_000, type: "ranges", items: 600, len: 20 }), {
      renderers: withDensity("ranges"),
    });
    const { section } = layerOf(t.ctx, "D");
    const vp = t.ctx.viewports.get(section);
    if (!vp) throw new Error("no viewport");
    const sweep: number[] = [];
    for (let w = 400; w <= 1200; w += 6) sweep.push(w, w - 4);
    // One composer's view throughout, as a surface keeps its pass's `DrawView` (the mode is kept per view and layer).
    const frame = record(t.ctx, "D").frame;
    const rec = new RecordingBackend();
    const run = async (widths: number[]) => {
      const modes: number[] = [];
      const naive: number[] = [];
      for (const w of widths) {
        vp.setSize(w, 400, 1);
        vp.fit();
        await t.ctx.graph.update();
        frame.draw(rec);
        modes.push(lastDraw.mode);
        naive.push(lastDraw.count > lastDraw.columns ? 1 : 0);
      }
      const flips = (m: number[]) => m.slice(1).filter((x, i) => x !== m[i]).length;
      return { flips: flips(modes), naive: flips(naive), first: modes[0], last: modes[modes.length - 1] };
    };
    const up = await run(sweep);
    expect(up).toMatchObject({ flips: 1, first: 1, last: 0 });
    // Without the hysteresis the jitter flips the mode back and forth.
    expect(up.naive).toBeGreaterThan(1);
    const down = await run([...sweep].reverse());
    expect(down).toMatchObject({ flips: 1, first: 0, last: 1 });
    expect(densityColumns).toBeTypeOf("function");
  });
});

describe("the hysteresis through the composer (CS1a: one `DrawView` per composer, kept across frames)", () => {
  it("a jittered resize sweep flips the mode once each way, frame after frame of the section", async () => {
    const t = await opened(densityDoc({ n: 20_000, type: "ranges", items: 600, len: 20 }), {
      renderers: withDensity("ranges"),
    });
    const { section } = layerOf(t.ctx, "D");
    const vp = t.ctx.viewports.get(section);
    if (!vp) throw new Error("no viewport");
    const composer = new SectionComposer(t.ctx, section);
    const rec = new RecordingBackend();
    const sweep: number[] = [];
    for (let w = 400; w <= 1200; w += 6) sweep.push(w, w - 4);
    const run = async (widths: number[]) => {
      const modes: number[] = [];
      for (const w of widths) {
        vp.setSize(w, 400, 1);
        vp.fit();
        await t.ctx.graph.update();
        rec.reset();
        lastDraw.mode = -1;
        composer.frame(rec, Dirty.tracks);
        // A frame whose size did not change draws nothing.
        if (lastDraw.mode >= 0) modes.push(lastDraw.mode);
      }
      expect(modes.length).toBeGreaterThan(widths.length - 2);
      return { flips: modes.slice(1).filter((x, i) => x !== modes[i]).length, first: modes[0], last: modes.at(-1) };
    };
    expect(await run(sweep)).toEqual({ flips: 1, first: 1, last: 0 });
    expect(await run([...sweep].reverse())).toEqual({ flips: 1, first: 0, last: 1 });
    composer.dispose();
  });
});

describe("the interval tree and the pairs' indexes (K4's counters)", () => {
  it("overlap queries equal a brute force, in sorted order, beside a long interval", () => {
    const rnd = rng(3);
    const n = 3000;
    const iv = Array.from({ length: n }, () => {
      const a = Math.floor(rnd() * 1e6);
      return [a, a + Math.floor(rnd() * 500)];
    });
    iv.push([0, 2e6]);
    iv.sort((x, y) => x[0] - y[0]);
    const lo = Float64Array.from(iv, (x) => x[0]);
    const t = tree(Float64Array.from(iv, (x) => x[1]));
    const out = new Int32Array(iv.length);
    const q = new Float64Array(2);
    for (let r = 0; r < 300; r++) {
      q[0] = Math.floor(rnd() * 1.1e6);
      q[1] = q[0] + (r % 3 === 0 ? 0 : Math.floor(rnd() * 3000));
      const n0 = IndexVisits.nodes;
      const k = overlapping(lo, t, q, out);
      const want = iv.flatMap((x, i) => (x[0] <= q[1] && x[1] >= q[0] ? [i] : []));
      expect(Array.from(out.subarray(0, k))).toEqual(want);
      // Output-sensitive: about 2 log n nodes per interval reported (a prefix maximum would visit every earlier one).
      expect(IndexVisits.nodes - n0).toBeLessThanOrEqual(4 * 12 * (k + 1));
    }
  });

  it("the arc hit visits the arcs over the point, never every earlier pair (one long pair spans them all)", async () => {
    const t = await opened(pairsDoc(), { renderers: [...KernelRenderers] as LayerRenderer[] });
    const { k, section } = layerOf(t.ctx, "Contacts");
    const st = t.ctx.graph.assembled.value?.view?.prepared[k].state as PairsState;
    t.ctx.viewports.get(section)?.setView(4, 4 * 15_000, "api");
    await t.ctx.graph.update();
    const e = entryOf(t.ctx, "Contacts");
    const h = Hit.create();
    // The first hit after the view moved finds the widest arc in view (`Arcs.frame`): the pairs in view, once.
    t.ctx.layout.hitTest(section, e.left + 1, e.top + e.height / 2, h);
    let worst = 0;
    for (let x = e.left + 5; x < e.right; x += 17) {
      const v0 = IndexVisits.intervals;
      t.ctx.layout.hitTest(section, x, e.top + e.height / 2, h);
      worst = Math.max(worst, IndexVisits.intervals - v0);
    }
    // The long pair and the few short ones over a position; the old prefix maximum walked about 3,700 pairs here.
    expect(worst).toBeGreaterThan(0);
    expect(worst).toBeLessThan(40);
    // A stab equals the brute force.
    const out = new Int32Array(st.count);
    const n = pairsAt(st, 0, 15_100, out);
    const key = keyOf(0, 15_100);
    const want = [];
    for (let p = 0; p < st.count; p++) if (st.lo[st.order.indexOf(p)] <= key && st.hi[p] >= key) want.push(p);
    expect(Array.from(out.subarray(0, n)).sort((a, b) => a - b)).toEqual(want.sort((a, b) => a - b));
  });

  it("the arc hit picks the arc drawn on top: the overlay's last lit one, else the draw's last over the point", async () => {
    const { lit, overlayOf } = await import("./rn2b-fixtures.ts");
    const t = await opened(pairsDoc({ count: 12_000, colorBy: "category", long: false }), {
      renderers: [...KernelRenderers] as LayerRenderer[],
    });
    const { k, section } = layerOf(t.ctx, "Contacts");
    const st = t.ctx.graph.assembled.value?.view?.prepared[k].state as PairsState;
    t.ctx.viewports.get(section)?.setView(3, 3 * 15_000, "api");
    await t.ctx.graph.update();
    const e = entryOf(t.ctx, "Contacts");
    const { rec, frame } = record(t.ctx, "Contacts");
    const v = frame.v;
    // Each recorded arc's pieces (an arc beyond the window drawn as its pieces in view) and its pair, by its ends
    // (pairs sharing both ends, and pieces of an arc cut by the window or the band, left out: -1), in the order drawn:
    // the draw's, then the overlay's over them (a few lit ends).
    const byEnds = new Map<string, number>();
    for (let p = 0; p < st.count; p++) {
      if (!Arcs.ends(v, st, p)) continue;
      const key = `${Math.round(st.xs[0])}>${Math.round(st.xs[1])}`;
      byEnds.set(key, byEnds.has(key) ? -1 : p);
    }
    const order: { k: number; x0: number; x1: number; y0: number; cy: number; y1: number }[] = [];
    const arcsOf = (r: RecordingBackend) => {
      const ps = r.pathItems;
      for (let p = 0; p < ps.n; p++) {
        const c = ps.coords.subarray(ps.cstart[p], p + 1 < ps.n ? ps.cstart[p + 1] : ps.nc);
        if (ps.mode[p] !== 1 || ps.verbs[ps.start[p] + 1] !== 2) continue;
        const whole = c.length === 6 && c[1] === c[5];
        for (let o = 0; o < c.length; o += 6)
          order.push({
            k: whole ? (byEnds.get(`${Math.round(c[0])}>${Math.round(c[4])}`) ?? -1) : -1,
            x0: c[o],
            x1: c[o + 4],
            y0: c[o + 1],
            cy: c[o + 3],
            y1: c[o + 5],
          });
      }
    };
    arcsOf(rec);
    const drawnArcs = order.length;
    const runs: [number, number][] = [];
    for (let j = 0; j < drawnArcs; j += Math.ceil(drawnArcs / 6))
      if (order[j].k >= 0) runs.push([st.start[2 * order[j].k], st.stop[2 * order[j].k]]);
    arcsOf(overlayOf(t.ctx, "Contacts", lit(MarkBit.selection, runs)));
    expect(order.length).toBeGreaterThan(drawnArcs);
    const h = Hit.create();
    let crossed = 0;
    let checked = 0;
    for (let x = e.left + 3.5; x < e.right; x += 3)
      for (let y = e.top + 6.5; y < e.top + e.height - 6; y += 2) {
        // The brute force: the last drawn curve within 4 px of the point (the quad's y at the point's x, linear in t).
        let top = -1;
        let under = 0;
        let largest = -1;
        for (let j = 0; j < order.length; j++) {
          const a = order[j];
          const tt = (x - a.x0) / (a.x1 - a.x0);
          const at = (1 - tt) * (1 - tt) * a.y0 + 2 * tt * (1 - tt) * a.cy + tt * tt * a.y1;
          if (!(tt >= 0 && tt <= 1) || Math.abs(y - at) > 4) continue;
          under++;
          top = j;
          largest = Math.max(largest, a.k);
        }
        // Within 5 px of the ends' line the end under the point wins.
        if (top < 0 || order[top].k < 0 || Math.abs(y - st.frame[0]) <= 5) continue;
        t.ctx.layout.hitTest(section, x, y, h);
        const want = st.set.annotations[st.annotation[order[top].k]]?.item[st.row[order[top].k]];
        expect(h.item, `at ${x}, ${y}`).toBe(want);
        checked++;
        if (under > 1 && largest !== order[top].k) crossed++;
      }
    // Points under two arcs or more decide by z-order: the largest pair index (the old pick) differs at some.
    expect(checked).toBeGreaterThan(200);
    expect(crossed).toBeGreaterThan(20);
  });

  it("the draw strokes the pairs in view once each (those whose curve shows), grouped by color, visiting only them", async () => {
    const t = await opened(pairsDoc(), { renderers: [...KernelRenderers] as LayerRenderer[] });
    const { k, section } = layerOf(t.ctx, "Contacts");
    const st = t.ctx.graph.assembled.value?.view?.prepared[k].state as PairsState;
    t.ctx.viewports.get(section)?.setView(4, 4 * 15_000, "api");
    await t.ctx.graph.update();
    const v0 = IndexVisits.intervals;
    const { rec, frame } = record(t.ctx, "Contacts");
    const visited = IndexVisits.intervals - v0;
    // In view, brute force: a pair whose keys meet the view's.
    const v = frame.v;
    const [a, b] = [v.indexOf(0, v.visLo[0] as never), v.indexOf(0, (v.visHi[0] - 1) as never)];
    let inView = 0;
    for (let p = 0; p < st.count; p++)
      if (st.hi[p] >= keyOf(0, a) && st.lo[st.order.indexOf(p)] <= keyOf(0, b)) inView++;
    expect(inView).toBeLessThan(st.count / 5);
    const scratch = new Int32Array(st.count);
    expect(pairsInView(v, st, scratch)).toBe(inView);
    // The frame's query (the pairs in view) and `Arcs.frame`'s widest (once per view: the pairs over the window and its
    // margins, half a window more): never a scan of every pair.
    expect(visited).toBeLessThanOrEqual(3 * inView);
    // Arcs: one subpath each, but a pair whose curve passes above the band throughout the window (its ends beyond
    // either edge); each color one contiguous run of subpaths.
    const p = rec.pathItems;
    let strokes = 0;
    const runs: number[] = [];
    for (let i = 0; i < p.n; i++) {
      if (p.mode[i] !== 1) continue;
      strokes++;
      if (runs[runs.length - 1] !== p.color[i]) runs.push(p.color[i]);
    }
    let shows = 0;
    for (let q = 0; q < st.count; q++)
      if (st.hi[q] >= keyOf(0, a) && st.lo[st.order.indexOf(q)] <= keyOf(0, b))
        shows += pushArc(v, st, q, 2, 0, new PathBatch(), new PathBatch()) ? 1 : 0;
    expect(shows).toBeGreaterThan(inView - 5);
    expect(strokes).toBe(shows);
    expect(new Set(runs).size).toBe(runs.length);
    expect(rec.count(Op.paths)).toBeGreaterThan(0);
    void Arcs;
  });

  it("the overlay's lit pairs: through the channels' runs and the end index, each once", async () => {
    const t = await opened(pairsDoc({ count: 2000 }), { renderers: [...KernelRenderers] as LayerRenderer[] });
    const { k, section } = layerOf(t.ctx, "Contacts");
    const st = t.ctx.graph.assembled.value?.view?.prepared[k].state as PairsState;
    const v = t.ctx.layout.axis(section);
    if (!v) throw new Error("no axis");
    const sets = [new MarkSet(), new MarkSet()];
    const lit: [number, number][][] = [
      [
        [100, 160],
        [5000, 5003],
      ],
      [[150, 400]],
    ];
    sets.forEach((m, j) => {
      m.reset(1);
      m.begin(0);
      for (const [a, b] of lit[j]) m.add(a, b);
      m.close(0, 20_000);
    });
    const marks = {
      bitsAt: () => 0,
      any: () => false,
      runs: (bit: number) => sets[bit === MarkBit.selection ? 0 : 1],
      primary: () => -1,
      item: () => -1,
      itemAnnotation: () => -1,
      element: () => false,
    };
    const out = new Int32Array(st.count);
    const e0 = IndexVisits.ends;
    const n = litPairs(v, st, marks, MarkBit.selection | MarkBit.highlight, out);
    const on = (i0: number, i1: number) => lit.flat().some(([a, b]) => i0 < b && i1 > a);
    const want = [];
    for (let p = 0; p < st.count; p++)
      if (on(st.start[2 * p], st.stop[2 * p]) || on(st.start[2 * p + 1], st.stop[2 * p + 1])) want.push(p);
    expect(Array.from(out.subarray(0, n)).sort((a, b) => a - b)).toEqual(want);
    expect(want.length).toBeGreaterThan(5);
    // The ends in the lit runs (and up to `maxEnd` before them), not every end.
    expect(IndexVisits.ends - e0).toBeLessThan(4 * want.length + 20);
  });

  it("`tiny` below the row's minimum: the ends in view as density from the end index", async () => {
    const t = await opened(pairsDoc({ count: 400 }), { renderers: [...KernelRenderers] as LayerRenderer[] });
    const { k } = layerOf(t.ctx, "Contacts");
    const entry = KernelRenderers.find((r) => r.name === "pairs") as LayerRenderer;
    const st = t.ctx.graph.assembled.value?.view?.prepared[k].state as PairsState;
    const v0 = IndexVisits.nodes;
    const { rec, frame } = record(t.ctx, "Contacts", new RecordingBackend(), true);
    expect(frame.v.height).toBe((entry.minHeight ?? 6) - 1);
    // Every end in view counted once: the bars' total area is at least the columns that hold an end.
    const b = rec.rectItems;
    expect(b.n).toBeGreaterThan(50);
    const cols = new Set<number>();
    for (let i = 0; i < b.n; i++) cols.add(b.x[i]);
    const ends = new Set<number>();
    const v = frame.v;
    for (let e = 0; e < 2 * st.count; e++) {
      const q = v.shownOf(st.block[e], st.start[e]);
      if (q >= v.visLo[0] && q < v.visHi[0])
        ends.add(Math.round(((v.x(0, q as never) + v.x(0, (q + 1) as never)) / 2) * v.dpr));
    }
    for (const x of ends) expect(cols.has(x) || cols.has(x - 1) || cols.has(x + 1)).toBe(true);
    // No arc or glyph below the minimum; the tree not walked for ends.
    expect([rec.pathItems.n, rec.glyphItems.n]).toEqual([0, 0]);
    expect(IndexVisits.nodes - v0).toBeLessThan(st.count);
    // At the minimum: the draw, glyphs and arcs.
    const full = record(t.ctx, "Contacts").rec;
    expect(full.pathItems.n).toBeGreaterThan(0);
  });
});

describe("matrix tiles (K4)", () => {
  it("40,000 positions: tiles of at most 4,096 samples, one run, nothing rebuilt or scanned on the next frame", async () => {
    const t = await opened(wideMatrix(), { renderers: [...KernelRenderers] as LayerRenderer[] });
    const { section } = layerOf(t.ctx, "M");
    const e = entryOf(t.ctx, "M");
    for (const cell of [1, 0.4, 0.02]) {
      t.ctx.viewports.get(section)?.setView(cell, cell * 1000, "api");
      await t.ctx.graph.update();
      const m0 = { ...MatrixVisits };
      const { rec } = record(t.ctx, "M");
      expect(MatrixVisits.runs - m0.runs).toBe(1);
      expect(MatrixVisits.scans).toBe(0);
      const images = rec.count(Op.image);
      expect(images).toBe(MatrixVisits.images - m0.images);
      // Every image a tile's part, the band's rows the bands `matrixBandAt` hits.
      for (let k = 0; k < rec.opCount; k++) {
        if (rec.ops[4 * k] !== Op.image) continue;
        const a = rec.ops[4 * k + 3];
        const [sw, sh, dy, dh] = [rec.args[a + 2], rec.args[a + 3], rec.args[a + 5], rec.args[a + 7]];
        expect(sw).toBeLessThanOrEqual(MATRIX_TILE);
        expect([sh, dy, dh]).toEqual([
          3,
          Math.round(e.top * e.dpr),
          Math.round((e.top + e.height) * e.dpr) - Math.round(e.top * e.dpr),
        ]);
      }
      const built = MatrixVisits.tiles;
      record(t.ctx, "M", rec);
      expect(MatrixVisits.tiles).toBe(built);
    }
  });
});
