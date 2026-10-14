/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for the axis and the viewport at pointer rate (VC3; seqspec:viewer#performance), on `heapScaling` (a step
// is one move): per-move axis lookups (`indexOf`, `shownOf`, `breakAt`, the view's `x`) and both run converters on a
// forward ranged unit, a filtered, a reversed and a cyclic piece; viewport pans and anchored zooms with a live
// subscriber; programmatic (`api`) moves, each committed at once, with no `committed` subscriber (a subscriber's
// records are its own, and show). The positive control runs the axis lookups' step allocating an object each
// (`leaking`) through the same measurement, which must catch it.

import { beforeAll, describe, expect, it } from "vitest";
import { type Axis, AxisGeometry, SectionAxisView, type ShownIndex } from "../src/axis.js";
import { axisOf, open, rangedDoc, resolved, setPiece, vp } from "./axis-fixtures.js";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from "./heap.js";

const N = 50_000;
let axes: Axis[] = [];
let sink = 0;
const out = new Int32Array(16);

beforeAll(async () => {
  const [v, f, w, c] = await Promise.all([
    resolved(rangedDoc(), "v"),
    resolved(rangedDoc(), "f"),
    resolved(rangedDoc(), "w"),
    resolved(rangedDoc(), "w"),
  ]);
  setPiece(w.view, 0, 0, 20, -1, 20);
  setPiece(c.view, 0, 15, 5, 1, 20);
  axes = [axisOf(v), axisOf(f), axisOf(w), axisOf(c)];
}, 60_000);

describe("the axis allocates nothing per move", () => {
  it("indexOf, shownOf, breakAt and x over forward, filtered, reversed and cyclic pieces", () => {
    const views = axes.map((ax) => {
      const v = new SectionAxisView();
      const l = { geometry: AxisGeometry.linear, dataLeft: 10, dataRight: 300, scrolledTop: 0, height: 100 };
      const s = { cell: 9, scrollX: 3, scrollY: 0, dpr: 2, lod: 0, contentHeight: 0 };
      v.sync(ax, l, { ...s, visLo: new Int32Array(1), visHi: new Int32Array(1) } as never);
      return v;
    });
    const step = (i: number) => {
      for (let a = 0; a < views.length; a++) {
        const v = views[a];
        const k = (i % 12) as ShownIndex;
        sink = (sink + v.indexOf(0, k) + v.shownOf(0, i % 21) + (v.breakAt(0, k) ? 1 : 0) + v.x(0, k)) | 0;
      }
    };
    const r = heapScaling(N, step);
    expect(sink).not.toBe(0.5);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
  });

  it("runToShown and shownToRuns", () => {
    const r = heapScaling(N, (i) => {
      for (let a = 0; a < axes.length; a++) {
        const ax = axes[a];
        sink = (sink + ax.runToShown(0, i % 20, (i % 20) + 4, out) + ax.shownToRuns(0, i % 8, (i % 8) + 3, out)) | 0;
      }
    });
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });
});

describe("the viewport allocates nothing per move", () => {
  it("pans and anchored zooms at pointer rate, with a live subscriber", async () => {
    const { ctx } = await open();
    ctx.views.setCurrent("g");
    await ctx.graph.update();
    const v = vp(ctx, "g.one");
    v.zoomTo(0, 8, 12); // cell 28 (the most), centered: every pan and zoom below changes the view
    v.live.version.subscribe(() => {
      sink = (sink + v.live.current.scrollX + v.live.current.visLo[0]) & 0xffff;
    });
    const before = v.live.version.value;
    const r = heapScaling(N, (k) => {
      // Constants per branch: a double chosen by `?:` is boxed by the caller before the call.
      if (!(k & 64)) v.panBy(k & 1 ? 3 : -3, 0, "pointer");
      else if (k & 1) v.zoomAt(1.25, 40, "pointer");
      else v.zoomAt(0.8, 40, "pointer");
    });
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(v.live.version.value - before).toBeGreaterThan(N);
    expect(v.live.current.interacting).toBe(true); // one settle timer for the whole gesture
  });
});

describe("an api move allocates nothing while nobody reads `committed`", () => {
  async function api() {
    const { ctx } = await open();
    ctx.views.setCurrent("g");
    await ctx.graph.update();
    const v = vp(ctx, "g.one");
    v.zoomTo(0, 8, 12);
    return v;
  }
  it("pans and views set with origin api: committed each, built only when read", async () => {
    const v = await api();
    const before = v.committed.value;
    const r = heapScaling(N, (k) => {
      if (k & 1) v.panBy(k & 2 ? 3 : -3, 0, "api");
      else v.setView(k & 2 ? 20 : 24, 40, "api");
    });
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    // Read, the last move's record; still at once for api (no interaction).
    expect(v.committed.value).not.toBe(before);
    expect([v.committed.value.origin, v.committed.value.scrollX, v.live.current.interacting]).toEqual([
      "api",
      v.live.current.scrollX,
      false,
    ]);
  });
  it("with a subscriber, each move's record is delivered and shows", async () => {
    const v = await api();
    let n = 0;
    v.committed.subscribe(() => n++);
    // Short runs: the records fill the young generation within a few thousand moves.
    const M = 1_000;
    const r = heapScaling(M, (k) => v.panBy(k & 1 ? 3 : -3, 0, "api"));
    expect(n).toBeGreaterThan(M);
    expect(isFlat(r), describeScaling(r, M)).toBe(false);
  });
});
