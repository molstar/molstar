/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The navigator's per-frame path (RN3 hot paths; seqspec:viewer#extension-points; K5), on the recording backend: the
// map as at the fit (per-block extents equal across zooms, the window mapped through the blocks) and the gestures
// through its inverse; the stacked navigator's right edge and marks per lane; the channels' runs as marks, accent
// included; a navigator on the ruler draws only its hairline; horizontal and vertical grips; the keyed map image built
// once.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { AxisGeometry, navUnits, SectionAxisView } from "../src/axis.js";
import { KeyedImage, SkinToken } from "../src/draw/api.js";
import { drawGrips, drawNavigator, drawNavigatorWindow, gripAt, NavView } from "../src/draw/nav.js";
import { DefaultNavigators, MarkBit, navigatorFrom } from "../src/index.js";
import { SectionMarks } from "../src/marks.js";
import { fixtureMeasurer, Op, RecordingBackend } from "../src/testing/index.js";
import { TextManager } from "../src/text.js";
import { WindowKind } from "../src/windows.js";
import { axisOf, blocksDoc, resolved } from "./axis-fixtures.js";
import { driven } from "./input-fixtures.js";

const L = { geometry: AxisGeometry.linear, dataLeft: 40, dataRight: 440, scrolledTop: 20, height: 120 };

/** A section's axis view over [40, 440), synced at a cell and scroll, its fit cell as the viewport computes it. */
async function synced(doc = blocksDoc(3, 50)) {
  const r = await resolved(doc, "v");
  const ax = axisOf(r);
  const v = new SectionAxisView();
  const gaps = ax.stacked ? 0 : (ax.count - 1) * ax.gap;
  const minCell = (400 - gaps) / ax.cells;
  const at = (cell: number, scrollX: number) => {
    const n = ax.count;
    const s = { cell, minCell, scrollX, scrollY: 0, dpr: 1, lod: 0, contentHeight: 0 };
    v.sync(ax, L, { ...s, visLo: new Int32Array(n), visHi: new Int32Array(n) } as never);
    return v;
  };
  return { ax, v, at, minCell };
}

function stackedDoc() {
  const b = createSeqViewSpec();
  const v = b.view({ ref: "v", label: "V", arrangement: "stacked" });
  for (const [i, len] of [80, 120, 40].entries())
    v.block({ sequence: b.data.sequence({ ref: `c${i}`, value: "M".repeat(len) }) });
  v.track({ ref: "v.l" }).letters();
  return b.getDocument();
}

/** A navigator view over the axis view, its row [0, 20), the window as the axis view maps it. */
function navView(v: SectionAxisView, parts = DefaultNavigators[2].parts, dpr = 1) {
  const n = new NavView().from(v);
  n.backend = new RecordingBackend();
  n.text = new TextManager(fixtureMeasurer(), () => 64);
  n.parts = parts;
  n.dpr = dpr;
  n.top = 0;
  n.height = 20;
  n.x0 = v.navX0;
  n.x1 = v.navX1;
  n.gripTop = 0;
  n.gripBottom = 20;
  return n;
}
const rec = (n: NavView) => n.backend as RecordingBackend;
/** The rects drawn in a color: [x, y, w, h] (device px). */
function rects(b: RecordingBackend, color: number) {
  const r = b.rectItems;
  const out: number[][] = [];
  for (let i = 0; i < r.n; i++) if (r.color[i] === color) out.push([r.x[i], r.y[i], r.w[i], r.h[i]]);
  return out;
}

describe("the map as at the fit", () => {
  it("per-block extents and the map's width are equal across zooms; the window maps through the blocks", async () => {
    const { ax, v, at, minCell } = await synced();
    at(minCell, 0);
    const x0s = Array.from(v.navX0s);
    const x1s = Array.from(v.navX1s);
    // At the fit the map is the content: three blocks of 50 and two 28 px gaps over 400 px.
    expect([x0s[0], x1s[2], v.navWidth, v.navX0, v.navX1].map((x) => +x.toFixed(9))).toEqual([40, 440, 400, 40, 440]);
    for (const cell of [minCell, Math.sqrt(minCell * 16), 4, 16]) {
      at(cell, 0);
      expect([Array.from(v.navX0s), Array.from(v.navX1s)]).toEqual([x0s, x1s]);
      expect(v.navWidth).toBeCloseTo(400, 9);
    }
    // The view's left edge at block 1's start: the window starts at its extent's start, at any zoom.
    for (const cell of [4, 16]) {
      at(cell, ax.offset[1] * cell + ax.gap);
      expect(v.navX0).toBeCloseTo(x0s[1], 9);
    }
    // At 16 px its right edge is 25 of block 1's positions further, at the fit's share of them.
    expect(v.navX1).toBeCloseTo(x0s[1] + (25 * (x1s[1] - x0s[1])) / 50, 9);
    // At 4 px it is past block 1 and the gap: 43 positions into block 2.
    at(4, ax.offset[1] * 4 + ax.gap);
    expect(v.navX1).toBeCloseTo(x0s[2] + (43 * (x1s[2] - x0s[2])) / 50, 9);
    // In the gap after block 0: proportionally across the map's gap.
    at(16, 50 * 16 + 14);
    expect(v.navX0).toBeCloseTo((x1s[0] + x0s[1]) / 2, 9);
  });

  it("navUnits inverts the window's map: units at any cell give back the content x", async () => {
    const { ax, v, at } = await synced();
    at(8, 0);
    const io = new Float64Array(2);
    for (const cx of [0, 120, 50 * 8 + 10, 50 * 8 + 28, 90 * 8 + 28 + 5, 150 * 8 + 56 - 1]) {
      io[0] = cx;
      // Content → map: the axis view's own map (`sync` with that scroll), then back.
      at(8, cx);
      io[0] = v.navX0;
      navUnits(v, io);
      expect(io[0] * 8 + io[1] * ax.gap).toBeCloseTo(cx, 6);
    }
  });

  it("the gestures go through it: a press jumps the window's center to the pointer; an end dragged holds the other", async () => {
    const t = await driven(blocksDoc(3, 400), { size: {} });
    const sec = t.ctx.viewports.sections.value[0];
    t.ctx.ui.setStyle(sec, "tracks");
    await t.ctx.graph.update();
    const vp = t.ctx.viewports.get(sec);
    vp?.zoomTo(0, 100, 160);
    const nav = () => t.frame(sec, WindowKind.navigator);
    const av = () => t.ctx.layout.axis(sec) as SectionAxisView;
    const y = 10;
    // A press in block 1's first quarter, right of the window: the window jumps, centered on it.
    const x = av().navX0s[1] + (av().navX1s[1] - av().navX0s[1]) / 4;
    t.pointer("down", sec, x, y);
    t.pointer("up", sec, x, y);
    expect(Math.abs((nav()[8] + nav()[9]) / 2 - x)).toBeLessThan(1);
    // Its right end dragged across a gap: the left end held, the right one under the pointer.
    const [x0, x1] = [nav()[8], nav()[9]];
    const to = av().navX0s[2] + 10;
    t.drag(sec, x1 - 2, y, to - 2, y);
    expect(Math.abs(nav()[8] - x0)).toBeLessThan(0.5);
    expect(Math.abs(nav()[9] - to)).toBeLessThan(0.5);
    // Its body dragged back by 30 px: the window follows the pointer.
    const m = (nav()[8] + nav()[9]) / 2;
    const before = nav()[8];
    t.drag(sec, m, y, m - 30, y);
    expect(Math.abs(nav()[8] - (before - 30))).toBeLessThan(1);
    // Its left end dragged into block 0: the right end held, the left one under the pointer.
    const [l0, l1] = [nav()[8], nav()[9]];
    const from = av().navX1s[0] - 20;
    t.drag(sec, l0 + 2, y, from + 2, y);
    expect(Math.abs(nav()[9] - l1)).toBeLessThan(0.5);
    expect(Math.abs(nav()[8] - from)).toBeLessThan(0.5);
  });

  it("a floored window (its share under 16 px, drawn wider) follows the pointer as drawn when its body is dragged", async () => {
    const t = await driven(blocksDoc(3, 20_000), { size: {} });
    const sec = t.ctx.viewports.sections.value[0];
    t.ctx.ui.setStyle(sec, "tracks");
    await t.ctx.graph.update();
    t.ctx.viewports.get(sec)?.zoomTo(1, 100, 110);
    const nav = () => t.frame(sec, WindowKind.navigator);
    const av = () => t.ctx.layout.axis(sec) as SectionAxisView;
    const [f0, f1] = [nav()[8], nav()[9]];
    expect(f1 - f0).toBeGreaterThan(av().navX1 - av().navX0 + 4);
    // From block 1's start across the gap after it into block 2.
    const to = av().navX0s[2] + 20;
    expect(f0).toBeLessThan(av().navX1s[1]);
    t.drag(sec, (f0 + f1) / 2, 10, (f0 + f1) / 2 + to - f0, 10);
    expect(Math.abs(nav()[8] - to)).toBeLessThan(1);
  });
});

describe("the navigator's draws", () => {
  it("stacked: the fade and grips span the map to the widest lane; marks per lane", async () => {
    const { v, at, minCell } = await synced(stackedDoc());
    at(minCell * 4, 0);
    const n = navView(v, navigatorFrom({ lanes: true, height: 30, thickness: 30 }).parts);
    n.height = 30;
    n.gripBottom = 30;
    expect([v.navX1s[2], v.navWidth].map((x) => +x.toFixed(9))).toEqual([+(40 + (40 / 120) * 400).toFixed(9), 400]);
    const marks = new SectionMarks();
    marks.selection.reset(3);
    marks.selection.begin(2);
    marks.selection.add(10, 20);
    marks.selection.close(0, 40);
    marks.updatePrimary(v, 3);
    n.marks = marks;
    drawNavigatorWindow(n);
    const b = rec(n);
    const fade = rects(b, SkinToken.surface).filter((r) => r[1] === 0);
    // The faded outside reaches the map's right edge (440), not the last block's (173).
    expect(Math.max(...fade.map((r) => r[0] + r[2]))).toBe(440);
    // Block 2's selection in lane 2 (rows [20, 30)), at its positions on its extent.
    expect(rects(b, SkinToken.selectionBorder).filter((r) => r[3] === 2)).toEqual([[40 + 33, 28, 34, 2]]);
  });

  it("the channels' runs as marks: highlight, accent, then selection, each its color, through the run converter", async () => {
    const { v, at, minCell } = await synced();
    at(minCell, 0);
    const n = navView(v);
    const marks = new SectionMarks();
    for (const [set, b, s, t] of [
      [marks.highlight, 0, 0, 10],
      [marks.accent, 1, 5, 6],
      [marks.selection, 2, 40, 50],
    ] as const) {
      set.reset(3);
      set.begin(b);
      set.add(s, t);
      set.close(0, 50);
    }
    marks.updatePrimary(v, 3);
    n.marks = marks;
    drawNavigatorWindow(n);
    const b = rec(n);
    const per = (v.navX1s[0] - v.navX0s[0]) / 50;
    const x = (bl: number, k: number) => Math.round(v.navX0s[bl] + k * per);
    expect(rects(b, SkinToken.highlight)).toEqual([[x(0, 0), 14, x(0, 10) - x(0, 0), 2]]);
    expect(rects(b, SkinToken.accent)).toEqual([[x(1, 5), 14, Math.max(1, x(1, 6) - x(1, 5)), 2]]);
    expect(rects(b, SkinToken.selectionBorder).filter((r) => r[3] === 2)).toEqual([
      [x(2, 40), 14, x(2, 50) - x(2, 40), 2],
    ]);
    // The marks' batch is handed after the window's (drawn over it), in channel order.
    const ops = Array.from({ length: b.opCount }, (_, k) => b.ops[k * 4]);
    expect(ops.filter((o) => o === Op.rects).length).toBeGreaterThanOrEqual(2);
  });

  it("on the ruler: only the hairline, even with a map to draw; off it, the map image", async () => {
    const { v, at, minCell } = await synced();
    at(minCell, 0);
    const rail = navView(v, DefaultNavigators[0].parts);
    rail.map = 7;
    drawNavigator(rail);
    expect([rec(rail).count(Op.image), rects(rec(rail), SkinToken.inkFaint).length, rail.overMap]).toEqual([
      0,
      3,
      false,
    ]);
    const strip = navView(v, DefaultNavigators[2].parts);
    strip.map = 7;
    drawNavigator(strip);
    expect([rec(strip).count(Op.image), strip.overMap]).toEqual([1, true]);
    const a = rec(strip).args;
    const o = rec(strip).ops[(rec(strip).opCount - 1) * 4 + 3];
    expect(Array.from(a.subarray(o, o + 8))).toEqual([0, 0, 400, 12, 40, 4, 400, 12]);
    // Over the map the window washes it and outlines itself in two tones.
    drawNavigatorWindow(strip);
    expect(rects(rec(strip), 0xffffff).length + rects(rec(strip), 0x111317).length).toBe(2);
  });
});

describe("grips", () => {
  /** The grips' pills in a color as [x, y, w, h] (device px), from the path store's round rects. */
  function pills(b: RecordingBackend, color: number) {
    const p = b.pathItems;
    const out: number[][] = [];
    for (let i = 0; i < p.n; i++) {
      if (p.color[i] !== color) continue;
      const c = p.cstart[i];
      out.push(Array.from(p.coords.subarray(c, c + 4)));
    }
    return out;
  }
  const grip = (vertical: boolean, lit: number, ends: number[], dpr = 1) => {
    const b = new RecordingBackend();
    gripAt.set([...ends, dpr]);
    drawGrips(b, vertical, lit);
    return b;
  };

  it("horizontal: centered on the window's ends, kept in the span, on a halo; the lit end in the hover tone", () => {
    const b = grip(false, 2, [100, 200, 40, 440, 10, 0, 20]);
    expect(pills(b, SkinToken.grip)).toEqual([[99, 5, 3, 10]]);
    expect(pills(b, SkinToken.gripHover)).toEqual([[199, 5, 3, 10]]);
    expect(pills(b, SkinToken.surface)).toEqual([
      [98, 4, 5, 12],
      [198, 4, 5, 12],
    ]);
    // At the full view, at the span's far ends.
    const full = grip(false, 0, [40, 440, 40, 440, 10, 0, 20]);
    expect(pills(full, SkinToken.grip)).toEqual([
      [40, 5, 3, 10],
      [437, 5, 3, 10],
    ]);
  });

  it("vertical: the overview's top and bottom grips along y, centered across; a short window draws them apart", () => {
    const b = grip(true, 3, [60, 160, 20, 220, 450, 440, 460], 2);
    // At DPR 2: 6 device px along, 20 across, centered on x 450.
    expect(pills(b, SkinToken.gripHover)).toEqual([
      [440, 57, 20, 6],
      [440, 157, 20, 6],
    ]);
    const short = grip(true, 0, [100, 102, 20, 220, 450, 440, 460], 2);
    const [g0, g1] = pills(short, SkinToken.grip);
    expect(g1[1] - (g0[1] + g0[3])).toBeGreaterThanOrEqual(4);
  });
});

describe("the keyed map image", () => {
  it("built once across pans; a key change builds it again, the old disposed", () => {
    const b = new RecordingBackend();
    const k = new KeyedImage(4);
    const frame = (layout: number, width: number) => {
      k.next[0] = layout;
      k.next[1] = width;
      k.next[2] = 1;
      if (k.stale(b.images)) {
        const layer = b.images.beginLayer(width, 12);
        layer.beginPass(0, width, 12, 1);
        layer.endPass();
        k.keep(b.images, b.images.endLayer());
      }
      return k.id;
    };
    const first = frame(1, 400);
    for (let i = 0; i < 100; i++) frame(1, 400);
    expect([k.builds, k.id]).toEqual([1, first]);
    frame(2, 400);
    expect([k.builds, b.images.content[first]]).toEqual([2, null]);
    // Another backend's images: built there.
    const other = new RecordingBackend();
    k.next[0] = 2;
    expect(k.stale(other.images)).toBe(true);
    k.dispose();
    expect(k.id).toBe(-1);
  });
});

it("the strip at DPR 2 draws its grounds per block on device pixels", async () => {
  const { v, at, minCell } = await synced();
  at(minCell, 0);
  const n = navView(v, DefaultNavigators[2].parts, 2);
  drawNavigator(n);
  const p = rec(n).pathItems;
  expect(p.n).toBe(3);
  expect(Array.from(p.coords.subarray(0, 5))).toEqual([80, 8, Math.round(v.navX1s[0] * 2) - 80, 24, 4]);
  void MarkBit;
});
