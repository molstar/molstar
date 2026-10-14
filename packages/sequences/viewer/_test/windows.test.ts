/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Window classification (VC3; seqspec:viewer#behaviors, seqspec:interaction#selection), headless: the frame a section's
// axis view and layout give, the drawn window's floor and its hit area, an end dragged (it stops at the least window,
// the floor and the span's end, the other end held, nothing sliding), the part under a point (ends just outside and a
// quarter inside,
// row ends only with `verticalZoom`), and the grips a part lights. The gestures that drive them are IN2's.

import { describe, expect, it } from 'vitest';
import { AxisGeometry, SectionAxisView } from '../src/axis.js';
import {
  windowEnds as ends,
  gripsLit,
  WINDOW_EDGE,
  WINDOW_FLOOR,
  WINDOW_FRAME,
  WINDOW_HIT,
  WindowKind,
  WindowPart,
  windowFrame,
  windowPart,
  windowRect,
} from '../src/windows.js';
import { axisOf, blocksDoc, resolved } from './axis-fixtures.js';

/** A frame over a span [x0, x0 + w) × [y0, y0 + h) with window fractions along x and y. */
function frame(x0: number, w: number, fx: [number, number], y0 = 0, h = 100, fy: [number, number] = [0, 1]) {
  const g = new Float64Array(WINDOW_FRAME);
  g.set([x0, w, y0, h, fx[0], fx[1], fy[0], fy[1]]);
  windowRect(g);
  return g;
}
const at = (x: number, y: number) => ({ x, y });

describe('the drawn window and its hit area', () => {
  it('at its share of the span, at least WINDOW_FLOOR (at most half the span), its position over the travel', () => {
    const half = frame(10, 200, [0.25, 0.75]);
    expect([half[8], half[9]]).toEqual([60, 160]);
    // 1% of 200 px is 2 px: floored to 16, a thumb whose position maps the scroll over the span less the thumb.
    const tiny = frame(10, 200, [0.495, 0.505]);
    expect(tiny[9] - tiny[8]).toBe(WINDOW_FLOOR);
    expect(tiny[8]).toBeCloseTo(10 + (0.495 / 0.99) * (200 - WINDOW_FLOOR), 9);
    const end = frame(10, 200, [0.99, 1]);
    expect(end[9]).toBe(210);
    // A span of 20 px: the floor is half of it.
    expect(frame(0, 20, [0.4, 0.45])[9] - frame(0, 20, [0.4, 0.45])[8]).toBe(10);
  });

  it('hit over at least WINDOW_HIT px, centered on it and inside the span', () => {
    const g = frame(0, 400, [0.5, 0.505]);
    expect(g[15] - g[14]).toBe(WINDOW_HIT >= WINDOW_FLOOR ? WINDOW_HIT : WINDOW_FLOOR);
    const big = frame(0, 400, [0.25, 0.75]);
    expect([big[14], big[15]]).toEqual([big[8], big[9]]);
    const atStart = frame(0, 400, [0, 0.001]);
    expect(atStart[14]).toBe(0);
  });
});

describe('an end dragged (windowEnds)', () => {
  const out = new Float64Array(2);
  /** The travel and least in, the fractions out (`windowEnds` reads its doubles from the buffer). */
  const windowEnds = (
    g: Float64Array,
    vertical: boolean,
    start: boolean,
    total: number,
    least: number,
    io: Float64Array,
  ) => {
    io[0] = total;
    io[1] = least;
    return ends(g, vertical, start, io);
  };
  it("moves only that end by the travel; stops at the least window and at the span's end, never sliding", () => {
    const g = frame(0, 400, [0.25, 0.5]);
    expect(windowEnds(g, false, false, 20, 0.01, out)).toBe(true);
    expect([out[0], out[1]]).toEqual([0.25, 0.55]);
    // Past the left end: stops at the least window (5% of the span), the left end held.
    expect(windowEnds(g, false, false, -300, 0.05, out)).toBe(true);
    expect([out[0], out[1]]).toEqual([0.25, 0.3]);
    // ... or at the floor, when the least is shorter.
    windowEnds(g, false, false, -300, 0.001, out);
    expect(out[1] * 400 - out[0] * 400).toBeCloseTo(WINDOW_FLOOR, 9);
    // Past the span's end: clamped there.
    windowEnds(g, false, false, 1000, 0.01, out);
    expect([out[0], out[1]]).toEqual([0.25, 1]);
    // The start's end past the span's start, and past the other end.
    windowEnds(g, false, true, -1000, 0.01, out);
    expect([out[0], out[1]]).toEqual([0, 0.5]);
    windowEnds(g, false, true, 1000, 0.05, out);
    expect([out[0], out[1]]).toEqual([0.45, 0.5]);
  });

  it('a floored thumb dragged in keeps the view; dragged out, zooms out from its drawn ends', () => {
    const g = frame(0, 400, [0.5, 0.501]);
    expect(windowEnds(g, false, false, -10, 0.001, out)).toBe(false);
    expect(windowEnds(g, false, false, 40, 0.001, out)).toBe(true);
    expect(out[1] - out[0]).toBeCloseTo((WINDOW_FLOOR + 40) / 400, 9);
  });

  it("along y: the overview's bottom end dragged past the top stops at the floor, the top held", () => {
    const g = frame(0, 50, [0, 1], 20, 200, [0.2, 0.6]);
    expect(windowEnds(g, true, false, -500, 0, out)).toBe(true);
    expect(out[0]).toBe(0.2);
    expect((out[1] - out[0]) * 200).toBeCloseTo(WINDOW_FLOOR, 9);
  });
});

describe('the part under a point (windowPart)', () => {
  it('the navigator: its ends just outside and at most a quarter inside; its body pans; outside jumps', () => {
    const g = frame(0, 400, [0.25, 0.5]);
    const [x0, x1] = [g[8], g[9]];
    const nav = (x: number) => windowPart(WindowKind.navigator, at(x, 50), g, true);
    expect(nav((x0 + x1) / 2)).toBe(WindowPart.pan);
    expect([nav(x0 - WINDOW_EDGE), nav(x0 + WINDOW_EDGE - 1), nav(x0 - WINDOW_EDGE - 1)]).toEqual([
      WindowPart.left,
      WindowPart.left,
      WindowPart.page,
    ]);
    expect([nav(x1 + WINDOW_EDGE), nav(x1 - WINDOW_EDGE + 1)]).toEqual([WindowPart.right, WindowPart.right]);
    // A short window: a quarter of it inside each end, so its middle still pans.
    const s = frame(0, 400, [0.5, 0.55]);
    const mid = (s[8] + s[9]) / 2;
    expect(windowPart(WindowKind.navigator, at(mid, 0), s, false)).toBe(WindowPart.pan);
    expect(windowPart(-1, at(mid, 0), s, false)).toBe(WindowPart.none);
  });

  it('the scrollbar and the overview: row ends only with verticalZoom, else the window only moves', () => {
    const g = frame(300, 20, [0, 1], 0, 200, [0.4, 0.6]);
    const [y0, y1] = [g[10], g[11]];
    const bar = (y: number, rows: boolean) => windowPart(WindowKind.scrollbar, at(310, y), g, rows);
    expect([bar(y0, true), bar(y1, true), bar((y0 + y1) / 2, true), bar(10, true)]).toEqual([
      WindowPart.top,
      WindowPart.bottom,
      WindowPart.pan,
      WindowPart.page,
    ]);
    expect([bar(y0, false), bar(y1, false)]).toEqual([WindowPart.pan, WindowPart.pan]);
    // The overview: x ends and y ends; a corner is the x end's.
    const o = frame(300, 60, [0.25, 0.75], 0, 200, [0.4, 0.6]);
    const ov = (x: number, y: number) => windowPart(WindowKind.overview, at(x, y), o, true);
    const [cx, cy] = [(o[8] + o[9]) / 2, (o[10] + o[11]) / 2];
    expect([ov(cx, cy), ov(o[8], cy), ov(cx, o[10]), ov(cx, o[11]), ov(o[8], o[10]), ov(cx, 5)]).toEqual([
      WindowPart.pan,
      WindowPart.left,
      WindowPart.top,
      WindowPart.bottom,
      WindowPart.left,
      WindowPart.page,
    ]);
  });

  it('the grips a part lights', () => {
    expect([
      gripsLit(WindowPart.pan, false),
      gripsLit(WindowPart.left, false),
      gripsLit(WindowPart.right, false),
    ]).toEqual([3, 1, 2]);
    expect([gripsLit(WindowPart.top, false), gripsLit(WindowPart.top, true), gripsLit(WindowPart.page, true)]).toEqual([
      0, 1, 0,
    ]);
  });
});

describe("a section's frame from its axis view (windowFrame)", () => {
  it("the navigator maps the whole axis; zoomed, the window is the view's share of it", async () => {
    const r = await resolved(blocksDoc(2, 50), 'v');
    const ax = axisOf(r);
    const v = new SectionAxisView();
    const l = { geometry: AxisGeometry.linear, dataLeft: 40, dataRight: 240, scrolledTop: 20, height: 120 };
    const state = {
      ...{ cell: 2, scrollX: 0, scrollY: 0, dpr: 1, lod: 0, contentHeight: 0 },
      ...{ visLo: new Int32Array(2), visHi: new Int32Array(2) },
    };
    v.sync(ax, l, state as never);
    const g = new Float64Array(WINDOW_FRAME);
    windowFrame(v, { width: 260, height: 120, scrolledTop: 20, vnavWidth: 20 }, WindowKind.navigator, g);
    // 100 cells at 2 px and one 28 px gap: 228 px of content, 200 of data area.
    // Its map as at the fit: here the cell itself, so linear (to rounding).
    expect([g[0], g[1], g[4], +g[5].toFixed(12)]).toEqual([40, 200, 0, +(200 / 228).toFixed(12)]);
    v.sync(ax, l, { ...state, scrollX: 28 } as never);
    windowFrame(v, { width: 260, height: 120, scrolledTop: 20, vnavWidth: 20 }, WindowKind.navigator, g);
    expect([g[4], g[5]].map((f) => f * 228)).toEqual([28, 228]);
    windowFrame(v, { width: 260, height: 120, scrolledTop: 20, vnavWidth: 20 }, WindowKind.overview, g);
    expect([g[0], g[1], g[2], g[3]]).toEqual([240, 20, 20, 100]);
  });

  it('no data area (the header and data columns fill the width): the windows floored at the start, never NaN', async () => {
    const r = await resolved(blocksDoc(2, 50), 'v');
    const ax = axisOf(r);
    const v = new SectionAxisView();
    // The data area empty at x 200; 400 px of rows, 100 in view, scrolled by 100.
    const l = { geometry: AxisGeometry.linear, dataLeft: 200, dataRight: 200, scrolledTop: 20, height: 120 };
    const state = {
      ...{ cell: 2, minCell: 2, scrollX: 0, scrollY: 100, dpr: 1, lod: 0, contentHeight: 400 },
      ...{ visLo: new Int32Array(2), visHi: new Int32Array(2) },
    };
    v.sync(ax, l, state as never);
    const g = new Float64Array(WINDOW_FRAME);
    windowFrame(v, { width: 260, height: 120, scrolledTop: 20, vnavWidth: 60 }, WindowKind.overview, g);
    expect(g.every(Number.isFinite)).toBe(true);
    // Across the overview: a floored window at its start, 16 px of its 60; down it: the rows in view, at their share.
    expect([g[8], g[9]]).toEqual([200, 200 + WINDOW_FLOOR]);
    expect([g[10], g[11]]).toEqual([20 + 25, 20 + 50]);
    windowFrame(v, { width: 260, height: 120, scrolledTop: 20, vnavWidth: 60 }, WindowKind.navigator, g);
    expect(g.every(Number.isFinite)).toBe(true);
  });
});
