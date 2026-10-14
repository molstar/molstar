/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The default behaviors through input records (seqspec:viewer#behaviors, seqspec:interaction#selection): drag to pan
// with its locks, the selection's click policy and drag ranges, the navigator, the vertical overview, clades, the
// viewport's emissions during a drag, the drag default and its two-call flip, the wheel, presets' overrides, and
// behaviors removable without leaving subscriptions, and the keyboard's first case (keyboard.test.ts has the rest,
// IN4). The commands' part is CM1's; the actions', VC2's.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { OrderedSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { ARG_COLUMN, case2, case3, case10 } from "../../core/_test/format/cases.ts";
import { TRNA_STO } from "../../core/_test/format/load-fixtures.ts";
import {
  DefaultSeqViewerSpec,
  Hit,
  HitRegion,
  HoverKind,
  type LabelProvider,
  type PresetEntry,
  type Selection,
} from "../src/index.ts";
import { gripsLit, railGrown, WindowKind, WindowPart } from "../src/windows.ts";
import { type Driven, driven, subscriptions } from "./input-fixtures.ts";
import { PRESETS } from "./layout-fixtures.ts";

const S = "family.main";
const labeled: LabelProvider = {
  name: "lab:any",
  priority: 1000,
  applies: (h) => h.kind !== HoverKind.none,
  label: (h, _lc, out) => out.title(`kind ${h.kind}`),
};
async function family() {
  const t = await driven(case3(), { size: { [S]: [1000, 150] }, registry: [{ labels: [labeled] }] });
  const l = t.layout(S);
  return { ...t, row: (d: number) => l.scrolledTop + d * 18 + 9, columnY: l.entryTop[1] + 9 };
}
const columns = (s: Selection) => {
  const out: number[] = [];
  if (s.columns) OrderedSet.forEach(s.columns.columns, (i: number) => out.push(i), null);
  return s.columns ? out : null;
};

describe("drag to pan, with axis locks", () => {
  it("pans x only under a y lock; a bound dominantAxis keeps the dominant axis; cursor, ui.drag, hover and label follow", async () => {
    const t = await family();
    const vp = t.vp(S);
    vp.zoomTo(0, 5, 15);
    t.ctx.behaviors.get("dragPan")?.setParams({ x: true, y: false });
    const s = vp.live.current;
    const [x0, y0] = [s.scrollX, s.scrollY];
    t.pointer("move", S, 500, t.row(1));
    expect([t.ctx.ui.cursor(S).value, t.ctx.hover.live.current.owner, t.ctx.labels.current.value?.title]).toEqual([
      "grab",
      "pointer",
      `kind ${HoverKind.aligned}`,
    ]);
    t.pointer("down", S, 500, t.row(1));
    t.pointer("move", S, 480, t.row(1) - 10);
    t.pointer("move", S, 450, t.row(1) - 20);
    expect([s.scrollX - x0, s.scrollY - y0]).toEqual([50, 0]);
    expect([t.ctx.ui.drag.value?.owner, t.ctx.ui.cursor(S).value]).toEqual(["dragPan", "grabbing"]);
    expect([t.ctx.hover.live.current.owner, t.ctx.labels.current.value]).toEqual(["", null]);
    t.pointer("up", S, 450, 400, { onDocument: true }); // outside the surface: the drag ends there
    expect([t.ctx.ui.drag.value, t.ctx.ui.cursor(S).value]).toEqual([null, "default"]);
    t.pointer("move", S, 450, t.row(0));
    expect([t.ctx.hover.live.current.owner, t.ctx.ui.cursor(S).value]).toEqual(["pointer", "grab"]);
    expect(t.ctx.labels.current.value).not.toBeNull();

    // The dominant-axis lock has no default modifier (Shift-drag selects columns): bound here to Shift.
    t.ctx.behaviors.get("select")?.setParams({ dragColumns: [] });
    t.ctx.behaviors.get("dragPan")?.setParams({ y: true, dominantAxis: { shift: true } });
    const [x1, y1] = [s.scrollX, s.scrollY];
    t.drag(S, 500, t.row(1), 460, t.row(1) - 12, { shift: true });
    expect([s.scrollX - x1, s.scrollY - y1]).toEqual([40, 0]);
    t.drag(S, 500, t.row(1), 496, t.row(1) - 20, { shift: true });
    expect([s.scrollX - x1, s.scrollY - y1]).toEqual([40, 20]);
  });
});

describe("selection toggle", () => {
  it("a click on a column entry selects it, again clears it (a double click keeps it); Shift extends, Ctrl adds; a drag sets the range once", async () => {
    const t = await family();
    const sel = t.ctx.selection;
    const origins: string[] = [];
    sel.current.subscribe((s) => origins.push(s.origin));
    origins.length = 0;
    const x = (i: number) => t.xOf(S, 0, i);
    t.pointer("move", S, x(10), t.columnY);
    expect([t.ctx.hover.live.current.region, t.ctx.hover.live.current.kind]).toEqual([
      HitRegion.column,
      HoverKind.column,
    ]);
    t.click(S, x(10), t.columnY);
    expect(columns(sel.current.value)).toEqual([10]);
    // A double click keeps what its first click selected; a click again later clears it.
    t.click(S, x(10), t.columnY);
    expect(columns(sel.current.value)).toEqual([10]);
    t.advance(500);
    t.click(S, x(10), t.columnY);
    expect(sel.current.value.columns).toBeNull();
    // Shift on columns extends the range (from the last click's column), Ctrl/⌘ adds one.
    t.click(S, x(12), t.columnY, { shift: true });
    t.advance(500);
    t.click(S, x(10), t.columnY, { shift: true });
    expect(columns(sel.current.value)).toEqual([10, 11, 12]);
    t.advance(500);
    t.click(S, x(15), t.columnY, { ctrl: true });
    expect(columns(sel.current.value)).toEqual([10, 11, 12, 15]);
    t.advance(500);

    const pending: number[] = [];
    sel.pending.version.subscribe(() => pending.push(sel.pending.current.endShown));
    const emitted = origins.length;
    t.drag(S, x(10), t.columnY, x(14), t.columnY, {}, 8);
    expect(pending.slice(1, -1)).toEqual([10, 11, 12, 13, 14]); // the press's column, each one crossed; then the release
    expect(sel.pending.current.active).toBe(false);
    expect(origins.length - emitted).toBe(1);
    expect(columns(sel.current.value)).toEqual([10, 11, 12, 13, 14]);

    // A drag on rows pans.
    t.vp(S).zoomTo(0, 5, 15);
    const before = sel.current.value;
    const sx = t.vp(S).live.current.scrollX;
    t.drag(S, 500, t.row(0), 460, t.row(0));
    expect([sel.current.value, t.vp(S).live.current.scrollX - sx]).toEqual([before, 40]);
    expect(new Set(origins)).toEqual(new Set(["pointer"]));
    sel.set(null, "host:x");
    expect(origins.at(-1)).toBe("host:x");
  });

  it("a click on nothing clears; presets turn select off per section", async () => {
    const t = await family();
    t.click(S, t.xOf(S, 0, 3), t.columnY);
    t.click(S, 500, 149); // data, below the rows: no element
    // On data a double click zooms: the click's change waits for the double-click time.
    expect(t.ctx.selection.current.value.columns).not.toBeNull();
    t.advance(500);
    expect(t.ctx.selection.current.value.columns).toBeNull();
    const msa = PRESETS.find((p) => p.name === "msa") as PresetEntry;
    t.ctx.register({ presets: [{ ...msa, name: "lab:read", behaviors: { select: false, dragPan: { x: false } } }] });
    t.ctx.ui.setStyle(S, "lab:read");
    await t.ctx.graph.update();
    const s = t.ctx.layout.sectionIndex(S);
    expect(t.ctx.behaviors.get("select")?.forSection(s)).toBeNull();
    expect(t.ctx.behaviors.get<{ x: boolean; y: boolean }>("dragPan")?.forSection(s)).toMatchObject({
      x: false,
      y: true,
    });
    expect(t.ctx.behaviors.get("dragPan")?.forSection(s)).toBe(t.ctx.behaviors.get("dragPan")?.forSection(s));
    t.click(S, t.xOf(S, 0, 3), t.columnY);
    expect(t.ctx.selection.current.value.columns).toBeNull();
  });
});

describe("the navigator's gestures", () => {
  it("a press jumps, dragging the window pans, dragging an edge zooms", async () => {
    const t = await family();
    const vp = t.vp(S);
    vp.zoomTo(0, 0, 10);
    const s = vp.live.current;
    const nav = () => t.ctx.layout.axis(S) as NonNullable<ReturnType<typeof t.ctx.layout.axis>>;
    // At the most zoom the window is 79% of the axis: a press right of it jumps to the end.
    t.pointer("down", S, 900, 10);
    expect([t.ctx.ui.drag.value, s.scrollX]).toEqual([null, s.contentWidth - 700]);
    t.pointer("up", S, 900, 10);
    const sx = s.scrollX;
    const mid = (nav().navX0 + nav().navX1) / 2;
    t.drag(S, mid, 10, mid - 20, 10);
    expect(Math.abs(s.scrollX - sx - (-20 * s.contentWidth) / 700)).toBeLessThan(2); // each step snapped to a pixel
    expect(t.ctx.ui.drag.value).toBeNull();
    const cell = s.cell;
    const left = nav().navX0;
    t.drag(S, left, 10, left - 40, 10);
    expect(s.cell).toBeLessThan(cell);
    expect(Math.abs(nav().navX0 - (left - 40))).toBeLessThan(1);
  });
});

describe("a window's end dragged (windows.ts `windowEnds`)", () => {
  const frame = (t: Driven, sec: string, kind: number) => t.frame(sec, kind);
  /** A long sequence (case2) with a 20 px navigator row. */
  const F = "fv.main";
  it("moves only that end by the pointer's travel; stops at the least window and at the axis's end, never sliding", async () => {
    const t = await driven(case2(), { size: { [F]: [1000, 200] } });
    t.ctx.ui.setStyle(F, "tracks");
    await t.ctx.graph.update();
    const vp = t.vp(F);
    const s = vp.live.current;
    const nav = () => frame(t, F, 0);
    vp.zoomTo(0, 100, 200);
    const g = nav();
    const [lo, hi, x0, x1] = [g[0], g[0] + g[1], g[8], g[9]];
    expect([lo < x0 - 50, x1 + 50 < hi]).toEqual([true, true]);
    // Pressed 2 px inside the right end, then a 1 px move: the end moves 1 px, nothing jumps, the left end stays.
    t.pointer("down", F, x1 - 2, 10);
    t.pointer("move", F, x1 - 1, 10);
    expect(t.ctx.ui.drag.value?.owner).toBe("navigator");
    expect(Math.abs(nav()[9] - (x1 + 1))).toBeLessThan(0.5);
    expect(Math.abs(nav()[8] - x0)).toBeLessThan(0.5);
    // Past the left end: the window stops exactly at its least (the zoom limit's share), the left end unchanged.
    t.pointer("move", F, x0 - 100, 10);
    const min = nav();
    expect(Math.abs(min[8] - x0)).toBeLessThan(0.5);
    expect(s.cell).toBe(s.maxCell);
    expect(Math.abs(min[9] - min[8] - (min[5] - min[4]) * min[1])).toBeLessThan(1e-9);
    // 200 px further, in steps: nothing moves (no pan, no shift of either end).
    const at = [min[8], min[9], s.cell, s.scrollX];
    for (let k = 1; k <= 4; k++) t.pointer("move", F, x0 - 100 - 50 * k, 10);
    const still = nav();
    expect([still[8], still[9], s.cell, s.scrollX]).toEqual(at);
    // Back past the press: the end follows again; past the axis's end it clamps there, the left end still held.
    t.pointer("move", F, x1 + 30, 10);
    expect(Math.abs(nav()[9] - (x1 + 32))).toBeLessThan(0.5);
    t.pointer("move", F, hi + 200, 10);
    t.pointer("up", F, hi + 200, 10);
    const end = nav();
    expect(Math.abs(end[9] - hi)).toBeLessThan(0.5);
    expect(Math.abs(end[8] - x0)).toBeLessThan(0.5);

    // The left end: past the axis's start it clamps there, the right end held; past the right end it stops.
    vp.zoomTo(0, 100, 200);
    const h = nav();
    t.drag(F, h[8] + 1, 10, lo - 200, 10);
    const a = nav();
    expect([Math.abs(a[8] - lo) < 0.5, Math.abs(a[9] - h[9]) < 0.5, s.scrollX]).toEqual([true, true, 0]);
    t.drag(F, a[8] + 1, 10, hi + 300, 10, {}, 12);
    const b = nav();
    expect(Math.abs(b[9] - h[9])).toBeLessThan(0.5);
    expect(s.cell).toBe(s.maxCell);
  });

  it("the vertical overview's bottom handle: dragged past the top it stops there, the top held, nothing slides", async () => {
    const t = await family();
    const vp = t.vp(S);
    const s = vp.live.current;
    vp.zoomTo(0, 0, 22);
    const g = frame(t, S, 1);
    const [y0, y1, cx] = [g[10], g[11], (g[8] + g[9]) / 2];
    t.pointer("down", S, cx, y1 - 1);
    t.pointer("move", S, cx, y1);
    expect(t.ctx.ui.drag.value?.owner).toBe("verticalOverview");
    expect(Math.abs(frame(t, S, 1)[11] - (y1 + 1))).toBeLessThan(0.5);
    expect(Math.abs(frame(t, S, 1)[10] - y0)).toBeLessThan(0.5);
    t.pointer("move", S, cx, y0 - 50);
    const m = frame(t, S, 1);
    const at = [m[10], m[11], s.rowScale, s.scrollY];
    expect(Math.abs(m[10] - y0)).toBeLessThan(0.5);
    for (let k = 1; k <= 4; k++) t.pointer("move", S, cx, y0 - 50 - 50 * k);
    const still = frame(t, S, 1);
    expect([still[10], still[11], s.rowScale, s.scrollY]).toEqual(at);
    t.pointer("up", S, cx, y0 - 250);
  });
});

/** The seq loci of sequence `a`'s positions [start, stop) in a driven viewer's set. */
function lociOfRange(t: Driven, start: number, stop: number) {
  const set = t.ctx.graph.assembled.value?.set;
  if (!set) throw new Error("no set");
  return { kind: "seq-loci", set, elements: [{ unit: 0, indices: OrderedSet.ofRange(start, stop) }] } as never;
}

describe("the navigator's double click", () => {
  /** Two clicks within `Input.DoubleClickMs`. */
  const dblclick = (t: Driven, sec: string, x: number, y: number) => {
    t.click(sec, x, y);
    t.advance(50);
    t.click(sec, x, y);
  };
  /** The navigator's map of the axis, x0 and width (windows.ts `windowFrame`). */
  const map = (t: Driven, sec: string) => {
    const g = t.frame(sec, WindowKind.navigator);
    return [g[0], g[1]];
  };
  /** The position (in cells along the whole axis, gaps out) at a fraction of the content, or at the view's center. */
  const position = (t: Driven, sec: string, f?: number) => {
    const s = t.vp(sec).live.current;
    const ax = t.vp(sec).axis as NonNullable<ReturnType<typeof t.vp>["axis"]>;
    const cx = f === undefined ? s.scrollX + (s.width - ax.insetX) / 2 : f * s.contentWidth;
    let b = 0;
    while (b + 1 < ax.count && ax.offset[b + 1] * s.cell + (b + 1) * ax.gap <= cx) b++;
    return (cx - b * ax.gap) / s.cell;
  };

  const F = "fv.main";
  /** A sequence long enough that its fit is below the letters threshold, with a 20 px navigator row. */
  const fv = async () => {
    const t = await driven(case2(), { size: { [F]: [1000, 200] } });
    t.ctx.ui.setStyle(F, "tracks");
    await t.ctx.graph.update();
    return t;
  };

  it("at the full view zooms in, animated, to the letters threshold centered there; zoomed, back to the fit", async () => {
    const t = await fv();
    const S = F;
    const vp = t.vp(S);
    const s = vp.live.current;
    expect(s.cell).toBe(s.minCell);
    const [x0, w] = map(t, S);
    const x = x0 + 0.4 * w;
    t.pointer("move", S, x, 10);
    const cursor = t.ctx.ui.cursor(S).value;
    const at = position(t, S, 0.4);
    dblclick(t, S, x, 10);
    expect(s.animating).toBe(true);
    t.advance(100);
    // About 60 positions across the data area: the cell follows the view's width (at least the letters threshold).
    const target = Math.max(vp.props.minLetterCell, (s.width - (vp.axis?.insetX ?? 0)) / 60);
    expect([s.cell > s.minCell, s.cell < target]).toEqual([true, true]);
    t.advance(300);
    expect([s.animating, s.cell, s.lod]).toEqual([false, target, 0]);
    expect(Math.abs(position(t, S) - at)).toBeLessThan(0.2); // the view snaps to device pixels
    expect(t.ctx.ui.cursor(S).value).toBe(cursor);

    // Panned off the fit (a press beside the window jumps), then a double click there: the full view.
    dblclick(t, S, x0 + 0.9 * w, 10);
    t.advance(400);
    expect([s.cell, s.scrollX, s.animating]).toEqual([s.minCell, 0, false]);
    expect(t.ctx.ui.cursor(S).value).toBe(cursor);
  });

  it("zoomed in part way (wider than the focus window), in to the focus window there; then back to the fit", async () => {
    const t = await fv();
    const vp = t.vp(F);
    const s = vp.live.current;
    const target = Math.max(vp.props.minLetterCell, (s.width - (vp.axis?.insetX ?? 0)) / 60);
    vp.setView(Math.sqrt(s.minCell * target), 0, "api");
    expect([s.cell > s.minCell, s.cell < target]).toEqual([true, true]);
    const [x0, w] = map(t, F);
    dblclick(t, F, x0 + 0.3 * w, 10);
    t.advance(400);
    expect([s.animating, s.cell]).toEqual([false, target]);
    expect(Math.abs(position(t, F) - position(t, F, 0.3))).toBeLessThan(2);
    dblclick(t, F, x0 + 0.3 * w, 10);
    t.advance(400);
    expect([s.cell, s.scrollX]).toEqual([s.minCell, 0]);
  });

  it("on the data (`doubleClickZoom`): in to the focus window there, from part way in too; then back to the fit", async () => {
    const t = await fv();
    const vp = t.vp(F);
    const s = vp.live.current;
    const v = t.ctx.layout.axis(F);
    const l = t.layout(F);
    if (!v) throw new Error("no axis");
    const target = Math.max(vp.props.minLetterCell, (s.width - (vp.axis?.insetX ?? 0)) / 60);
    const y = l.scrolledTop + 4;
    const x = v.left + 0.25 * (s.width - (vp.axis?.insetX ?? 0));
    const at = position(t, F, (s.scrollX + (x - v.left)) / s.contentWidth);
    dblclick(t, F, x, y);
    t.advance(400);
    expect([s.animating, s.cell]).toEqual([false, target]);
    expect(Math.abs(position(t, F) - at)).toBeLessThan(2);
    dblclick(t, F, x, y);
    t.advance(400);
    expect([s.cell, s.scrollX]).toEqual([s.minCell, 0]);
    // Part way in: in to the focus window, not out.
    vp.setView(Math.sqrt(s.minCell * target), 0, "api");
    dblclick(t, F, x, y);
    t.advance(400);
    expect(s.cell).toBe(target);
    // Off: nothing moves.
    vp.fit({ animate: false, origin: "api" });
    t.ctx.behaviors.get("doubleClickZoom")?.setParams({ doubleClick: "none" });
    dblclick(t, F, x, y);
    t.advance(400);
    expect(s.cell).toBe(s.minCell);
  });

  it("on an element: its whole extent, then the focus window there, then the full view", async () => {
    const b = createSeqViewSpec({ title: "dbl" });
    const a = b.data.sequence({ ref: "a", value: "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ".repeat(10) });
    const d = b.data.annotation({ ref: "d", geometry: "ranges", target: a, value: { start: [100], stop: [180] } });
    const v = b.view({ ref: "v", label: "V" }).block({ sequence: a });
    v.track({ ref: "v.feat", label: "Features" }).ranges({ data: d });
    const t = await driven(b.getDocument(), { size: { "v.main": [1000, 200] } });
    const S = "v.main";
    const vp = t.vp(S);
    const s = vp.live.current;
    const l = t.layout(S);
    const ax = t.ctx.layout.axis(S);
    if (!ax) throw new Error("no axis");
    // The one track's entry (a scrolled entry: below the layout's scrolled top).
    const y = (0 < l.pinnedCount ? 0 : l.scrolledTop) + l.entryTop[0] + l.entryHeight[0] / 2;
    const xAt = (p: number) => ax.x(0, p as never) + s.cell / 2;
    const dbl = (x: number) => {
      t.click(S, x, y);
      t.advance(50);
      t.click(S, x, y);
      t.advance(600);
    };
    // 1: the range [100, 180) fills the view (with its padding).
    dbl(xAt(140));
    expect(vp.showsFocus(lociOfRange(t, 100, 180), 0.05)).toBe(true);
    expect(s.visLo[0] <= 100 && s.visHi[0] >= 180 && s.visHi[0] - s.visLo[0] < 100).toBe(true);
    // 2: already shown exactly: the focus window at the position.
    dbl(xAt(140));
    const target = Math.max(vp.props.minLetterCell, (s.width - (vp.axis?.insetX ?? 0)) / 60);
    expect(s.cell).toBeCloseTo(Math.min(s.maxCell, target), 5);
    // 3: at the focus window: the full view.
    dbl(xAt(Math.floor((s.visLo[0] + s.visHi[0]) / 2)));
    expect([s.cell, s.scrollX]).toEqual([s.minCell, 0]);
  });

  it('a single click still jumps; `doubleClick: "none"` turns it off; its zoom cell is a param', async () => {
    const t = await fv();
    const S = F;
    const vp = t.vp(S);
    const s = vp.live.current;
    const [x0, w] = map(t, S);
    vp.zoomTo(0, 0, 10);
    const cell = s.cell;
    t.click(S, x0 + 0.95 * w, 10);
    t.advance(1000);
    expect([s.cell, s.scrollX]).toEqual([cell, s.contentWidth - (s.width - (vp.axis?.insetX ?? 0))]);

    vp.fit();
    t.ctx.behaviors.get("navigator")?.setParams({ doubleClick: "none" });
    dblclick(t, S, x0 + 0.4 * w, 10);
    t.advance(1000);
    expect([s.cell, s.scrollX]).toEqual([s.minCell, 0]);

    t.ctx.behaviors.get("navigator")?.setParams({ doubleClick: "toggle", zoomCell: 14 });
    t.advance(1000);
    const at = position(t, S, 0.4);
    dblclick(t, S, x0 + 0.4 * w, 10);
    t.advance(1000);
    expect(s.cell).toBe(14);
    expect(Math.abs(position(t, S) - at)).toBeLessThan(0.2);
  });

  it("on the rulerRail's grab zone too", async () => {
    const R = "fv.main";
    const t = await driven(case2(), { size: { [R]: [1000, 200] } });
    t.ctx.ui.setStyle(R, "railed");
    await t.ctx.graph.update();
    const s = t.vp(R).live.current;
    const y = t.layout(R).entryHeight[0] - 3;
    const [x0, w] = map(t, R);
    dblclick(t, R, x0 + 0.5 * w, y);
    t.advance(1000);
    expect(s.cell).toBe(Math.max(t.vp(R).props.minLetterCell, (s.width - (t.vp(R).axis?.insetX ?? 0)) / 60));
    dblclick(t, R, x0 + 0.5 * w, y);
    t.advance(1000);
    expect(s.cell).toBe(s.minCell);
  });
});

describe("the cursor shows the gesture a press starts (one classification, windows.ts)", () => {
  it("navigator, vertical overview and data: hover cursor, the drag's cursor and what the drag does agree", async () => {
    const t = await family();
    const vp = t.vp(S);
    vp.zoomTo(0, 0, 22);
    const s = vp.live.current;
    const l = t.layout(S);
    const v = () => t.ctx.layout.axis(S) as NonNullable<ReturnType<typeof t.ctx.layout.axis>>;
    const cursor = () => t.ctx.ui.cursor(S).value;
    /** Hover, press, move, look, release: the cursors before and during, and what moved. */
    const at = (x: number, y: number, dx: number, dy: number) => {
      t.pointer("move", S, x, y);
      const before = cursor();
      const [sx, sy, cell, rows] = [s.scrollX, s.scrollY, s.cell, s.rowScale];
      t.pointer("down", S, x, y);
      const jumped = s.scrollX !== sx || s.scrollY !== sy;
      for (let i = 1; i <= 4; i++) t.pointer("move", S, x + (dx * i) / 4, y + (dy * i) / 4);
      const during = cursor();
      const did = jumped
        ? "jump"
        : s.cell !== cell
          ? "zoom"
          : s.rowScale !== rows
            ? "rows"
            : s.scrollX !== sx || s.scrollY !== sy
              ? "pan"
              : "none";
      t.pointer("up", S, x + dx, y + dy);
      vp.zoomTo(0, 0, 22);
      vp.setRowScale(1, 0);
      vp.scrollTo(s.scrollX, 0, "api");
      return [before, during, did];
    };
    const mid = () => (v().navX0 + v().navX1) / 2;
    expect(at(mid(), 10, 20, 0)).toEqual(["grab", "grabbing", "pan"]);
    // At the zoom limit and the axis's start, the left end stays dragged either way (clamped); the right end zooms out.
    expect(at(v().navX0, 10, -30, 0)).toEqual(["ew-resize", "ew-resize", "none"]);
    expect(at(v().navX0, 10, 30, 0)).toEqual(["ew-resize", "ew-resize", "none"]);
    expect(at(v().navX1, 10, 30, 0)).toEqual(["ew-resize", "ew-resize", "zoom"]);
    expect(at(v().left + (v().right - v().left) * 0.95, 10, 0, 0)).toEqual(["pointer", "pointer", "jump"]);
    // The windows as drawn (their true share) and hit (windows.ts: at least WINDOW_HIT long).
    const f = (kind: number) => t.frame(S, kind);
    const wx = (x: number) => f(1)[8] + (f(1)[9] - f(1)[8]) * x;
    expect(at(wx(0.5), f(1)[11], 0, 10)).toEqual(["ns-resize", "ns-resize", "rows"]);
    expect(at(wx(1), (f(1)[10] + f(1)[11]) / 2, 4, 0)).toEqual(["ew-resize", "ew-resize", "zoom"]);
    expect(at(500, t.row(1), -40, 0)).toEqual(["grab", "grabbing", "pan"]);
    // The scrollbar (`verticalNav: "both"`, rows overflowing): its thumb's ends zoom rows with `verticalZoom`.
    const bar = l.width - l.vnavWidth - 4;
    expect(at(bar, (f(2)[10] + f(2)[11]) / 2, 0, 6)).toEqual(["grab", "grabbing", "pan"]);
    expect(at(bar, f(2)[11], 0, 6)).toEqual(["ns-resize", "ns-resize", "rows"]);
  });

  it("the overview's and the scrollbar's lit parts (ui.window), and the grips each lights", async () => {
    const t = await family();
    const l = t.layout(S);
    const f = (kind: number) => t.frame(S, kind);
    const lit = () => [t.ctx.ui.window.current.kind, t.ctx.ui.window.current.part];
    const o = f(1);
    const [cx, cy] = [(o[8] + o[9]) / 2, (o[10] + o[11]) / 2];
    t.pointer("move", S, cx, o[11] - 1);
    expect(lit()).toEqual([1, WindowPart.bottom]);
    t.pointer("move", S, o[9] - 1, cy);
    expect(lit()).toEqual([1, WindowPart.right]);
    t.pointer("move", S, cx, cy);
    expect(lit()).toEqual([1, WindowPart.pan]);
    const b = f(2);
    t.pointer("move", S, l.width - l.vnavWidth - 4, b[10] + 1);
    expect(lit()).toEqual([2, WindowPart.top]);
    // The grips lit, as bits (1 the start's, 2 the end's) of the ends along y (`vertical`) or x.
    const W = WindowPart;
    expect([W.left, W.right, W.top, W.bottom, W.pan, W.page].map((p) => gripsLit(p, false))).toEqual([
      1, 2, 0, 0, 3, 0,
    ]);
    expect([W.left, W.right, W.top, W.bottom, W.pan, W.page].map((p) => gripsLit(p, true))).toEqual([0, 0, 1, 2, 3, 0]);
  });

  it("without verticalZoom the windows' ends only move them (and say so)", async () => {
    const t = await family();
    t.ctx.ui.setStyle(S, "msa-fixed");
    const l = t.layout(S);
    const v = () => t.ctx.layout.axis(S) as NonNullable<ReturnType<typeof t.ctx.layout.axis>>;
    const c = l.width - l.vnavWidth;
    const wx = c + (((v().navX0 + v().navX1) / 2 - v().left) / (v().right - v().left)) * l.vnavWidth;
    for (const [x, y] of [
      [wx, v().vnavY1 - 1],
      [l.width - l.vnavWidth - 4, v().vnavY1 - 1],
    ]) {
      t.pointer("move", S, x, y);
      const r = t.vp(S).live.current.rowScale;
      t.drag(S, x, y, x, y + 6);
      expect([t.ctx.ui.cursor(S).value === "ns-resize", t.vp(S).live.current.rowScale === r]).toEqual([false, true]);
    }
  });
});

const rail = (t: Driven, R: string, i: number) => railGrown(t.ctx.ui.drag.value, t.layout(R), t.ctx.pointer.current, i);

describe("the rulerRail navigator", () => {
  it("is the pinned ruler's bottom border: no row of its own; it hits as the navigator, with its cursors and gestures", async () => {
    const R = "fv.main";
    const t = await driven(case2(), { size: { [R]: [1000, 200] } });
    t.ctx.ui.setStyle(R, "railed");
    const l = t.layout(R);
    expect([l.navHeight, l.railEntry, l.entryTop[0]]).toEqual([0, 0, 0]);
    const vp = t.vp(R);
    vp.zoomTo(0, 10, 40);
    const v = t.ctx.layout.axis(R) as NonNullable<ReturnType<typeof t.ctx.layout.axis>>;
    const y = l.entryHeight[0] - 3;
    const mid = (v.navX0 + v.navX1) / 2;
    t.pointer("move", R, mid, y);
    expect(t.ctx.ui.cursor(R).value).toBe("grab");
    const sx = vp.live.current.scrollX;
    t.drag(R, mid, y, mid + 30, y);
    expect(vp.live.current.scrollX).toBeGreaterThan(sx);
    t.pointer("move", R, v.navX1, y);
    expect(t.ctx.ui.cursor(R).value).toBe("ew-resize");
  });

  it("stays grown for its whole drag, wherever the pointer goes; shrinks after the release off it", async () => {
    const R = "fv.main";
    const t = await driven(case2(), { size: { [R]: [1000, 200] } });
    t.ctx.ui.setStyle(R, "railed");
    t.vp(R).zoomTo(0, 10, 40);
    const v = t.ctx.layout.axis(R) as NonNullable<ReturnType<typeof t.ctx.layout.axis>>;
    const i = t.ctx.layout.sectionIndex(R);
    const y = t.layout(R).entryHeight[0] - 3;
    const mid = (v.navX0 + v.navX1) / 2;
    t.pointer("move", R, mid, 120);
    expect(rail(t, R, i)).toBe(false);
    t.pointer("move", R, mid, y);
    expect(rail(t, R, i)).toBe(true);
    t.pointer("down", R, mid, y);
    t.pointer("move", R, mid + 10, 150); // far below the ruler row
    t.pointer("move", R, mid + 20, 400, { onDocument: true }); // outside the surface
    expect([t.ctx.ui.drag.value?.owner, rail(t, R, i)]).toEqual(["navigator", true]);
    t.pointer("up", R, mid + 20, 150);
    expect(rail(t, R, i)).toBe(false);
  });

  it("names the lit window part with the cursor (ui.window): a grip, the body; the press's for its whole drag", async () => {
    const R = "fv.main";
    const t = await driven(case2(), { size: { [R]: [1000, 200] } });
    t.ctx.ui.setStyle(R, "railed");
    t.vp(R).zoomTo(0, 10, 40);
    const v = t.ctx.layout.axis(R) as NonNullable<ReturnType<typeof t.ctx.layout.axis>>;
    const i = t.ctx.layout.sectionIndex(R);
    const y = t.layout(R).entryHeight[0] - 3;
    const lit = () => {
      const w = t.ctx.ui.window.current;
      return [w.section, w.kind, w.part, t.ctx.ui.cursor(R).value];
    };
    const [x0, x1] = [v.navX0, v.navX1];
    t.pointer("move", R, x0 + 1, y);
    expect(lit()).toEqual([i, 0, WindowPart.left, "ew-resize"]);
    t.pointer("move", R, x1 - 1, y);
    expect(lit()).toEqual([i, 0, WindowPart.right, "ew-resize"]);
    t.pointer("move", R, (x0 + x1) / 2, y);
    expect(lit()).toEqual([i, 0, WindowPart.pan, "grab"]);
    t.pointer("move", R, x1 + 100, y);
    expect(lit()).toEqual([i, 0, WindowPart.page, "pointer"]);
    // Off the windows: none (the version bumps only on a change).
    const version = t.ctx.ui.window.version.value;
    t.pointer("move", R, 500, 150);
    t.pointer("move", R, 520, 160);
    expect([lit().slice(0, 3), t.ctx.ui.window.version.value - version]).toEqual([[-1, -1, WindowPart.none], 1]);
    // A drag of the right end keeps that grip lit, wherever the pointer goes, until the release.
    t.pointer("move", R, x1 - 1, y);
    t.pointer("down", R, x1 - 1, y);
    t.pointer("move", R, x1 + 40, 150);
    t.pointer("move", R, x1 + 400, 400, { onDocument: true });
    expect(lit().slice(0, 3)).toEqual([i, 0, WindowPart.right]);
    t.pointer("up", R, 500, 150);
    expect(lit().slice(0, 3)).toEqual([-1, -1, WindowPart.none]);
    // A press outside the window jumps there and pans it: the body is lit for the drag.
    t.pointer("move", R, 900, y);
    t.pointer("down", R, 900, y);
    t.pointer("move", R, 880, 150);
    expect(lit().slice(0, 3)).toEqual([i, 0, WindowPart.pan]);
    t.pointer("up", R, 880, 150);
  });

  const region = (t: Driven, sec: string, x: number, y: number) => {
    const h = Hit.create();
    t.ctx.layout.hitTest(sec, x, y, h);
    return h.region;
  };
  it("its grab zone: the whole ruler row on a blocks axis; on columns its lower 40%, the rest selects columns", async () => {
    const R = "fv.main";
    const t = await driven(case2(), { size: { [R]: [1000, 200] } });
    t.ctx.ui.setStyle(R, "railed");
    const l = t.layout(R);
    const h = l.entryHeight[0];
    for (const y of [0, h * 0.3, h * 0.7, h - 0.5]) expect([y, region(t, R, 500, y)]).toEqual([y, HitRegion.navigator]);
    expect(region(t, R, 500, h + 1)).toBe(HitRegion.data);
    expect(region(t, R, 50, h / 2)).toBe(HitRegion.header);

    const b = createSeqViewSpec({ title: "columns" });
    const aln = b.data.alignment({ ref: "aln", source: { value: TRNA_STO, format: "stockholm" } });
    const v = b.view({ ref: "av", label: "Alignment" }).block({ alignment: aln });
    v.track({ label: "" }).ruler();
    v.track({ label: "Rows" }).rows();
    const A = "av.main";
    const c = await driven(b.getDocument(), { size: { [A]: [1000, 200] } });
    c.ctx.ui.setStyle(A, "railed");
    const m = c.layout(A);
    const rh = m.entryHeight[0];
    expect([m.railEntry, m.navHeight]).toEqual([0, 0]);
    expect([region(c, A, 500, rh * 0.2), region(c, A, 500, rh * 0.5)]).toEqual([HitRegion.column, HitRegion.column]);
    expect([region(c, A, 500, rh * 0.7), region(c, A, 500, rh - 0.5)]).toEqual([
      HitRegion.navigator,
      HitRegion.navigator,
    ]);
  });
});

describe("clades", () => {
  it("a clade toggled collapses rows into a wedge; it survives an unrelated patch; a click on its header expands it", async () => {
    const t = await family();
    const rows = () => t.layout(S).rows.count;
    expect(rows()).toBe(4);
    t.ctx.ui.toggleClade("sh2", 1);
    expect([rows() < 4, Array.from(t.ctx.ui.clades.value.get("sh2") ?? [])]).toEqual([true, [1]]);
    await t.ctx.applyPatch({
      format: "seq_view_spec_patch",
      ops: [{ op: "set", ref: "family.logo", params: { label: "L" } }],
    });
    expect(Array.from(t.ctx.ui.clades.value.get("sh2") ?? [])).toEqual([1]);
    t.click(S, 50, t.row(0));
    expect([rows(), Array.from(t.ctx.ui.clades.value.get("sh2") ?? [])]).toEqual([4, []]);
  });
});

describe("the vertical overview", () => {
  it("claims the press; its window pans both axes; handles zoom rows and columns; elsewhere jumps", async () => {
    const t = await family();
    const vp = t.vp(S);
    vp.zoomTo(0, 0, 22);
    const s = vp.live.current;
    const l = t.layout(S);
    const f = () => t.frame(S, WindowKind.overview);
    const wx = () => [f()[8], 0] as const;
    const ww = () => f()[9] - f()[8];
    const owners: string[] = [];
    t.ctx.interaction.pointer.subscribe((e) => owners.push(e.owner), { priority: -1 });
    const [x0, y0] = [s.scrollX, s.scrollY];
    const cx = wx()[0] + ww() / 2;
    const cy = (f()[10] + f()[11]) / 2;
    const perPx = f()[12];
    t.pointer("down", S, cx, cy);
    t.pointer("move", S, cx + 2, cy + 5);
    t.pointer("up", S, cx + 2, cy + 5);
    expect(owners[0]).toBe("verticalOverview");
    expect(s.scrollX - x0).toBeCloseTo(2 * perPx * s.contentWidth, 0);
    expect(s.scrollY - y0).toBeGreaterThan(0);

    vp.scrollTo(0, 0, "api");
    const r = s.rowScale;
    t.drag(S, wx()[0] + ww() / 2, f()[11], wx()[0] + ww() / 2, f()[11] + 10);
    expect(s.rowScale).toBeLessThan(r);
    const cell = s.cell;
    t.drag(S, wx()[0] + ww(), cy, wx()[0] + ww() + 4, cy);
    expect(s.cell).toBeLessThan(cell);

    vp.setRowScale(1, 0);
    vp.scrollTo(s.scrollX, 0, "api");
    t.click(S, wx()[0] + ww() / 2, l.height - 3); // below the window: jumps there
    expect(s.scrollY).toBeGreaterThan(0);
  });
});

describe("subscribing to the viewport", () => {
  it("live emits per step of a drag; committed once, on release, with origin pointer and model blocks", async () => {
    const t = await driven(case10());
    const vp = t.vp("chains.main");
    vp.zoomTo(0, 200, 260, { origin: "api" });
    const committed: string[] = [];
    vp.committed.subscribe((c) => {
      if (c.origin !== "host:app") committed.push(`${c.origin} ${c.blocks.map((b) => `${b.ref}:${b.start}`).join()}`);
    });
    let lives = 0;
    vp.live.version.subscribe(() => lives++);
    [committed.length, lives] = [0, 0];
    t.pointer("down", "chains.main", 500, 100);
    for (let i = 1; i <= 60; i++) t.pointer("move", "chains.main", 500 - i, 100);
    expect(committed).toEqual([]);
    t.pointer("up", "chains.main", 440, 100);
    expect(lives).toBeGreaterThanOrEqual(56);
    expect(committed).toHaveLength(1);
    expect(committed[0]).toMatch(/^pointer 1iep\/A:\d+/);
    vp.zoomTo(0, 10, 20, { origin: "host:app" });
    expect([committed.length, vp.committed.value.origin]).toEqual([1, "host:app"]);
  });
});

describe("Shift-drag on data selects columns", () => {
  const cols = (t: Driven) => {
    const c = t.ctx.selection.current.value.columns;
    const got: number[] = [];
    if (c) OrderedSet.forEach(c.columns, (k, _i, out: number[]) => void out.push(k), got);
    return got;
  };
  it("whole columns without panning; a second Shift-drag extends from the anchor (`shiftPositions: extend`)", async () => {
    const t = await family();
    t.vp(S).zoomTo(0, 5, 20);
    const sx = t.vp(S).live.current.scrollX;
    t.drag(S, t.xOf(S, 0, 6), t.row(0), t.xOf(S, 0, 8), t.row(2), { shift: true });
    t.drag(S, t.xOf(S, 0, 12), t.row(1), t.xOf(S, 0, 13), t.row(1), { shift: true });
    expect(cols(t)).toEqual([6, 7, 8, 9, 10, 11, 12, 13]);
    expect(t.vp(S).live.current.scrollX).toBe(sx);
  });

  it("`shiftPositions: toggle`: a Shift-drag adds its columns, or removes them when pressed on a selected one", async () => {
    const t = await family();
    t.ctx.behaviors.get("select")?.setParams({ shiftPositions: "toggle" });
    t.vp(S).zoomTo(0, 5, 20);
    t.drag(S, t.xOf(S, 0, 6), t.row(0), t.xOf(S, 0, 8), t.row(2), { shift: true });
    t.drag(S, t.xOf(S, 0, 12), t.row(1), t.xOf(S, 0, 13), t.row(1), { shift: true });
    expect(cols(t)).toEqual([6, 7, 8, 12, 13]);
    t.drag(S, t.xOf(S, 0, 7), t.row(1), t.xOf(S, 0, 12), t.row(1), { shift: true });
    expect(cols(t)).toEqual([6, 13]);
  });
});

describe("policy by params and behaviors", () => {
  it("two calls flip the drag default; the two never both own a drag", async () => {
    const t = await family();
    t.vp(S).zoomTo(0, 5, 20);
    t.ctx.behaviors.get("select")?.setParams({ dragRegions: ["data", "column"], dragOn: [] });
    t.ctx.behaviors
      .get("dragPan")
      ?.setParams({ dragRegions: ["header", "gutter"], dragOn: [{ button: "primary", alt: true }] });
    const sx = t.vp(S).live.current.scrollX;
    t.drag(S, t.xOf(S, 0, 8), t.row(0), t.xOf(S, 0, 11), t.row(1));
    const a = t.ctx.selection.current.value.aligned;
    expect(a && !Array.isArray(a.elements)).toBe(true);
    expect(t.vp(S).live.current.scrollX).toBe(sx);
    t.drag(S, 600, t.row(0), 560, t.row(0), { alt: true });
    expect(t.vp(S).live.current.scrollX - sx).toBe(40);

    // 300 random drags: one owner at most, and only the owner acts.
    let seed = 7;
    const rnd = (n: number) => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed % n;
    };
    const owners = new Set<string>();
    t.ctx.interaction.drag.subscribe((e) => owners.add(e.owner), { priority: -100 });
    for (let k = 0; k < 300; k++) {
      t.vp(S).zoomTo(0, 5, 20);
      const [x, y] = [100 + rnd(860), 15 + rnd(130)];
      const sel = t.ctx.selection.current.value;
      const before = t.vp(S).live.current.scrollX;
      let owner = "";
      const sub = t.ctx.ui.drag.subscribe((d) => {
        if (d) owner = d.owner;
      });
      t.drag(S, x, y, x - 30 - rnd(40), y, { alt: rnd(2) === 1 });
      sub.unsubscribe();
      if (owner === "select") expect(t.vp(S).live.current.scrollX).toBe(before);
      if (owner === "dragPan") expect(t.ctx.selection.current.value).toBe(sel);
      t.advance(400);
    }
    expect(
      [...owners].every((o) => ["", "select", "dragPan", "navigator", "verticalOverview", "scrollbar"].includes(o)),
    ).toBe(true);
  });

  it("without select a click selects nothing but still emits; Focus moves to the FLVR arginine column", async () => {
    const t = await family();
    const clicks: number[] = [];
    t.ctx.interaction.click.subscribe((e) => clicks.push(e.hit.block));
    t.ctx.behaviors.get("select")?.remove();
    t.click(S, t.xOf(S, 0, 4), t.columnY);
    expect([t.ctx.selection.current.value.columns, clicks]).toEqual([null, [0]]);
    expect(t.ctx.focus.fromLoci(t.ctx.refs.toLoci("flvr-arg") as never)).toBe(true);
    t.advance(1000);
    const s = t.vp(S).live.current;
    expect(s.visLo[0] <= ARG_COLUMN && ARG_COLUMN < s.visHi[0] && s.cell > 22).toBe(true);
  });

  it("a behavior that throws in a subscriber is reported and stays; the others still run", async () => {
    const t = await family();
    const issues: string[] = [];
    t.ctx.events.issue.subscribe((i) => issues.push(i.code));
    t.ctx.behaviors.add({
      name: "lab:bad",
      display: { name: "Bad" },
      category: "interaction",
      priority: 60,
      params: {},
      create: (ctx) => {
        let sub: { unsubscribe(): void } | null = null;
        return {
          register: () => {
            sub = ctx.interaction.click.subscribe(() => {
              throw new Error("boom");
            });
          },
          dispose: () => sub?.unsubscribe(),
        };
      },
    });
    t.click(S, t.xOf(S, 0, 4), t.columnY);
    expect(issues).toEqual(["seqv.behavior.error"]);
    expect(t.ctx.behaviors.get("lab:bad")).not.toBeNull();
    expect(t.ctx.selection.current.value.columns).not.toBeNull();
  });
});

describe("the wheel and row zoom", () => {
  it("Ctrl+wheel zooms anchored at the cursor; a horizontal wheel pans; a box section's rows scroll", async () => {
    const t = await family();
    const s = t.vp(S).live.current;
    const l = t.layout(S);
    const at = (x: number) => (s.scrollX + x - l.dataLeft) / s.cell;
    const x = 600;
    const before = at(x);
    expect(t.wheel(S, x, t.row(0), 0, -200, { ctrl: true })).toBe(true);
    expect(s.cell).toBeGreaterThan(22);
    expect(at(x)).toBeCloseTo(before, 1);
    const sx = s.scrollX;
    expect(t.wheel(S, x, t.row(0), 30, 2)).toBe(true);
    expect(s.scrollX - sx).toBe(30);
    expect(t.wheel(S, x, t.row(0), 0, 20)).toBe(true);
    expect(s.scrollY).toBe(20);
    expect(t.wheel(S, x, t.row(0), 0, 100)).toBe(true); // clamped at the end
    expect(t.wheel(S, x, t.row(0), 0, 100)).toBe(false); // past the end: the page's
    const r = s.rowScale;
    expect(t.wheel(S, x, t.row(0), 0, -100, { alt: true })).toBe(true);
    expect(s.rowScale).toBeGreaterThan(r);
  });

  it("a content section passes the vertical wheel to the page", async () => {
    const t = await driven(case10());
    t.ctx.views.setCurrent("strip");
    await t.ctx.graph.update();
    expect(t.wheel("strip.main", 500, 30, 0, 50)).toBe(false);
  });
});

describe("the keyboard", () => {
  it("arrows move the keyboard's hover; Enter selects; Escape clears; Tab leaves the surface; + zooms", async () => {
    const t = await family();
    // The first key makes the cursor (the first position in view of the first scrolled row) and then acts.
    expect(t.key(S, "ArrowRight")).toBe(true);
    const h = t.ctx.hover.live.current;
    expect([h.owner, h.kind, h.aligned.column, h.dataColumn]).toEqual(["keyboard", HoverKind.aligned, 1, -1]);
    t.key(S, "ArrowUp");
    expect([h.kind, h.column.column]).toEqual([HoverKind.column, 1]); // onto the pinned columns above
    t.key(S, "Enter");
    expect(columns(t.ctx.selection.current.value)).toEqual([1]);
    t.key(S, "Escape");
    expect(t.ctx.selection.current.value.columns).toBeNull();
    // Tab is the browser's: it leaves the surface; there is no header focus to walk.
    expect([t.key(S, "Tab"), t.key(S, "Tab", { shift: true })]).toEqual([false, false]);
    const cell = t.vp(S).live.current.cell;
    t.key(S, "+");
    expect(t.vp(S).live.current.cell).toBeGreaterThan(cell);
    expect(t.key(S, "q")).toBe(false);
  });
});

describe("behaviors are removable", () => {
  it("removing every default behavior leaves no subscription they held", async () => {
    const bare = await driven(case3(), { size: { [S]: [1000, 150] } });
    const t = await driven(case3(), { size: { [S]: [1000, 150] } });
    for (const b of bare.ctx.behaviors.list.value) bare.ctx.behaviors.get(b.name)?.remove();
    expect(t.ctx.behaviors.list.value.map((b) => b.name)).toEqual(
      DefaultSeqViewerSpec().behaviors.map((b) => b.behavior.name),
    );
    expect(subscriptions(t.ctx)).toBeGreaterThan(subscriptions(bare.ctx));
    for (const b of t.ctx.behaviors.list.value) t.ctx.behaviors.get(b.name)?.remove();
    expect(subscriptions(t.ctx)).toBe(subscriptions(bare.ctx));
  });
});
