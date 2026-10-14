/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A navigator with nowhere to go (seqspec:interaction#viewports): the whole axis in view at the zoom
// limit is drawn quiet (no grips, no hover), its hit still reports the region, and its gestures and cursor are inert;
// it is live again as soon as the view can zoom (a narrower view), with no relayout.

import { describe, expect, it } from "vitest";
import type { SectionAxisView } from "../src/axis.ts";
import { SkinToken } from "../src/draw/api.ts";
import { drawNavigatorWindow, NavView } from "../src/draw/nav.ts";
import { DefaultNavigators, Hit, HitRegion } from "../src/index.ts";
import { fixtureMeasurer, RecordingBackend } from "../src/testing/index.ts";
import { TextManager } from "../src/text.ts";
import { windowIdle } from "../src/windows.ts";
import { blocksDoc } from "./axis-fixtures.ts";
import { driven } from "./input-fixtures.ts";

/** A 20-position section in a wide view: the cell is at its limit and the whole axis in view. */
async function idle() {
  const t = await driven(blocksDoc(1, 20), { size: {} });
  const sec = t.ctx.viewports.sections.value[0];
  t.ctx.ui.setStyle(sec, "tracks");
  await t.ctx.graph.update();
  const vp = t.vp(sec);
  return { t, sec, vp, s: vp.live.current };
}

describe("a navigator with nowhere to go", () => {
  it("windowIdle: the whole axis in view and zoom-in at its limit; a fixed cell; not without an axis or width", () => {
    const s = { cell: 28, maxCell: 28, contentWidth: 840, width: 1000 };
    const ax = { insetX: 0, fixedCell: 0 };
    expect(windowIdle(s, ax)).toBe(true);
    expect(windowIdle({ ...s, cell: 20 }, ax)).toBe(false);
    expect(windowIdle({ ...s, contentWidth: 1200 }, ax)).toBe(false);
    expect(windowIdle({ ...s, contentWidth: 1000 }, ax)).toBe(true);
    expect(windowIdle({ ...s, contentWidth: 1200 }, { insetX: 0, fixedCell: 28 })).toBe(true);
    expect(windowIdle({ ...s, width: 0 }, ax)).toBe(false);
    expect(windowIdle(s, null)).toBe(false);
  });

  it("gestures and the cursor are inert at max zoom and the full window; live again in a narrower view", async () => {
    const { t, sec, vp, s } = await idle();
    expect(windowIdle(s, vp.axis)).toBe(true);
    const y = 10;
    const hit = Hit.create();
    t.ctx.layout.hitTest(sec, 500, y, hit);
    expect(hit.region).toBe(HitRegion.navigator);
    const state = () => [s.cell, s.scrollX, s.animating, t.ctx.ui.drag.value?.owner ?? ""];
    const before = state();
    // Hover: the default cursor and no window part lit.
    t.pointer("move", sec, 500, y);
    expect([t.ctx.ui.cursor(sec).value, t.ctx.ui.window.current.part]).toEqual(["default", 0]);
    // A press and a drag (body, an end), a click to jump and a double click do nothing.
    t.pointer("down", sec, 500, y);
    t.pointer("move", sec, 400, y);
    expect(t.ctx.ui.drag.value?.owner).not.toBe("navigator");
    t.pointer("up", sec, 400, y);
    t.drag(sec, 500, y, 300, y);
    t.drag(sec, s.width - 1, y, s.width - 200, y);
    t.click(sec, 100, y);
    t.advance(50);
    t.click(sec, 100, y);
    t.advance(1000);
    expect(state()).toEqual(before);
    t.pointer("move", sec, 500, y);
    expect(t.ctx.ui.cursor(sec).value).toBe("default");
    // Data keeps its own cursor.
    t.pointer("move", sec, 500, 100);
    expect(t.ctx.ui.cursor(sec).value).not.toBe("pointer");

    // A narrower view can zoom out of the limit: the navigator is live (no relayout of its own).
    vp.setSize(300, 400, 1);
    t.advance(1000);
    expect(windowIdle(s, vp.axis)).toBe(false);
    expect(s.minCell).toBeLessThan(s.maxCell);
    t.pointer("move", sec, 250, y);
    expect(t.ctx.ui.cursor(sec).value).not.toBe("default");
    t.pointer("down", sec, 250, y);
    t.pointer("move", sec, 240, y);
    expect(t.ctx.ui.drag.value?.owner).toBe("navigator");
    t.pointer("up", sec, 240, y);
    t.click(sec, 250, y);
    t.advance(50);
    t.click(sec, 250, y);
    t.advance(1000);
    expect(s.cell).toBeGreaterThan(s.minCell);

    // Wide again: idle again, the cursor default.
    vp.setSize(1000, 400, 1);
    vp.fit();
    t.advance(1000);
    expect(windowIdle(s, vp.axis)).toBe(true);
    t.pointer("move", sec, 501, y);
    expect(t.ctx.ui.cursor(sec).value).toBe("default");
  });
});

describe("the idle window's draw", () => {
  const view = async (idleOn: boolean, hot = false) => {
    const t = await driven(blocksDoc(2, 30), { size: {} });
    const sec = t.ctx.viewports.sections.value[0];
    const v = t.ctx.layout.axis(sec) as SectionAxisView;
    const n = new NavView().from(v);
    n.backend = new RecordingBackend();
    n.text = new TextManager(fixtureMeasurer(), () => 64);
    n.parts = DefaultNavigators[0].parts;
    n.top = 0;
    n.height = 8;
    n.x0 = v.navX0;
    n.x1 = v.navX1;
    n.gripTop = 0;
    n.gripBottom = 8;
    n.lit = 3;
    n.hot = hot;
    n.idle = idleOn;
    drawNavigatorWindow(n);
    return n.backend as RecordingBackend;
  };
  const colors = (b: RecordingBackend, which: "rect" | "path") => {
    const items = which === "rect" ? b.rectItems : b.pathItems;
    return Array.from({ length: items.n }, (_, i) => items.color[i]);
  };

  it("no grips, no window line over the track, no hover growth; live, the same call draws them", async () => {
    const on = await view(true, true);
    const off = await view(false, true);
    expect(colors(on, "path")).toEqual([]);
    expect(colors(off, "path").length).toBe(4);
    expect(colors(off, "rect")).toContain(SkinToken.inkQuiet);
    expect(colors(on, "rect")).not.toContain(SkinToken.inkQuiet);
    // The track keeps its space and tone: the quiet line, thinner than the grown one.
    expect(colors(on, "rect")).toContain(SkinToken.inkFaint);
    const h = (b: RecordingBackend) => b.rectItems.h[colors(b, "rect").indexOf(SkinToken.inkFaint)];
    expect(h(on)).toBeLessThan(h(off));
  });
});
