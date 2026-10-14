/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A resize draws what a fresh mount at the new size draws, pixel for pixel, on every canvas
// (seqspec:viewer#the-render-backend): the tracks, the overlay and the navigator's row, whose band is drawn only on a
// layout or size change (a pan does not redraw it), so a band drawn from the previous size's axis view stays wrong
// until the next layout change. The strip navigator: its miniature, the window's fade and ticks, over two blocks side
// by side. On `SectionSurface` with `readback: true`; the resize comes through the `ResizeObserver`, as a rail drag
// makes it.
// Found porting it: the composer sized the navigator's map image from the view's band as the last pass left it, so a
// resized surface drew its miniature 12 px thick and a fresh one 8 (other bar positions and widths); it sets the row first.

import { afterEach, describe, expect, it } from "vitest";
import type { PresetEntry } from "../../src/index.js";
import { twoChains } from "./chains.js";
import { canvasesOf, close, diff, type Opened, opened, pixelsOf, settle } from "./open.js";

const open: Opened[] = [];
afterEach(() => {
  for (const t of open.splice(0)) close(t);
});

const STRIP: PresetEntry = { name: "lab:strip", geometry: "linear", navigator: "strip" };
const H = 160;

async function strip(width: number, dpr: number) {
  const t = await opened(twoChains({ length: 300, tracks: 1, layers: ["ranges"] }), {
    width,
    height: H,
    dpr,
    registry: [{ presets: [STRIP] }],
  });
  open.push(t);
  t.ctx.ui.setStyle(t.section, STRIP.name);
  await t.ctx.graph.update();
  settle(t);
  return t;
}
/** The ResizeObserver's callback after a layout. */
const observed = () => new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r)));
const NAMES = ["tracks", "overlay", "navigator", "navigator overlay", "vertical overview", "its overlay"];

describe("a resize equals a fresh mount, pixel for pixel", () => {
  for (const dpr of [1, 2])
    it(`the strip navigator's miniature, fade and ticks over two blocks, zoomed and at the fit (${dpr} device px per CSS px)`, async () => {
      const t = await strip(900, dpr);
      const vp = t.ctx.viewports.get(t.section);
      if (!vp) throw new Error("no viewport");
      let compared = 0;
      for (const [w, zoom] of [
        [640, true],
        [1000, true],
        [700, false],
        [900, false],
      ] as const) {
        if (zoom) vp.zoomTo(0, 40, 160);
        else vp.fit();
        settle(t);
        t.host.style.width = `${w}px`;
        await observed();
        settle(t);
        const resized = canvasesOf(t).map(pixelsOf);
        const s = vp.live.current;
        expect(s.width).toBe(w);
        // A fresh mount at the new size, at the same view.
        const f = await strip(w, dpr);
        f.ctx.viewports.get(f.section)?.setView(s.cell, s.scrollX, "api");
        settle(f);
        const fresh = canvasesOf(f).map(pixelsOf);
        for (let k = 0; k < fresh.length; k++) {
          const c = canvasesOf(f)[k];
          const d = diff(resized[k], fresh[k], c.width);
          expect(d.count, `${w} px${zoom ? ", zoomed" : ""}: ${NAMES[k]} at (${d.x}, ${d.y})`).toBe(0);
          compared++;
        }
        close(f);
        open.splice(open.indexOf(f), 1);
      }
      expect(compared).toBe(4 * 6);
    });
});
