/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The dot plot as pixels (seqspec:format#layer-types `dot_plot`), in Chromium: a grid's block of high values over a
// low ground, each cell where x and y place it (y the section's window down the track), at 1 and 2 device px per CSS px
// and zoomed in, so a cell covers exactly its device pixels; the hover's row across the overlay.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { afterEach, describe, expect, it } from "vitest";
import type { SectionAxisView } from "../../src/axis.js";
import { Hit } from "../../src/index.js";
import { canvasesOf, close, type Opened, opened, pixelsOf, settle } from "./open.js";

let t: Opened | null = null;
afterEach(() => {
  if (t) close(t);
  t = null;
});

const N = 40;
/** Low (0) everywhere but i in [10, 20) × j in [25, 30): 10. */
function doc() {
  const b = createSeqViewSpec({ title: "pixels" });
  const s = b.data.sequence({ ref: "s", value: "MKTAYIAKQR".repeat(N / 10) });
  const value = Array.from({ length: N * N }, (_, x) => {
    const [i, j] = [Math.floor(x / N), x % N];
    return i >= 10 && i < 20 && j >= 25 && j < 30 ? 10 : 0;
  });
  const g = b.data.annotation({
    ref: "g",
    geometry: "grid",
    value: { target: ["s"], b_target: ["s"], value: [value] },
  });
  const palette = {
    kind: "discrete",
    mode: "absolute",
    colors: [
      ["#ff0000", 0],
      ["#0000ff", 5],
    ],
  };
  const r = b.data.annotation({
    ref: "r",
    geometry: "ranges",
    target: s,
    value: { start: [2], stop: [8], label: ["x"] },
  });
  const sec = b.view({ ref: "v", label: "V" }).section({ ref: "v.s" }).block({ sequence: s });
  // A feature drawn first leaves its row in the shared theme location: the plot's cells must not read it.
  sec.track({ ref: "t-r", label: "Features" }).ranges({ data: r, labels: true });
  sec
    .track({ ref: "t", label: "Plot" })
    .dot_plot({ data: g, color: { theme: "field_palette", params: { field: "value", palette } } });
  return b.getDocument();
}

describe("the dot plot's pixels", () => {
  for (const [dpr, cell] of [
    [1, 0],
    [2, 0],
    [2, 30],
  ] as const)
    it(`a cell where x and y place it, its device pixels exactly (${dpr} device px per CSS px, cell ${cell || "fit"})`, async () => {
      t = await opened(doc(), { width: 900, height: 400, dpr });
      if (cell) {
        t.ctx.viewports.get(t.section)?.setView(cell, 8 * cell, "api");
        settle(t);
      }
      const main = canvasesOf(t)[0];
      const px = pixelsOf(main);
      const v = t.ctx.layout.axis(t.section) as SectionAxisView;
      const l = t.ctx.layout.get(t.section)?.current;
      if (!l) throw new Error("no layout");
      const top = l.scrolledTop + l.entryTop[l.entryCount - 1];
      const h = l.entryHeight[l.entryCount - 1];
      const w = v.right - v.left;
      /** Cell (i, j)'s device pixel at a fraction (fx, fy) into it. */
      const at = (i: number, j: number, fx = 0.5, fy = 0.5) => {
        const x0 = v.x(0, i as never);
        const x1 = v.x(0, (i + 1) as never);
        const yj = (k: number) => top + ((v.x(0, k as never) - v.left) / w) * h;
        const x = Math.floor((x0 + fx * (x1 - x0)) * dpr);
        const y = Math.floor((yj(j) + fy * (yj(j + 1) - yj(j))) * dpr);
        const k = (y * main.width + x) * 4;
        return [px[k], px[k + 1], px[k + 2], px[k + 3]];
      };
      const blue = [0, 0, 255, 255];
      const red = [255, 0, 0, 255];
      expect(at(12, 27)).toEqual(blue);
      expect(at(19, 25)).toEqual(blue);
      expect(at(9, 27)).toEqual(red);
      expect(at(12, 30)).toEqual(red);
      expect(at(20, 29)).toEqual(red);
      if (cell) {
        // Zoomed in, a cell's corners are its own: the block's first cell and its neighbors, each its own color.
        expect(at(10, 25, 0.1, 0.1)).toEqual(blue);
        expect(at(10, 25, 0.9, 0.9)).toEqual(blue);
        expect(at(9, 25, 0.9, 0.5)).toEqual(red);
        expect(at(10, 24, 0.5, 0.9)).toEqual(red);
      }
    });

  it("the hover's crosshair: a column and a row through its cell on the overlay", async () => {
    t = await opened(doc(), { width: 900, height: 400, dpr: 1 });
    const v = t.ctx.layout.axis(t.section) as SectionAxisView;
    const l = t.ctx.layout.get(t.section)?.current;
    if (!l) throw new Error("no layout");
    const top = l.scrolledTop + l.entryTop[l.entryCount - 1];
    const h = l.entryHeight[l.entryCount - 1];
    const cx = (k: number) => (v.x(0, k as never) + v.x(0, (k + 1) as never)) / 2;
    const y = top + ((cx(30) - v.left) / (v.right - v.left)) * h;
    const hit = Hit.create();
    t.ctx.layout.hitTest(t.section, cx(5), y, hit);
    t.ctx.hover.setFromHit(hit, "pointer");
    settle(t);
    const overlay = pixelsOf(canvasesOf(t)[1]);
    const W = canvasesOf(t)[1].width;
    const alpha = (x: number, yy: number) => overlay[(Math.floor(yy) * W + Math.floor(x)) * 4 + 3];
    // Lit: position 5's column and 30's row (and 30's column, 5's row) stronger than the cells between them (which a
    // data hover's row band, the track's, may tint).
    const between = alpha(cx(17), top + ((cx(17) - v.left) / (v.right - v.left)) * h);
    expect(alpha(cx(5), top + 2)).toBeGreaterThan(between);
    expect(alpha(v.left + 2, y)).toBeGreaterThan(between);
  });
});
