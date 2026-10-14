/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The fingerprint as pixels in Chromium (seqspec:viewer#extension-points): the scene's rows through the Canvas2D
// backend. At a cell of a hydrogen bond the hydrogen bond's stripe is its kind's color and the stripes of the kinds the
// cell lacks are the ground, whatever the backend's color grouping; where a pose has no contact the cell is empty; the
// rows sit two device pixels apart. At DPR 1, 2 and 3.

import { afterEach, describe, expect, it } from "vitest";
import { Kind, KindStyle } from "../../../src/interactions/index.ts";
import { SectionLayouts } from "../../../src/layout.ts";
import { fingerprintFrames } from "./fingerprint-fixtures.ts";

let dispose: (() => void) | null = null;
afterEach(() => {
  dispose?.();
  dispose = null;
});

/** A CSS color as RGBA, through a canvas (`rgb(…)`, `#rrggbb`). */
function rgba(css: string): number[] {
  const c = document.createElement("canvas").getContext("2d") as CanvasRenderingContext2D;
  c.fillStyle = css;
  c.fillRect(0, 0, 1, 1);
  return Array.from(c.getImageData(0, 0, 1, 1).data);
}
const near = (a: number[], b: number[]) => a.every((x, i) => Math.abs(x - b[i]) <= 2);

describe("the unit scene's fingerprint on Canvas2D", () => {
  it.each([1, 2, 3])("stripes, grounds and empty cells (DPR %i)", async (dpr) => {
    const t = await fingerprintFrames("unit", undefined, { dpr, back: 0 });
    dispose = t.dispose;
    t.step();
    const g = t.canvas.getContext("2d") as CanvasRenderingContext2D;
    const s = t.ctx.graph.assembled.value?.view?.prepared[t.k].state as {
      kinds: Uint8Array;
      contacted: Int32Array[];
      colOf: Int32Array[];
      cells: Uint16Array[];
      width: Int32Array;
    };
    const l = t.ctx.layout.get(t.scene.section)?.current;
    const ax = t.ctx.layout.axis(t.scene.section);
    if (!l || !ax) throw new Error("no layout");
    const e = l.entryCount - 1;
    const top = l.scrolledTop + l.entryTop[e];
    // The row's content stands above its border (`SectionLayouts.entryBorder`).
    const border = SectionLayouts.entryBorder(l, e, t.ctx.layout.props(t.scene.section).rowBorders);
    const height = l.entryHeight[e] - (border ? Math.max(1, Math.round(dpr)) / dpr : 0);
    const a = Math.round((top + 1) * dpr);
    const H = Math.round((top + height - 1) * dpr) - a;
    const S = 4;
    const K = s.kinds.length;
    const G = H >= 6 * S ? 2 : 0;
    const row = (q: number) => {
      const y0 = a + Math.round((q * (H + G)) / S);
      return [y0, a + Math.round(((q + 1) * (H + G)) / S) - G];
    };
    const stripe = (q: number, j: number) => {
      const [y0, y1] = row(q);
      return [y0 + Math.round((j * (y1 - y0)) / K), y0 + Math.round(((j + 1) * (y1 - y0)) / K)];
    };
    const px = (index: number, y: number) => {
      const k = ax.shownOf(0, index);
      const x = Math.floor(((ax.x(0, k) + ax.x(0, (k + 1) as never)) / 2) * dpr);
      return Array.from(g.getImageData(x, y, 1, 1).data);
    };
    const ground = rgba(t.ctx.skin.value.grid);
    const hbond = rgba(KindStyle.hbond[1]);
    const at = (q: number, j: number, index: number) => {
      const [y0, y1] = stripe(q, j);
      const first = px(index, y0);
      expect(px(index, y1 - 1)).toEqual(first);
      return first;
    };
    // Residue 7: poses 1 and 2 hold the hydrogen bond (the top stripe), the other stripes are the ground.
    for (const q of [0, 1]) {
      expect(near(at(q, 0, 7), hbond)).toBe(true);
      for (let j = 1; j < K; j++) expect(near(at(q, j, 7), ground)).toBe(true);
    }
    // Pose 3 has no contact there: the cell is empty (the canvas was never cleared: transparent).
    for (let j = 0; j < K; j++) expect(at(2, j, 7)[3]).toBe(0);
    // The rows two device pixels apart: the gap above row 1 is no cell's.
    expect(G).toBe(2);
    expect(px(7, row(1)[0] - 1)[3]).toBe(0);
    expect(s.kinds[0]).toBe(Kind.hbond);
  });
});
