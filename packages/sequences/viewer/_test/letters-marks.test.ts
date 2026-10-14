/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A text-only letters row (no cell color: case 15's design rows) where text does not fit marks which positions hold a
// residue (seqspec:viewer#the-render-backend): a dot in each cell while a cell holds one, then a thin bar over each
// run of residues, so a deletion stays a break; a row with cell colors, or where text fits, marks nothing.

import { describe, expect, it } from "vitest";
import { SkinToken } from "../src/draw/backend.ts";
import { drawn, layerOf, record, rects, texts } from "./draw-fixtures.ts";
import { marksDoc } from "./letters-marks-fixtures.ts";

async function marks(color: boolean, cell: number, dpr: number, label = "Design") {
  const t = await drawn(marksDoc(color), { width: 600, height: 200, dpr });
  const { section } = layerOf(t.ctx, label);
  t.ctx.viewports.get(section)?.setView(cell, 0, "api");
  await t.ctx.graph.update();
  const { rec } = record(t.ctx, label);
  const g = rec.glyphItems;
  const dots = Array.from({ length: g.n }, (_, i) => ({ x: g.x[i], r: g.radius[i], c: g.color[i], a: g.alpha[i] }));
  const bars = rects(rec).filter((x) => x[4] === SkinToken.ink);
  const strings = texts(rec).map((d) => d.text);
  t.ctx.dispose();
  return { dots, bars, strings };
}

describe("a text-only letters row where text does not fit", () => {
  it("dots in the residues' cells while a cell holds a dot; none over the deletion", async () => {
    // 5 CSS px a cell at 2x: 10 device px, a dot of 4 (40%) across.
    const r = await marks(false, 5, 2);
    expect(r.strings).toEqual([]);
    expect(r.bars).toEqual([]);
    const step = r.dots[1].x - r.dots[0].x;
    expect(step).toBeGreaterThan(9);
    expect(r.dots.length).toBeGreaterThan(60);
    expect(r.dots.every((d) => d.c === SkinToken.ink && d.a === 128 && Math.abs(d.r - 0.2 * step) < 0.01)).toBe(true);
    // One cell apart, but across the deletion's ten.
    const gaps = r.dots.slice(1).map((d, i) => Math.round((d.x - r.dots[i].x) / step));
    expect(gaps.filter((x) => x !== 1)).toEqual([11]);
  });

  it("a bar per run of residues below that: the deletion a break", async () => {
    for (const [cell, dpr] of [
      [5, 1],
      [2, 2],
      [1.2, 1],
    ] as const) {
      const r = await marks(false, cell, dpr);
      expect(r.dots, `${cell} at ${dpr}x`).toEqual([]);
      expect(r.bars.length, `${cell} at ${dpr}x`).toBe(2);
      // Thin: a few device px, centered; the break as wide as the deletion's ten cells (30 are before it).
      expect(r.bars[0][3]).toBeLessThanOrEqual(3 * dpr);
      expect(r.bars[0][5]).toBe(128);
      const br = r.bars[1][0] - (r.bars[0][0] + r.bars[0][2]);
      expect(Math.abs(br - (r.bars[0][2] / 30) * 10), `${cell} at ${dpr}x`).toBeLessThanOrEqual(3);
    }
  });

  it("a row with cell colors, or where text fits, marks nothing; a plain row is one bar", async () => {
    for (const [color, cell, dpr] of [
      [true, 5, 2],
      [true, 1.2, 2],
      [false, 20, 1],
    ] as const) {
      const r = await marks(color, cell, dpr);
      expect(r.dots, `${color} ${cell}`).toEqual([]);
      expect(r.bars, `${color} ${cell}`).toEqual([]);
    }
    const plain = await marks(false, 5, 1, "Plain");
    expect(plain.dots).toEqual([]);
    expect(plain.bars.length).toBe(1);
    expect((await marks(false, 5, 2, "Plain")).dots.length).toBeGreaterThan(60);
  });
});
